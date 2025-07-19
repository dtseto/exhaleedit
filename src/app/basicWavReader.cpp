/* basicWavReader.cpp - source file for class with basic WAVE file reading capability
 * written by C. R. Helmrich, last modified in 2024 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2024 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleAppPch.h"
#include "basicWavReader.h"

// static helper functions
static unsigned reverseFourBytes (const uint8_t* b)
{
  return ((unsigned) b[3] << 24) | ((unsigned) b[2] << 16) | ((unsigned) b[1] << 8) | (unsigned) b[0];
}

static int64_t fourBytesToLength (const uint8_t* b, const int64_t lengthLimit)
{
  int64_t chunkLength = (int64_t) reverseFourBytes (b);

  chunkLength += chunkLength & 1;  // make sure it is even

  return __min (lengthLimit, chunkLength); // for security
}

static inline unsigned getFrames (const int fileHandle, void* dataBuf, const unsigned frameCount,
                                  const bool roundSize, const unsigned bytesPerFrame)
{
  const int size = ((frameCount + (roundSize ? 1 << (BWR_READ_FRACT - 1) : 0)) >> BWR_READ_FRACT) * bytesPerFrame;
  int  bytesRead = _READ (fileHandle, dataBuf, size); // 1

  if ((bytesRead = __max (0, bytesRead)) < size)
  {
    const int br = _READ (fileHandle, (uint8_t*) dataBuf + bytesRead, size - bytesRead);

    bytesRead += __max (0, br); // bytes of read attempt 2
  }
  return bytesRead / bytesPerFrame;  // num of frames read
}

static const uint8_t allowedChMask[8] = {0x00, 0x04, 0xC0, 0x00, 0x00, 0x37, 0xCF, 0xFF};

// private reader functions
bool BasicWavReader::readRiffHeader ()
{
  uint8_t b[FILE_HEADER_SIZE] = {0};  // temp. byte buffer

  m_bytesRead = _READ (m_fileHandle, b, 8);
  if ((m_bytesRead += _READ (m_fileHandle, &b[8], FILE_HEADER_SIZE - 8)) != FILE_HEADER_SIZE) return false; // error
  m_bytesRemaining -= m_bytesRead;
  m_chunkLength = fourBytesToLength (&b[4], m_bytesRemaining) - 4; // minus 4 bytes for WAVE tag

  return (b[0] == 'R' && b[1] == 'I' && b[2] == 'F' && b[3] == 'F' &&
          b[8] == 'W' && b[9] == 'A' && b[10]== 'V' && b[11]== 'E' &&
          m_bytesRemaining > 32);  // true: RIFF supported
}

bool BasicWavReader::readFormatChunk ()
{
  uint8_t b[CHUNK_FORMAT_MAX] = {0};  // temp. byte buffer

  if (!seekToChunkTag (b, 0x20746D66 /*fmt */) || (m_chunkLength < CHUNK_FORMAT_SIZE) || (m_chunkLength > CHUNK_FORMAT_MAX))
  {
    return false; // fmt_ chunk invalid or read incomplete
  }
  m_bytesRead = _READ (m_fileHandle, b, 8);
  for (int64_t i = 8; i < m_chunkLength; i += 2) m_bytesRead += _READ (m_fileHandle, &b[i], 2);
  if (m_bytesRead != m_chunkLength) return false; // error
  m_bytesRemaining -= m_bytesRead;

  m_waveChannels  = (unsigned (b[3]) << 8) | b[2]; // <64k
  if ((b[0] == 0xFE) && (b[1] == 0xFF) && (m_chunkLength == CHUNK_FORMAT_MAX) && (b[16] == CHUNK_FORMAT_MAX - CHUNK_FORMAT_SIZE - 2) &&
      (b[17] == 0) && (b[18] <= b[14]) && ((b[19] | b[25] | b[26] | b[27] | b[28] | b[29] | b[31] | b[33] | b[34] | b[36]) == 0))
  {
    m_waveDataType = WAV_TYPE (b[24]-1); // extensible WAV
    if (b[21] == 6) b[20] |= 3 << (m_waveChannels - 2);
    if ((b[18] + 8u <= b[14]) || (b[20] != 0 && b[20] != allowedChMask[__min (7, m_waveChannels)] &&
        b[20] != (1u << __min (8, m_waveChannels)) - 1u)) return false;
    b[14] = b[18];
    b[ 1] = 0;
  }
  else
  m_waveDataType  = WAV_TYPE (b[0]-1); // 1: PCM, 3: float
  m_waveFrameRate = reverseFourBytes (&b[4]);  // frames/s
  m_waveBitRate   = reverseFourBytes (&b[8]) * 8; // bit/s
  m_waveFrameSize = b[12];  // bytes/s divided by frames/s
  m_waveBitDepth  = b[14]; // only 8, 16, 24, 32 supported

  return ((m_waveDataType == WAV_PCM || (m_waveDataType == WAV_FLOAT && (m_waveBitDepth & 15) == 0)) &&
          (m_waveChannels > 0 && m_waveChannels <= 63) && isSamplingRateSupported (m_waveFrameRate) &&
          (m_waveBitRate == 8 * m_waveFrameRate * m_waveFrameSize) && (b[ 1] == 0) && (b[ 3] == 0) &&
          (m_waveFrameSize * 8 == m_waveBitDepth * m_waveChannels) && (b[13] == 0) && (b[15] == 0) &&
          (m_waveBitDepth >= 8 && m_waveBitDepth <= 32 && (m_waveBitDepth & 7) == 0) &&
          m_bytesRemaining > 8); // true: format supported
}

bool BasicWavReader::readDataHeader ()
{
  uint8_t b[CHUNK_HEADER_SIZE] = {0}; // temp. byte buffer

  if (!seekToChunkTag (b, 0x61746164 /*data*/))
  {
    return false; // data chunk invalid or read incomplete
  }
  return (m_chunkLength > 0); // true: WAVE data available
}

// private helper function
bool BasicWavReader::seekToChunkTag (uint8_t* const buf, const uint32_t tagID)
{
  if ((m_bytesRead = _READ (m_fileHandle, buf, CHUNK_HEADER_SIZE)) != CHUNK_HEADER_SIZE) return false; // error
  m_bytesRemaining -= m_bytesRead;
  m_chunkLength = fourBytesToLength (&buf[4], m_bytesRemaining);

  while ((*((uint32_t* const) buf) != tagID) &&
         (m_bytesRemaining > 0)) // seek until tagID found
  {
    if ((m_bytesRemaining > LLONG_MAX - USHRT_MAX) || (m_readOffset = _SEEK (m_fileHandle, m_chunkLength, 1 /*SEEK_CUR*/)) == -1)
    {
      // for stdin compatibility, don't abort, try reading
      for (int64_t i = m_chunkLength >> 1; i > 0; i--) m_bytesRead = _READ (m_fileHandle, buf, 2);
    }
    m_bytesRemaining -= m_chunkLength;

    if (m_bytesRemaining <= 0) return false;  // unlikely!

    if ((m_bytesRead = _READ (m_fileHandle, buf, CHUNK_HEADER_SIZE)) != CHUNK_HEADER_SIZE) return false; // error
    m_bytesRemaining -= m_bytesRead;
    m_chunkLength = fourBytesToLength (&buf[4], m_bytesRemaining);
  }
  return (m_bytesRemaining > 0);
}

// static reading functions
unsigned BasicWavReader::readDataFloat16 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                          const unsigned chanCount, void* tempBuf)
{
  const unsigned rest = ((frameCount >> (BWR_READ_FRACT - 1)) << (BWR_READ_FRACT - 1)) < frameCount ? 1 : 0;
  unsigned framesRead = 0;

  for (unsigned fract = 0; fract < (1 << BWR_READ_FRACT) + rest; fract++)
  {
    const int16_t* fBuf = (const int16_t*) tempBuf; // words
    const unsigned read = getFrames (fileHandle, tempBuf, frameCount, (fract & 1), chanCount * 2);

    for (unsigned i = read * chanCount; i > 0; i--)
    {
      const int16_t i16 = *(fBuf++);
      const int32_t e = ((i16 & 0x7C00) >> 10) - 18; // exp.
      // an exponent e <= -12 will lead to zero-quantization
      *frameBuf = int32_t (e < 0 ? (1024 + (i16 & 0x03FF) + (1 << (-1 - e)) /*rounding offset*/) >> -e
                                 : (e > 12 ? MAX_VALUE_AUDIO24 /*inf*/ : (1024 + (i16 & 0x03FF)) << e));
      if ((i16 & 0x8000) != 0) *frameBuf *= -1; // neg. sign
      frameBuf++;
    }
    framesRead += read;
  }
  if (framesRead < frameCount) // zero out missing samples
  {
    memset (frameBuf, 0, (frameCount - framesRead) * chanCount * sizeof (int32_t));
  }
  return framesRead;
}

unsigned BasicWavReader::readDataFloat32 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                          const unsigned chanCount, void* tempBuf)
{
  const unsigned rest = ((frameCount >> (BWR_READ_FRACT - 1)) << (BWR_READ_FRACT - 1)) < frameCount ? 1 : 0;
  unsigned framesRead = 0;

  for (unsigned fract = 0; fract < (1 << BWR_READ_FRACT) + rest; fract++)
  {
    const float*   fBuf = (const float*) tempBuf; // 4 bytes
    const unsigned read = getFrames (fileHandle, tempBuf, frameCount, (fract & 1), chanCount * 4);

    for (unsigned i = read * chanCount; i > 0; i--)
    {
      const float   f32 = *fBuf * float (1 << 23); // * 2^23
      fBuf++;
      *frameBuf = int32_t (f32 + (f32 < 0.0 ? -0.5 : 0.5)); // rounding
      if (*frameBuf < MIN_VALUE_AUDIO24) *frameBuf = MIN_VALUE_AUDIO24;
      else
      if (*frameBuf > MAX_VALUE_AUDIO24) *frameBuf = MAX_VALUE_AUDIO24;
      frameBuf++;
    }
    framesRead += read;
  }
  if (framesRead < frameCount) // zero out missing samples
  {
    memset (frameBuf, 0, (frameCount - framesRead) * chanCount * sizeof (int32_t));
  }
  return framesRead;
}

unsigned BasicWavReader::readDataLnPcm08 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                          const unsigned chanCount, void* tempBuf)
{
  const unsigned rest = ((frameCount >> (BWR_READ_FRACT - 1)) << (BWR_READ_FRACT - 1)) < frameCount ? 1 : 0;
  unsigned framesRead = 0;

  for (unsigned fract = 0; fract < (1 << BWR_READ_FRACT) + rest; fract++)
  {
    const uint8_t* iBuf = (const uint8_t*) tempBuf; // 1b
    const unsigned read = getFrames (fileHandle, tempBuf, frameCount, (fract & 1), chanCount);

    for (unsigned i = read * chanCount; i > 0; i--)
    {
      *(frameBuf++) = ((int32_t) *(iBuf++) - 128) << 16; // * 2^16
    }
    framesRead += read;
  }
  if (framesRead < frameCount) // zero out missing samples
  {
    memset (frameBuf, 0, (frameCount - framesRead) * chanCount * sizeof (int32_t));
  }
  return framesRead;
}

unsigned BasicWavReader::readDataLnPcm16 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                          const unsigned chanCount, void* tempBuf)
{
  const unsigned rest = ((frameCount >> (BWR_READ_FRACT - 1)) << (BWR_READ_FRACT - 1)) < frameCount ? 1 : 0;
  unsigned framesRead = 0;

  for (unsigned fract = 0; fract < (1 << BWR_READ_FRACT) + rest; fract++)
  {
    const int16_t* iBuf = (const int16_t*) tempBuf; // words
    const unsigned read = getFrames (fileHandle, tempBuf, frameCount, (fract & 1), chanCount * 2);

    for (unsigned i = read * chanCount; i > 0; i--)
    {
      *(frameBuf++) = (int32_t) *(iBuf++) * (1 << 8); // * 2^8
    }
    framesRead += read;
  }
  if (framesRead < frameCount) // zero out missing samples
  {
    memset (frameBuf, 0, (frameCount - framesRead) * chanCount * sizeof (int32_t));
  }
  return framesRead;
}

unsigned BasicWavReader::readDataLnPcm24 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                          const unsigned chanCount, void* tempBuf)
{
  const unsigned rest = ((frameCount >> (BWR_READ_FRACT - 1)) << (BWR_READ_FRACT - 1)) < frameCount ? 1 : 0;
  unsigned framesRead = 0;

  for (unsigned fract = 0; fract < (1 << BWR_READ_FRACT) + rest; fract++)
  {
    const uint8_t* iBuf = (const uint8_t*) tempBuf; // 3b
    const unsigned read = getFrames (fileHandle, tempBuf, frameCount, (fract & 1), chanCount * 3);

    for (unsigned i = read * chanCount; i > 0; i--)
    {
      const int32_t i24 = (int32_t) iBuf[0] | ((int32_t) iBuf[1] << 8) | ((int32_t) iBuf[2] << 16);
      iBuf += 3;
      *(frameBuf++) = (i24 > MAX_VALUE_AUDIO24 ? i24 + 2 * MIN_VALUE_AUDIO24 : i24);
    }
    framesRead += read;
  }
  if (framesRead < frameCount) // zero out missing samples
  {
    memset (frameBuf, 0, (frameCount - framesRead) * chanCount * sizeof (int32_t));
  }
  return framesRead;
}

unsigned BasicWavReader::readDataLnPcm32 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                          const unsigned chanCount, void* tempBuf)
{
  const unsigned rest = ((frameCount >> (BWR_READ_FRACT - 1)) << (BWR_READ_FRACT - 1)) < frameCount ? 1 : 0;
  unsigned framesRead = 0;

  for (unsigned fract = 0; fract < (1 << BWR_READ_FRACT) + rest; fract++)
  {
    const int32_t* iBuf = (const int32_t*) tempBuf; // dword
    const unsigned read = getFrames (fileHandle, tempBuf, frameCount, (fract & 1), chanCount * 4);

    for (unsigned i = read * chanCount; i > 0; i--)
    {
      const int32_t i24 = ((*iBuf >> 1) + (1 << 6)) >> 7; // * 2^-8 with rounding, overflow-safe
      iBuf++;
      *(frameBuf++) = __min (MAX_VALUE_AUDIO24, i24);
    }
    framesRead += read;
  }
  if (framesRead < frameCount) // zero out missing samples
  {
    memset (frameBuf, 0, (frameCount - framesRead) * chanCount * sizeof (int32_t));
  }
  return framesRead;
}

// public functions
unsigned BasicWavReader::open (const int wavFileHandle, const uint16_t maxFrameRead, const int64_t fileLength /*= LLONG_MAX*/)
{
  m_bytesRemaining = fileLength;
  m_fileHandle     = wavFileHandle;

  if ((m_fileHandle == -1) || (fileLength <= 44))
  {
    return 1; // file handle invalid or file too small
  }
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
  if ((fileLength < LLONG_MAX) && (m_readOffset = _telli64 (m_fileHandle)) != 0)
#else // Linux, MacOS, Unix
  if ((fileLength < LLONG_MAX) && (m_readOffset = lseek (m_fileHandle, 0, 1 /*SEEK_CUR*/)) != 0)
#endif
  {
    m_readOffset = _SEEK (m_fileHandle, 0, 0 /*SEEK_SET*/);
  }
  if ((m_readOffset != 0) || !readRiffHeader ())
  {
    return 2; // file type invalid or file seek failed
  }
  if (!readFormatChunk ())
  {
    return 3; // audio format invalid or not supported
  }
  if (!readDataHeader ())
  {
    return 4; // WAVE data part invalid or unsupported
  }
  if ((m_byteBuffer = (char*) malloc (m_waveFrameSize * ((maxFrameRead + (1 << (BWR_READ_FRACT - 1))) >> BWR_READ_FRACT))) == nullptr)
  {
    return 5; // read-in byte buffer allocation failed
  }
  m_frameLimit = maxFrameRead;

  if (m_waveChMpegMap > 0)
  {
    if (m_waveChannels < 3) m_waveChMpegMap = 0; // mono and stereo inputs do not need to be remapped
    else
    if (m_waveChannels == 6) m_waveChMpegMap++; // when m_waveChMpegMap is even, relocate LFE channel
  }

  // ready to read audio data: initialize byte counter
  if (m_bytesRemaining > m_chunkLength)
  {
    m_bytesRemaining = m_chunkLength;
  }
  m_chunkLength = 0;

  if (m_waveDataType == WAV_PCM) // & function pointer
  {
    switch (m_waveBitDepth)
    {
      case 8:
        m_readDataFunc = readDataLnPcm08; break;
      case 16:
        m_readDataFunc = readDataLnPcm16; break;
      case 24:
        m_readDataFunc = readDataLnPcm24; break;
      default:
        m_readDataFunc = readDataLnPcm32; break;
    }
  }
  else  m_readDataFunc = (m_waveBitDepth == 16 ? readDataFloat16 : readDataFloat32);

  return (m_readDataFunc == nullptr ? 6 : 0); // 0: OK
}

unsigned BasicWavReader::read (int32_t* const frameBuf, const uint16_t frameCount)
{
  const unsigned framesTotal = __min (m_frameLimit, frameCount);
  unsigned framesRead;

  if ((frameBuf == nullptr) || (m_fileHandle == -1) || (framesTotal == 0) || (m_byteBuffer == nullptr) ||
      (m_bytesRemaining <= 0)) // end of chunk reached
  {
    if (frameBuf != nullptr) memset (frameBuf, 0, framesTotal * m_waveChannels * sizeof (int32_t));

    return 0; // invalid args or class not initialized
  }
  framesRead  = m_readDataFunc (m_fileHandle, frameBuf, framesTotal, m_waveChannels, m_byteBuffer);
  m_bytesRead = m_waveFrameSize * framesRead;
  if ((m_bytesRemaining -= m_bytesRead) < 0)
  {
    m_bytesRead = unsigned (m_bytesRead + m_bytesRemaining);
    framesRead  = m_bytesRead / m_waveFrameSize;
  }
  m_chunkLength += m_bytesRead;

  if (framesRead < framesTotal) eaExtrapolate (frameBuf, framesRead, framesTotal, m_waveChannels); // fade-out, for gapless playback on more content

  if (m_waveChMpegMap > 0)
  {
    int32_t* sampBuf = frameBuf; // remap multichannel PCM input to MPEG-4 or D channel configuration

    if (m_waveChMpegMap & 1)
    {
      for (unsigned i = framesTotal; i > 0; i--, sampBuf += m_waveChannels)
      {
        int32_t  c = sampBuf[2];
        sampBuf[2] = sampBuf[1];
        sampBuf[1] = sampBuf[0];
        sampBuf[0] = c; // cntr.
      }
    }
    else // m_waveChMpegMap even
    {
      for (unsigned i = framesTotal; i > 0; i--, sampBuf += m_waveChannels)
      {
        int32_t  c = sampBuf[2], /*int32*/ l = sampBuf[3];
        sampBuf[2] = sampBuf[1];  sampBuf[3] = sampBuf[4];
        sampBuf[1] = sampBuf[0];  sampBuf[4] = sampBuf[5];
        sampBuf[0] = c; /*cntr.*/ sampBuf[5] = l; // LFE
      }
    }
  }

  return framesRead;
}

void BasicWavReader::reset ()
{
  m_byteBuffer     = nullptr;
  m_bytesRead      = 0;
  m_bytesRemaining = 0;
  m_chunkLength    = 0;
  m_frameLimit     = 0;
  m_readDataFunc   = nullptr;
  m_readOffset     = 0;
  m_waveBitDepth   = 0;
  m_waveChannels   = 0;
  m_waveFrameRate  = 0;

  if (m_fileHandle != -1) _SEEK (m_fileHandle, 0, 0 /*SEEK_SET*/);
}
