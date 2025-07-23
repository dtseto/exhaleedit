/* exhaleEnc.h - header file for class providing Extended HE-AAC encoding capability
 * written by C. R. Helmrich, last modified in 2025 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2025 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _EXHALE_ENC_H_
#define _EXHALE_ENC_H_

#include "exhaleDecl.h"
#include "exhaleLibPch.h"
#include "bitAllocation.h"
#include "bitStreamWriter.h"
#include "entropyCoding.h"
#include "lappedTransform.h"
#include "linearPrediction.h"
#include "quantization.h"
#include "specAnalysis.h"
#include "specGapFilling.h"
#include "stereoProcessing.h"
#include "tempAnalysis.h"

// constant and experimental macro
#define WIN_SCALE double (1 << 23) //DS change to 24
#define EE_MORE_MSE              0 // DS change to 5 1-9: MSE optimized encoding with TNS disabled starting at bit-rate mode 1-9

// channelConfigurationIndex setup
typedef enum USAC_CCI : signed char
{
  CCI_UNDEF = -1,
  CCI_CONF  = 0,  // channel-to-speaker mapping defined in UsacChannelConfig() (not to be used here!)
  CCI_1_CH  = 1,  // 1.0: front-center
  CCI_2_CH  = 2,  // 2.0: front-left, front-right
  CCI_3_CH  = 3,  // 3.0: front-center, front-left, front-right
  CCI_4_CH  = 4,  // 4.0: front-center, front-left, front-right, back-center
  CCI_5_CH  = 5,  // 5.0: front-center, front-left, front-right, back-left, back-right
  CCI_6_CH  = 6,  // 5.1: front-center, front-left, front-right, back-left, back-right, LFE
  CCI_8_CH  = 7,  // 7.1: front-center, front-left, front-right, side-left, side-right, back-left, back-right, LFE
  CCI_2_CHM = 8,  // 2.0, dual-mono: channel1, channel2
  CCI_3_CHR = 9,  // 3.0, R-rotated: front-left, front-right, back-center
  CCI_4_CHR = 10, // 4.0, R-rotated: front-left, front-right, back-left, back-right
  CCI_7_CH  = 11, // 6.1: front-center, front-left, front-right, back-left, back-right, back-center, LFE
  CCI_8_CHS = 12  // 7.1, surround: front-center, front-L, front-R, surround-L, surround-R, back-L, back-R, LFE
} USAC_CCI;

// coreCoderFrameLength definition
typedef enum USAC_CCFL : short
{
  CCFL_UNDEF = -1,
#if !RESTRICT_TO_AAC
  CCFL_768   = 768, // LD
#endif
  CCFL_1024  = 1024 // LC
} USAC_CCFL;

// ==========================================================================================
// Static Constant Data (Moved from exhaleEnc.cpp)
// ==========================================================================================
// ISO/IEC 23003-3, Table 73
static const uint8_t numberOfChannels[USAC_MAX_NUM_ELCONFIGS] = {0, 1, 2, 3, 4, 5, 6, 8, 2, 3, 4, 7, 8};

// ISO/IEC 14496-3, Table 4.140
static const uint16_t sfbOffsetL0[42] = { // 88.2 and 96 kHz
    0,   4,   8,  12,  16,  20,  24,  28,  32,  36,  40,  44,  48,  52,  56,  64,  72,  80,  88,  96, 108,
  120, 132, 144, 156, 172, 188, 212, 240, 276, 320, 384, 448, 512, 576, 640, 704, 768, 832, 896, 960, 1024
};
// ISO/IEC 14496-3, Table 4.141
static const uint16_t sfbOffsetS0[13] = {
  0, 4, 8, 12, 16, 20, 24, 32, 40, 48, 64, 92, 128
};

// ISO/IEC 14496-3, Table 4.138
static const uint16_t sfbOffsetL1[48] = { // 64 kHz
    0,   4,   8,  12,  16,  20,  24,  28,  32,  36,  40,  44,  48,  52,  56,  64,  72,  80,  88, 100, 112, 124, 140, 156,
  172, 192, 216, 240, 268, 304, 344, 384, 424, 464, 504, 544, 584, 624, 664, 704, 744, 784, 824, 864, 904, 944, 984, 1024
};
// ISO/IEC 14496-3, Table 4.139
static const uint16_t sfbOffsetS1[13] = {
  0, 4, 8, 12, 16, 20, 24, 32, 40, 48, 64, 92, 128
};

// ISO/IEC 14496-3, Table 4.131
static const uint16_t sfbOffsetL2[52] = { // 32, 44.1, and 48 kHz
    0,   4,   8,  12,  16,  20,  24,  28,  32,  36,  40,  48,  56,  64,  72,  80,  88,  96, 108, 120, 132, 144, 160, 176, 196, 216, 240,
  264, 292, 320, 352, 384, 416, 448, 480, 512, 544, 576, 608, 640, 672, 704, 736, 768, 800, 832, 864, 896, 928, 960/*!*/, 992/*!*/, 1024
};
// ISO/IEC 14496-3, Table 4.130
static const uint16_t sfbOffsetS2[15] = {
  0, 4, 8, 12, 16, 20, 28, 36, 44, 56, 68, 80, 96, 112, 128
};

// ISO/IEC 14496-3, Table 4.136
static const uint16_t sfbOffsetL3[48] = { // 22.05 and 24 kHz
    0,   4,   8,  12,  16,  20,  24,  28,  32,  36,  40,  44,  52,  60,  68,  76,  84,  92, 100, 108, 116, 124, 136, 148,
  160, 172, 188, 204, 220, 240, 260, 284, 308, 336, 364, 396, 432, 468, 508, 552, 600, 652, 704, 768, 832, 896, 960, 1024
};
// ISO/IEC 14496-3, Table 4.137
static const uint16_t sfbOffsetS3[16] = {
  0, 4, 8, 12, 16, 20, 24, 28, 36, 44, 52, 64, 76, 92, 108, 128
};

// ISO/IEC 14496-3, Table 4.134
static const uint16_t sfbOffsetL4[44] = { // 11.025, 12, and 16 kHz
    0,   8,  16,  24,  32,  40,  48,  56,  64,  72,  80,  88, 100, 112, 124, 136, 148, 160, 172, 184, 196, 212,
  228, 244, 260, 280, 300, 320, 344, 368, 396, 424, 456, 492, 532, 572, 616, 664, 716, 772, 832, 896, 960, 1024
};
// ISO/IEC 14496-3, Table 4.135
static const uint16_t sfbOffsetS4[16] = {
  0, 4, 8, 12, 16, 20, 24, 28, 32, 40, 48, 60, 72, 88, 108, 128
};

// ISO/IEC 14496-3, Table 4.132
static const uint16_t sfbOffsetL5[41] = { // 8 kHz
    0,  12,  24,  36,  48,  60,  72,  84,  96, 108, 120, 132, 144, 156, 172, 188, 204, 220, 236, 252, 268,
  288, 308, 328, 348, 372, 396, 420, 448, 476, 508, 544, 580, 620, 664, 712, 764, 820, 880, 944, 1024
};
// ISO/IEC 14496-3, Table 4.133
static const uint16_t sfbOffsetS5[16] = {
  0, 4, 8, 12, 16, 20, 24, 28, 36, 44, 52, 60, 72, 88, 108, 128
};

// long-window SFB offset tables
static const uint16_t* swbOffsetsL[USAC_NUM_FREQ_TABLES] = {
  sfbOffsetL0, sfbOffsetL1, sfbOffsetL2, sfbOffsetL3, sfbOffsetL4, sfbOffsetL5
};
static const uint8_t numSwbOffsetL[USAC_NUM_FREQ_TABLES] = {42, 48, 52, 48, 44, 41};

// short-window SFB offset tables
static const uint16_t* swbOffsetsS[USAC_NUM_FREQ_TABLES] = {
  sfbOffsetS0, sfbOffsetS1, sfbOffsetS2, sfbOffsetS3, sfbOffsetS4, sfbOffsetS5
};
static const uint8_t numSwbOffsetS[USAC_NUM_FREQ_TABLES] = {13, 13, 15, 16, 16, 16};

// ISO/IEC 23003-3, Table 79
static const uint8_t freqIdxToSwbTableIdxAAC[USAC_NUM_SAMPLE_RATES + 2] = {
  /*96000*/ 0, 0, 1, 2, 2, 2,/*24000*/ 3, 3, 4, 4, 4, 5, 5, // AAC
  255, 255, 1, 2, 2, 2, 2, 2,/*25600*/ 3, 3, 3, 4, 4, 4, 4 // USAC
};
#if !RESTRICT_TO_AAC
static const uint8_t freqIdxToSwbTableIdx768[USAC_NUM_SAMPLE_RATES + 2] = {
  /*96000*/ 0, 0, 0, 1, 1, 2,/*24000*/ 2, 2, 3, 4, 4, 4, 4, // AAC
  255, 255, 0, 1, 2, 2, 2, 2,/*25600*/ 2, 3, 3, 3, 3, 4, 4 // USAC
};
#endif

// ISO/IEC 23003-3, Table 131
static const uint8_t tnsScaleFactorBandLimit[2 /*long/short*/][USAC_NUM_FREQ_TABLES] = { // TNS_MAX_BANDS
  {31, 34, 51 /*to be corrected to 42 (44.1) and 40 (48 kHz)!*/, 47, 43, 40}, {9, 10, 14, 15, 15, 15}
};

static const uint8_t sbrRateOffset[10] = {7, 6, 6, 8, 7, 8, 9, 9, 9, 9}; // used for scaleSBR

// scale_factor_grouping map
// group lengths based on transient location:  1133, 1115, 2114, 3113, 4112, 5111, 3311, 1331
static const uint8_t scaleFactorGrouping[8] = {0x1B, 0x0F, 0x47, 0x63, 0x71, 0x78, 0x6C, 0x36};

static const uint8_t windowGroupingTable[8][NUM_WINDOW_GROUPS] = { // for window_group_length
  {1, 1, 3, 3}, {1, 1, 1, 5}, {2, 1, 1, 4}, {3, 1, 1, 3}, {4, 1, 1, 2}, {5, 1, 1, 1}, {3, 3, 1, 1}, {1, 3, 3, 1}
};

// window_sequence equalizer
static const USAC_WSEQ windowSequenceSynch[5][5] = {  // 1st: chan index 0, 2nd: chan index 1
  {ONLY_LONG,   LONG_START,  EIGHT_SHORT, LONG_STOP,   STOP_START }, // left: ONLY_LONG
#if RESTRICT_TO_AAC
  {LONG_START,  LONG_START,  EIGHT_SHORT, EIGHT_SHORT, STOP_START }, // Left: LONG_START
#else
  {LONG_START,  LONG_START,  EIGHT_SHORT, STOP_START,  STOP_START }, // Left: LONG_START
#endif
  {EIGHT_SHORT, EIGHT_SHORT, EIGHT_SHORT, EIGHT_SHORT, EIGHT_SHORT}, // Left: EIGHT_SHORT
#if RESTRICT_TO_AAC
  {LONG_STOP,   EIGHT_SHORT, EIGHT_SHORT, LONG_STOP,   STOP_START }, // Left: LONG_STOP
#else
  {LONG_STOP,   STOP_START,  EIGHT_SHORT, LONG_STOP,   STOP_START }, // Left: LONG_STOP
#endif
  {STOP_START,  STOP_START,  EIGHT_SHORT, STOP_START,  STOP_START }  // Left: STOP_START
};

// overall BL USAC encoding class
class ExhaleEncoder : public ExhaleEncAPI
{
private:

  // member variables
  uint16_t        m_bandwidCurr[USAC_MAX_NUM_CHANNELS];
  uint16_t        m_bandwidPrev[USAC_MAX_NUM_CHANNELS];
  BitAllocator    m_bitAllocator; // for scale factor init
  uint8_t         m_bitRateMode;
  USAC_CCI        m_channelConf;
  int32_t*        m_coreSignals[USAC_MAX_NUM_CHANNELS];
  CoreCoderData*  m_elementData[USAC_MAX_NUM_ELEMENTS];
  EntropyCoder    m_entropyCoder[USAC_MAX_NUM_CHANNELS];
  uint32_t        m_frameCount;
  USAC_CCFL       m_frameLength;
  int8_t          m_frequencyIdx;
  bool            m_indepFlag; // usacIndependencyFlag bit
  uint32_t        m_indepPeriod;
  LinearPredictor m_linPredictor; // for pre-roll est, TNS
  uint8_t*        m_mdctQuantMag[USAC_MAX_NUM_CHANNELS];
  int32_t*        m_mdctSignals[USAC_MAX_NUM_CHANNELS];
  int32_t*        m_mdstSignals[USAC_MAX_NUM_CHANNELS];
  uint8_t         m_meanSpecPrev[USAC_MAX_NUM_CHANNELS]; // for
  uint8_t         m_meanTempPrev[USAC_MAX_NUM_CHANNELS]; // SBR
#if !RESTRICT_TO_AAC
  bool            m_noiseFilling[USAC_MAX_NUM_ELEMENTS];
#endif
  bool            m_nonMpegExt;
  uint8_t         m_numElements;
  uint8_t         m_numSwbLong;
  uint8_t         m_numSwbShort;
  unsigned char*  m_outAuData;
  BitStreamWriter m_outStream; // for access unit creation
  int32_t*        m_pcm24Data;
  uint8_t         m_perCorrHCurr[USAC_MAX_NUM_ELEMENTS];
  uint8_t         m_perCorrLCurr[USAC_MAX_NUM_ELEMENTS];
  uint8_t         m_priLength;
  uint32_t        m_rateFactor; // RC
  SfbGroupData*   m_scaleFacData[USAC_MAX_NUM_CHANNELS];
  uint16_t        m_sfbLoudMem[2][26][32]; // loudness mem
  SfbQuantizer    m_sfbQuantizer; // powerlaw quantization
  uint8_t         m_shiftValSBR; // SBR ratio for shifting
  SpecAnalyzer    m_specAnalyzer; // for spectral analysis
  uint32_t        m_specAnaCurr[USAC_MAX_NUM_CHANNELS];
  uint8_t         m_specFlatPrev[USAC_MAX_NUM_CHANNELS];
#if !RESTRICT_TO_AAC
  SpecGapFiller   m_specGapFiller;// for noise/gap filling
#endif
  StereoProcessor m_stereoCoder;  // for M/S stereo coding
  uint8_t         m_swbTableIdx;
  TempAnalyzer    m_tempAnalyzer; // for temporal analysis
  uint32_t        m_tempAnaCurr[USAC_MAX_NUM_CHANNELS];
  uint32_t        m_tempAnaNext[USAC_MAX_NUM_CHANNELS];
  uint8_t         m_tempFlatPrev[USAC_MAX_NUM_CHANNELS];
  int32_t*        m_tempIntBuf;  // temporary int32 buffer
  int32_t*        m_timeSignals[USAC_MAX_NUM_CHANNELS];
#if !RESTRICT_TO_AAC
  uint8_t         m_timeWarpTCX[USAC_MAX_NUM_ELEMENTS]; // for TW, TCX
#endif
  int32_t*        m_timeWindowL[2];  // long window halves
  int32_t*        m_timeWindowS[2]; // short window halves
  int16_t         m_tranLocCurr[USAC_MAX_NUM_CHANNELS];
  int16_t         m_tranLocNext[USAC_MAX_NUM_CHANNELS];
  LappedTransform m_transform; // time-frequency transform

    unsigned        m_targetBandwidth; // ADD THIS LINE

  // helper functions
  unsigned applyTnsToWinGroup (SfbGroupData& grpData, const uint8_t grpIndex, const uint8_t maxSfb, TnsData& tnsData,
                               const unsigned channelIndex, const unsigned n, const bool realOnlyCalc);
  unsigned eightShortGrouping (SfbGroupData& grpData, uint16_t* const grpOffsets,
                               int32_t* const mdctSignal, int32_t* const mdstSignal);
  unsigned getOptParCorCoeffs (const SfbGroupData& grpData, const uint8_t maxSfb, TnsData& tnsData,
                               const unsigned channelIndex, const uint8_t firstGroupIndexToTest = 0);
  uint32_t getThr             (const unsigned channelIndex, const unsigned sfbIndex);
  unsigned psychBitAllocation ();
  unsigned quantizationCoding ();
  unsigned spectralProcessing ();
  unsigned temporalProcessing ();

public:

  // constructor
  ExhaleEncoder (int32_t* const inputPcmData,       unsigned char* const outputAuData,
                 const unsigned sampleRate = 44100, const unsigned numChannels = 2,
                 const unsigned frameLength = 1024, const unsigned indepPeriod = 45,
                 const unsigned varBitRateMode = 3
#if !RESTRICT_TO_AAC
               , const bool useNoiseFilling = true, const bool useEcodisExt = false
#endif
    );
  // destructor
  virtual ~ExhaleEncoder ();
  // public functions
  unsigned encodeLookahead ();
  unsigned encodeFrame ();
  unsigned initEncoder (unsigned char* const audioConfigBuffer, uint32_t* const audioConfigBytes = nullptr);

}; // ExhaleEncoder

#endif // _EXHALE_ENC_H_
