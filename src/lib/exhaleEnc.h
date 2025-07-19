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
#define WIN_SCALE double (1 << 24) //DS change to 24
#define EE_MORE_MSE              5 // DS change to 5 1-9: MSE optimized encoding with TNS disabled starting at bit-rate mode 1-9

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
