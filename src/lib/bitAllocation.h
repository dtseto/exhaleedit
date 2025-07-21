/* bitAllocation.h - header file for class needed for psychoacoustic bit-allocation
 * written by C. R. Helmrich, last modified in 2025 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2025 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _BIT_ALLOCATION_H_
#define _BIT_ALLOCATION_H_

#include "exhaleLibPch.h"
#include "linearPrediction.h"

// constants, experimental macros
#define BA_EPS                  1
#define BA_MORE_CBR            0 // DS change 1 1: force more constant bit-rate (CBR, experimental!)

// class for audio bit-allocation
class BitAllocator
{
private:
    
    // member variables
    uint32_t m_avgStepSize[USAC_MAX_NUM_CHANNELS];
    uint8_t  m_avgSpecFlat[USAC_MAX_NUM_CHANNELS];
    uint8_t  m_avgTempFlat[USAC_MAX_NUM_CHANNELS];
    uint8_t  m_rateIndex; // preset
    uint8_t* m_tempSfbValue;
    LinearPredictor* m_tnsPredictor;
    
    // ADD: New member variable for tonality
    float    m_tonalityIndex[USAC_MAX_NUM_CHANNELS][MAX_NUM_SWB_LONG];
    
    // ADD: New member variable for Minimum SNR model
    double   m_sfbMinSnr[USAC_MAX_NUM_CHANNELS][MAX_NUM_SWB_LONG];

    // ADD: New private method for Minimum SNR calculation
    void initMinSnr (const unsigned nChannels, const unsigned samplingRate, const long bitrate, const SfbGroupData* const groupData[USAC_MAX_NUM_CHANNELS]);

public:
    
    // constructor
    BitAllocator ();
    // destructor
    ~BitAllocator () { MFREE (m_tempSfbValue); }
    // public functions
    void getChAverageSpecFlat (uint8_t meanSpecFlatInCh[USAC_MAX_NUM_CHANNELS], const unsigned nChannels);
    void getChAverageTempFlat (uint8_t meanTempFlatInCh[USAC_MAX_NUM_CHANNELS], const unsigned nChannels);
    uint16_t   getRateCtrlFac (const int32_t rateRatio, const unsigned samplingRate, const uint32_t specFlatness,
                               const bool prevEightShorts = false);
    uint8_t       getScaleFac (const uint32_t sfbStepSize, const int32_t* const sfbSignal, const uint8_t sfbWidth,
                               const uint32_t sfbRmsValue);
    unsigned initAllocMemory  (LinearPredictor* const linPredictor, const uint8_t numSwb, const uint8_t bitRateMode);
    
    // MODIFY: Add bitrate parameter for Min SNR calculation
    unsigned initSfbStepSizes (const SfbGroupData* const groupData[USAC_MAX_NUM_CHANNELS], const uint8_t numSwbShort,
                               const uint32_t specAnaStats[USAC_MAX_NUM_CHANNELS],
                               const uint32_t tempAnaStats[USAC_MAX_NUM_CHANNELS],
                               const unsigned nChannels, const unsigned samplingRate, uint32_t* const sfbStepSizes,
                               const unsigned lfeChannelIndex, const unsigned ad = 0u, const bool tnsDisabled = false, const long bitrate = 64000);
    unsigned imprSfbStepSizes (const SfbGroupData* const groupData[USAC_MAX_NUM_CHANNELS], const uint8_t numSwbShort,
                               const int32_t* const mdctSpec[USAC_MAX_NUM_CHANNELS], const unsigned nSamplesInFrame,
                               const unsigned nChannels, const unsigned samplingRate, uint32_t* const sfbStepSizes,
                               const unsigned firstChannelIndex, const uint8_t* const sfm, const bool commonWindow,
                               const uint8_t* const sfbStereoData = nullptr, const uint8_t stereoConfig = 0);



    // ADD: New method declaration
    float estimateTonality(const int32_t* mdctSpec, const uint16_t* sfbOffsets,
                          const unsigned sfb, const unsigned sfbWidth);

    
}; // BitAllocator

#endif // _BIT_ALLOCATION_H_
