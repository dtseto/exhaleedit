/* bitAllocation.cpp - source file for class needed for psychoacoustic bit-allocation
 * written by C. R. Helmrich, last modified in 2025 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2025 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "bitAllocation.h"

// static helper functions
static inline uint32_t intSqrt (const uint32_t val)
{
  return uint32_t (0.5 + sqrt ((double) val));
}

static inline uint32_t jndModel (const uint32_t val, const uint32_t mean,
                                 const unsigned expTimes512, const unsigned mulTimes512)
{
  const double exp = (double) expTimes512 / 512.0; // exponent
  const double mul = (double) mulTimes512 / 512.0; // a factor
  const double res = pow (mul * (double) val, exp) * pow ((double) mean, 1.0 - exp);

  return uint32_t (__min ((double) UINT_MAX, res + 0.5));
}

static inline uint32_t squareMeanRoot (const uint32_t value1, const uint32_t value2)
{
  const double meanRoot = (sqrt ((double) value1) + sqrt ((double) value2)) * 0.5;

  return uint32_t (meanRoot * meanRoot + 0.5);
}

static void jndPowerLawAndPeakSmoothing (uint32_t* const  stepSizes, const unsigned nStepSizes, const bool lowRateMode,
                                         const uint32_t avgStepSize, const uint8_t sfm, const uint8_t tfm)
{
  const unsigned  expTimes512 = 512u - sfm; // 1.0 - sfm / 2.0
  const unsigned  mulTimes512 = __min (expTimes512, 512u - tfm);
  uint32_t         stepSizeM3 = 0, stepSizeM2 = 0, stepSizeM1 = BA_EPS; // hearing threshold around zero Hz
  unsigned b;

  for (b = 0; b < __min (2, nStepSizes); b++)
  {
    stepSizeM3 = stepSizeM2;
    stepSizeM2 = stepSizeM1;
    stepSizeM1 = stepSizes[b] = jndModel (stepSizes[b], avgStepSize, expTimes512, mulTimes512);
  }
  stepSizes[0] = __min (stepSizeM1, stepSizes[0]); // `- becomes --
  for (/*b*/; b < nStepSizes; b++)
  {
    const uint64_t modifiedB = (lowRateMode ? 16 + b : b);
    const uint64_t oneMinusB = 128 - modifiedB;
    const uint32_t stepSizeB = jndModel (stepSizes[b], avgStepSize, expTimes512, mulTimes512);

    if ((stepSizeM3 <= stepSizeM2) && (stepSizeM3 <= stepSizeM1) && (stepSizeB <= stepSizeM2) && (stepSizeB <= stepSizeM1))
    {
      const uint32_t maxM3M0 = __max (stepSizeM3, stepSizeB); // smoothen local spectral peak of _�`- shape

      stepSizes[b - 2] = uint32_t ((modifiedB * stepSizes[b - 2] + oneMinusB * __min (maxM3M0, stepSizes[b - 2]) + 64) >> 7); // _-`-
      stepSizes[b - 1] = uint32_t ((modifiedB * stepSizes[b - 1] + oneMinusB * __min (maxM3M0, stepSizes[b - 1]) + 64) >> 7); // _---
    }
    stepSizeM3 = stepSizeM2;
    stepSizeM2 = stepSizeM1;
    stepSizeM1 = (stepSizes[b] = stepSizeB); // modified step-size may be smoothened in next loop iteration
  }
}
// ADD: Helper function inspired by libaac for psychoacoustic model
// Function to convert a spectral line index to the Bark scale
static double barcLineValue(int n, int samplerate) {
    // Assuming 1024-line transform for frequency calculation
    double center_freq = (double)n * samplerate / (2.0 * 1024.0);
    return 13.0 * atan(0.00076 * center_freq) + 3.5 * atan(pow(center_freq / 7500.0, 2.0));
}


// constructor
BitAllocator::BitAllocator ()
{
  for (unsigned ch = 0; ch < USAC_MAX_NUM_CHANNELS; ch++)
  {
    m_avgStepSize[ch] = 0;
    m_avgSpecFlat[ch] = 0;
    m_avgTempFlat[ch] = 0;
  }
  m_rateIndex    = 0;
  m_tempSfbValue = nullptr;
  m_tnsPredictor = nullptr;
}

// ADD: New private method implementation for Minimum SNR
/*
 * initMinSnr: Calculates the minimum signal-to-noise ratio for each SFB.
 * This is a floating-point adaptation of FDKaacEnc_initMinSnr.
 */
void BitAllocator::initMinSnr(const unsigned nChannels, const unsigned samplingRate, const long bitrate, const SfbGroupData* const groupData[USAC_MAX_NUM_CHANNELS])
{
    for (unsigned ch = 0; ch < nChannels; ch++)
    {
        const SfbGroupData& grpData = *groupData[ch];
        // Determine transform size based on window groups. Default to 1024 for long blocks.
        const int numLines = (grpData.numWindowGroups == 1) ? 1024 : 128;
        const int sfbActive = grpData.sfbsPerGroup;
        const uint16_t* sfbOffset = grpData.sfbOffsets;
        const bool isShortBlock = (grpData.numWindowGroups != 1);

        if (sfbActive == 0) continue;

        // Constants adapted from libaaccrop
        const double MAX_BARC = 24.0;
        const double MAX_BARCP1 = 25.0;
        const double BITS2PEFAC = 1.44269504089; // log2(e)
        const double PERS2P4 = 0.024;
        const double MAX_SNR = 0.8;
        const double MIN_SNR = 0.003;

        double barcFactor = std::min(barcLineValue(sfbOffset[sfbActive], samplingRate), MAX_BARC) / MAX_BARCP1;
        
        double pePerWindow = (double)bitrate / samplingRate;
        pePerWindow *= BITS2PEFAC * PERS2P4;

        // Adjust pePerWindow based on transform size
        switch (numLines) {
            case 1024: pePerWindow /= 1024.0; break;
            case 128:  pePerWindow /= 128.0;  break;
            case 512:  pePerWindow /= 512.0;  break;
            case 480:  pePerWindow /= 480.0;  break;
        }

        if (isShortBlock) {
            pePerWindow *= 1.5;
        }

        double pePart_const = pePerWindow / barcFactor;

        for (int sfb = 0; sfb < sfbActive; sfb++) {
            double barcWidth = barcLineValue(sfbOffset[sfb + 1], samplingRate) - barcLineValue(sfbOffset[sfb], samplingRate);
            double pePart = pePart_const * barcWidth;
            int sfbWidth = sfbOffset[sfb + 1] - sfbOffset[sfb];

            if (sfbWidth > 0) {
                pePart /= sfbWidth;
            }

            double tmp = pow(2.0, pePart);
            double snr = (tmp / 2.0) - 1.5;
            snr = std::max(1.0, snr);
            snr = 1.0 / snr;

            snr = std::max(MIN_SNR, std::min(MAX_SNR, snr));
            
            // Store the final SNR value for this band
            m_sfbMinSnr[ch][sfb] = snr;
        }
    }
}

// public functions
void BitAllocator::getChAverageSpecFlat (uint8_t meanSpecFlatInCh[USAC_MAX_NUM_CHANNELS], const unsigned nChannels)
{
  if ((meanSpecFlatInCh == nullptr) || (nChannels > USAC_MAX_NUM_CHANNELS))
  {
    return;
  }
  memcpy (meanSpecFlatInCh, m_avgSpecFlat, nChannels * sizeof (uint8_t));
}

void BitAllocator::getChAverageTempFlat (uint8_t meanTempFlatInCh[USAC_MAX_NUM_CHANNELS], const unsigned nChannels)
{
  if ((meanTempFlatInCh == nullptr) || (nChannels > USAC_MAX_NUM_CHANNELS))
  {
    return;
  }
  memcpy (meanTempFlatInCh, m_avgTempFlat, nChannels * sizeof (uint8_t));
}

uint16_t BitAllocator::getRateCtrlFac (const int32_t rateRatio, const unsigned samplingRate, const uint32_t specFlatness,
                                       const bool prevEightShorts /*= false*/)
{
#if BA_MORE_CBR
  const int32_t ratioFac = rateRatio * (40 - 5 * m_rateIndex);
  const uint32_t brRatio = __max ((prevEightShorts ? (ratioFac * ratioFac + (1 << 16)) >> 17 : 0) - SHRT_MIN, __min (USHRT_MAX, ratioFac)) -
                           (m_rateIndex == 2 ? 1 << 12 : 0);  // rate tuning
  const uint16_t mSfmSqr = (m_rateIndex <= 2 && samplingRate >= 27713 ? (specFlatness * specFlatness) >> m_rateIndex : 0);
#else
  const uint32_t brRatio = __max (1 << 15, __min ((1 << 16) - 1, rateRatio * (36 - 9 * m_rateIndex)));
  const uint16_t mSfmSqr = (m_rateIndex < 2 && samplingRate >= 27713 ? (specFlatness * specFlatness) >> m_rateIndex : 0);
#endif
  const uint16_t mSfmFac = 256 - (((32 + m_rateIndex) * (specFlatness << 4) - mSfmSqr + (1 << 9)) >> 10);

  return uint16_t ((brRatio * mSfmFac + (1 << 7)) >> 8);
}

uint8_t BitAllocator::getScaleFac (const uint32_t sfbStepSize, const int32_t* const sfbSignal, const uint8_t sfbWidth,
                                   const uint32_t sfbRmsValue)
{
  uint8_t sf;
  uint32_t u;
#if !RESTRICT_TO_AAC
  uint64_t meanSpecLoudness = 0;
  double d;
#endif

  if ((sfbSignal == nullptr) || (sfbWidth == 0) || (sfbRmsValue < 46))
  {
    return 0; // use lowest scale factor
  }
#if RESTRICT_TO_AAC
  u = 0;
  for (sf = 0; sf < sfbWidth; sf++)
  {
    u += intSqrt (abs (sfbSignal[sf]));
  }
  u = uint32_t ((u * 16384ui64 + (sfbWidth >> 1)) / sfbWidth);
  u = uint32_t (0.5 + sqrt ((double) u) * 128.0);

  if (u < 42567) return 0;

  u = uint32_t ((sfbStepSize * 42567ui64 + (u >> 1)) / u);
  sf = (u > 1 ? uint8_t (0.5 + 17.7169498394 * log10 ((double) u)) : 4);
#else
  for (sf = 0; sf < sfbWidth; sf++) // simple, low-complexity derivation method for USAC's arithmetic coder
  {
    const int64_t temp = ((int64_t) sfbSignal[sf] + 8) >> 4; // avoid overflow

    meanSpecLoudness += temp * temp;
  }
  meanSpecLoudness = uint64_t (0.5 + pow (256.0 * (double) meanSpecLoudness / sfbWidth, 0.25));

  u = uint32_t (0.5 + pow ((double) sfbRmsValue, 0.75) * 256.0);    // u = 2^8 * (sfbRmsValue^0.75)
  u = uint32_t ((meanSpecLoudness * sfbStepSize * 665 + (u >> 1)) / u); // u = sqrt(6.75) * m*thr/u
  d =  (u > 1 ? log10 ((double) u) : 0.25);

  u = uint32_t (0.5 + pow ((double) sfbRmsValue, 0.25) * 16384.0); // u = 2^14 * (sfbRmsValue^0.25)
  u = uint32_t (((uint64_t) sfbStepSize * 42567 + (u >> 1)) / u);         // u = sqrt(6.75) * thr/u
  d += (u > 1 ? log10 ((double) u) : 0.25);

  sf = uint8_t (0.5 + 8.8584749197 * d);  // sf = (8/3) * log2(u1*u2) = (8/3) * (log2(u1)+log2(u2))
#endif

  return __min (SCHAR_MAX, sf);
}

unsigned BitAllocator::initAllocMemory (LinearPredictor* const linPredictor, const uint8_t numSwb, const uint8_t bitRateMode)
{
  if (linPredictor == nullptr)
  {
    return 1; // invalid arguments error
  }
  m_rateIndex    = bitRateMode;
  m_tnsPredictor = linPredictor;

  if ((m_tempSfbValue = (uint8_t*) malloc (__max (MAX_PREDICTION_ORDER * sizeof (short), numSwb) * sizeof (uint8_t))) == nullptr)
  {
    return 2; // memory allocation error
  }

  return 0; // no error
}

// MODIFY: Add bitrate parameter to signature
unsigned BitAllocator::initSfbStepSizes (const SfbGroupData* const groupData[USAC_MAX_NUM_CHANNELS], const uint8_t numSwbShort,
                                         const uint32_t specAnaStats[USAC_MAX_NUM_CHANNELS],
                                         const uint32_t tempAnaStats[USAC_MAX_NUM_CHANNELS],
                                         const unsigned nChannels, const unsigned samplingRate, uint32_t* const sfbStepSizes,
                                         const unsigned lfeChannelIndex, const unsigned ad /*= 0u*/, const bool tnsDisabled /*= false*/, const long bitrate /*= 64000*/)
{
      // ADD: Initialize the Minimum SNR model
      initMinSnr(nChannels, samplingRate, bitrate, groupData);

      // equal-loudness weighting based on data from: K. Kurakata, T. Mizunami, and K. Matsushita, "Percentiles
      // of Normal Hearing-Threshold Distribution Under Free-Field Listening Conditions in Numerical Form," Ac.
      // Sci. Tech, vol. 26, no. 5, pp. 447-449, Jan. 2005, https://www.researchgate.net/publication/239433096.
        
        // Clear tonality indices at start
        memset(m_tonalityIndex, 0, sizeof(m_tonalityIndex));

        
      const unsigned HF/*idx*/= ((123456 - samplingRate) >> 11) + (samplingRate < 37566 ? 2 : ad); // start SFB
      const unsigned LF/*idx*/= 9;
      const unsigned MF/*idx*/= (samplingRate < 27713 ? HF : __min (HF, 30u));
      const unsigned msShift  = (samplingRate + 36736) >> 15; // TODO: 768 smp
      const unsigned msOffset = 1 << (msShift - 1);
      uint32_t nMeans = 0, sumMeans = 0;

      if ((groupData == nullptr) || (specAnaStats == nullptr) || (tempAnaStats == nullptr) || (sfbStepSizes == nullptr) ||
          (numSwbShort < MIN_NUM_SWB_SHORT) || (numSwbShort > MAX_NUM_SWB_SHORT) || (nChannels > USAC_MAX_NUM_CHANNELS) ||
          (samplingRate < 7350) || (samplingRate > 96000) || (lfeChannelIndex > USAC_MAX_NUM_CHANNELS))
      {
        return 1; // invalid arguments error
      }

      for (unsigned ch = 0; ch < nChannels; ch++)
      {
        const SfbGroupData& grpData = *groupData[ch];
        const uint32_t maxSfbInCh = __min (MAX_NUM_SWB_LONG, grpData.sfbsPerGroup);
        const uint32_t nBandsInCh = grpData.numWindowGroups * maxSfbInCh;
        const uint32_t* rms = grpData.sfbRmsValues;
        uint32_t* stepSizes = &sfbStepSizes[ch * numSwbShort * NUM_WINDOW_GROUPS];
    // --- apply INTRA-channel simultaneous masking, equal-loudness weighting, and thresholding to SFB RMS data
        uint32_t maskingSlope = 0, gr, b, elw = 58254; // = 64k*8/9
        uint32_t rmsEqualLoud = 0;
        uint32_t sumStepSizes = 0;

        m_avgStepSize[ch] = 0;

        b = ((specAnaStats[ch] >> 16) & UCHAR_MAX);
        b = __max (b * b, (tempAnaStats[ch] >> 24) * (tempAnaStats[ch] >> 24));
        m_avgSpecFlat[ch] = uint8_t ((b + (1 << 7)) >> 8); // max. of squared SFM from spec. and temp. analysis

        b = ((tempAnaStats[ch] >> 16) & UCHAR_MAX);
        b = __max (b * b, (specAnaStats[ch] >> 24) * (specAnaStats[ch] >> 24));
        m_avgTempFlat[ch] = uint8_t ((b + (1 << 7)) >> 8); // max. of squared TFM from spec. and temp. analysis

        if ((nBandsInCh == 0) || (grpData.numWindowGroups > NUM_WINDOW_GROUPS))
        {
          continue;
        }
        if ((ch == lfeChannelIndex) || (grpData.numWindowGroups != 1)) // LFE, SHORT windows: no masking or ELW
        {
          for (gr = 0; gr < grpData.numWindowGroups; gr++)
          {
            const uint32_t* gRms = &rms[numSwbShort * gr];
            uint32_t* gStepSizes = &stepSizes[numSwbShort * gr];

            for (b = numSwbShort - 1; b >= maxSfbInCh; b--)
            {
              gStepSizes[b] = 0;
            }
            for (/*b*/; b > 0; b--)
            {
              gStepSizes[b] = __max (gRms[b], BA_EPS);
              sumStepSizes += intSqrt (gStepSizes[b]);
            }
            gStepSizes[0]   = __max (gRms[0], BA_EPS);
            sumStepSizes   += intSqrt (gStepSizes[0]);
          } // for gr

          if (ch != lfeChannelIndex)
          {
    // --- SHORT windows: apply perceptual just noticeable difference (JND) model and local band-peak smoothing
            nMeans++;

            for (b = maxSfbInCh - 1; b > 0; b--) // gentle temporal band-peak smoothing; a spectral one follows
            {
              uint32_t maxGrpStep = stepSizes[b], stepSizeM1 = BA_EPS;

              for (gr = 1; gr < grpData.numWindowGroups; gr++)
              {
                const uint32_t curGrpStep = stepSizes[b + numSwbShort * gr];

                if (curGrpStep > maxGrpStep) maxGrpStep = curGrpStep;
              }
              for (gr = 0; gr < grpData.numWindowGroups; gr++)
              {
                const uint32_t newGrpStep = __max (stepSizeM1, (gr + 1 == grpData.numWindowGroups ? BA_EPS : stepSizes[b + numSwbShort * (gr + 1)]));

                stepSizeM1 = stepSizes[b + numSwbShort * gr];

                if ((stepSizeM1 == maxGrpStep) && (maxGrpStep > newGrpStep))
                {
                  const uint32_t sqrtOldStep = intSqrt (maxGrpStep);
                  const uint32_t sqrtNewStep = intSqrt (newGrpStep);
                  uint32_t& gStepSize = stepSizes[b + numSwbShort * gr];

                  sumStepSizes += (gStepSize = (sqrtOldStep + sqrtNewStep) >> 1) - sqrtOldStep;
                  gStepSize *= gStepSize; // for square-mean-root
                }
              }
            } // for b

            m_avgStepSize[ch] = __min (USHRT_MAX, (sumStepSizes + (nBandsInCh >> 1)) / nBandsInCh);
            sumMeans += m_avgStepSize[ch];
            m_avgStepSize[ch] *= m_avgStepSize[ch];

            for (gr = 0; gr < grpData.numWindowGroups; gr++) // separate spectral peak smoothing for each group
            {
              jndPowerLawAndPeakSmoothing (&stepSizes[numSwbShort * gr], maxSfbInCh, false,
                                           m_avgStepSize[ch], m_avgSpecFlat[ch], 0);
            }
          }
          continue;
        }

        stepSizes[0]   = __max (rms[0], BA_EPS);
        for (b = 1; b < __min (LF, maxSfbInCh); b++) // apply steeper low-frequency simultaneous masking slopes
        {
          maskingSlope = (stepSizes[b - 1] + (msOffset << (9u - b))) >> (msShift + 9u - b);
          stepSizes[b] = __max (rms[b], maskingSlope + BA_EPS);
        }
        for (/*b*/; b < __min (MF, maxSfbInCh); b++) // apply typical mid-frequency simultaneous masking slopes
        {
          maskingSlope = (stepSizes[b - 1] + msOffset) >> msShift;
          stepSizes[b] = __max (rms[b], maskingSlope + BA_EPS);
        }
        if ((samplingRate >= 27713) && (samplingRate < 75132))
        {
          for (/*b*/; b < __min (HF, maxSfbInCh); b++) // compensate high-frequency slopes for linear SFB width
          {
            maskingSlope = ((uint64_t) stepSizes[b - 1] * (9u + b - MF) + (msOffset << 3u)) >> (msShift + 3u);
            stepSizes[b] = __max (rms[b], maskingSlope + BA_EPS);
          }
          for (/*b = HF region*/; b < maxSfbInCh; b++) // apply extra high-frequency equal-loudness attenuation
          {
            for (unsigned d = b - HF; d > 0; d--)
            {
              elw = (elw * 52430 - SHRT_MIN) >> 16; // elw *= 4/5
            }
            rmsEqualLoud = uint32_t (((uint64_t) rms[b] * elw - SHRT_MIN) >> 16);   // equal loudness weighting
            maskingSlope = ((uint64_t) stepSizes[b - 1] * (9u + b - MF) + (msOffset << 3u)) >> (msShift + 3u);
            stepSizes[b] = __max (rmsEqualLoud, maskingSlope + BA_EPS);
          }
        }
        else // no equal-loudness weighting for low or high rates
        {
          for (/*b = MF region*/; b < maxSfbInCh; b++) // compensate high-frequency slopes for linear SFB width
          {
            maskingSlope = ((uint64_t) stepSizes[b - 1] * (9u + b - MF) + (msOffset << 3u)) >> (msShift + 3u);
            stepSizes[b] = __max (rms[b], maskingSlope + BA_EPS);
          }
        }
        stepSizes[b] = 0;
        for (b -= 1; b > __min (MF, maxSfbInCh); b--) // complete simultaneous masking by reversing the pattern
        {
          sumStepSizes += intSqrt (stepSizes[b]);
          maskingSlope     = ((uint64_t) stepSizes[b] * (8u + b - MF) + (msOffset << 3u)) >> (msShift + 3u);
          stepSizes[b - 1] = __max (stepSizes[b - 1], maskingSlope);
        }
        for (/*b*/; b > __min (LF, maxSfbInCh); b--)  // typical reversed mid-freq. simultaneous masking slopes
        {
          sumStepSizes += intSqrt (stepSizes[b]);
          maskingSlope     = (stepSizes[b] + msOffset) >> msShift;
          stepSizes[b - 1] = __max (stepSizes[b - 1], maskingSlope);
        }
        for (/*b = min (9, maxSfbInCh)*/; b > 0; b--) // steeper reversed low-freq. simultaneous masking slopes
        {
          sumStepSizes += intSqrt (stepSizes[b]);
          maskingSlope     = (stepSizes[b] + (msOffset << (10u - b))) >> (msShift + 10u - b);
          stepSizes[b - 1] = __max (stepSizes[b - 1], maskingSlope);
        }
        sumStepSizes   += intSqrt (stepSizes[0]);

    // --- LONG window: apply perceptual JND model and local band-peak smoothing, undo equal-loudness weighting
        nMeans++;
        m_avgStepSize[ch] = __min (USHRT_MAX, (sumStepSizes + (nBandsInCh >> 1)) / nBandsInCh);
        sumMeans += m_avgStepSize[ch];
        m_avgStepSize[ch] *= m_avgStepSize[ch];

        jndPowerLawAndPeakSmoothing (stepSizes, maxSfbInCh, (m_rateIndex == 0) && (samplingRate >= 27713),
                                     m_avgStepSize[ch], m_avgSpecFlat[ch], tnsDisabled ? m_avgTempFlat[ch] : 0);

        if ((samplingRate >= 27713) && (samplingRate < 75132))
        {
          elw = 36; // 36/32 = 9/8
          for (b = HF; b < maxSfbInCh; b++)  // undo above additional high-frequency equal-loudness attenuation
          {
            for (unsigned d = b - HF; d > 0; d--)
            {
              elw = (16u + elw * 40) >> 5; // inverse elw *= 5/4. NOTE: this may overflow for 64 kHz, that's OK
            }
            if (elw == 138 || elw >= 1024) elw--;
            rmsEqualLoud = uint32_t (__min (UINT_MAX, (16u + (uint64_t) stepSizes[b] * elw) >> 5));
            stepSizes[b] = rmsEqualLoud;
          }
        }
      } // for ch

      if ((nMeans < 2) || (sumMeans <= nMeans * BA_EPS)) // in case of one channel or low-RMS input, we're done
      {
        return 0; // no error
      }

      sumMeans = (sumMeans + (nMeans >> 1)) / nMeans;
      sumMeans *= sumMeans;  // since we've averaged square-roots

      for (unsigned ch = 0; ch < nChannels; ch++)
      {
        const SfbGroupData& grpData = *groupData[ch];
        const uint32_t maxSfbInCh = __min (MAX_NUM_SWB_LONG, grpData.sfbsPerGroup);
        const uint32_t nBandsInCh = grpData.numWindowGroups * maxSfbInCh;
        const uint32_t chStepSize = m_avgStepSize[ch];
        uint32_t* stepSizes = &sfbStepSizes[ch * numSwbShort * NUM_WINDOW_GROUPS];
    // --- apply INTER-channel simultaneous masking and JND modeling to calculated INTRA-channel step-size data
        uint64_t mAvgStepSize; // modified and averaged step-size

        if ((nBandsInCh == 0) || (grpData.numWindowGroups > NUM_WINDOW_GROUPS) || (ch == lfeChannelIndex))
        {
          continue;
        }
        mAvgStepSize = jndModel (chStepSize, sumMeans, 7 << 6 /*7/8*/, 512);

        for (unsigned gr = 0; gr < grpData.numWindowGroups; gr++)
        {
          uint32_t* gStepSizes = &stepSizes[numSwbShort * gr];

          for (unsigned b = 0; b < maxSfbInCh; b++)
          {
            gStepSizes[b] = uint32_t (__min (UINT_MAX, (mAvgStepSize * gStepSizes[b] + (chStepSize >> 1)) / chStepSize));
          }
        }

        m_avgStepSize[ch] = (uint32_t) mAvgStepSize;
      } // for ch

      return 0; // no error
    }

    unsigned BitAllocator::imprSfbStepSizes (const SfbGroupData* const groupData[USAC_MAX_NUM_CHANNELS], const uint8_t numSwbShort,
                                             const int32_t* const mdctSpec[USAC_MAX_NUM_CHANNELS], const unsigned nSamplesInFrame,
                                             const unsigned nChannels, const unsigned samplingRate, uint32_t* const sfbStepSizes,
                                             const unsigned firstChannelIndex, const uint8_t* const sfm, const bool commonWindow,
                                             const uint8_t* const sfbStereoData /*= nullptr*/, const uint8_t stereoConfig /*= 0*/)
    {
      const uint8_t maxSfbL16k = 16 + __min (4 + (samplingRate >> 10), (9 << 17) / __max (1, samplingRate)); // SFB index at 15.8 kHz
      const uint32_t redFactor = __max ((samplingRate < 25495 ? 2 : 1), __min (3, m_rateIndex)) - (stereoConfig >> 3);
      const uint32_t redWeight = __min (4, 9 - __min (9, m_rateIndex));
      short* const  tempCoeffs = (short* const) m_tempSfbValue;

      if ((groupData == nullptr) || (mdctSpec == nullptr) || (sfbStepSizes == nullptr) || (sfm == nullptr) || (nSamplesInFrame > 2048) ||
          (numSwbShort < MIN_NUM_SWB_SHORT) || (numSwbShort > MAX_NUM_SWB_SHORT) || (nChannels > USAC_MAX_NUM_CHANNELS) ||
          (samplingRate < 7350) || (samplingRate > 96000) || (firstChannelIndex > USAC_MAX_NUM_CHANNELS))
      {
        return 1; // invalid arguments error
      }

      for (unsigned ch = firstChannelIndex; ch < firstChannelIndex + nChannels; ch++)
      {
        const SfbGroupData& grpData = *groupData[ch];
        const uint32_t maxSfbInCh = __min (MAX_NUM_SWB_LONG, grpData.sfbsPerGroup);
        const bool    eightShorts = (grpData.numWindowGroups != 1);
        const bool  lowRateTuning = (m_rateIndex == 0) && (samplingRate >= 25495 && sfm[ch] <= (SCHAR_MAX >> 1));
        const bool undercodingRed = (m_rateIndex >  0) || (samplingRate >= 25495 && sfm[ch] * 8 > UCHAR_MAX * 7) || lowRateTuning;
        const uint32_t* rms = grpData.sfbRmsValues;
        uint32_t* stepSizes = &sfbStepSizes[ch * numSwbShort * NUM_WINDOW_GROUPS];

        if ((grpData.numWindowGroups * maxSfbInCh == 0) || (grpData.numWindowGroups > NUM_WINDOW_GROUPS))
        {
          continue;
        }
        for (unsigned gr = 0; gr < grpData.numWindowGroups; gr++)
        {
          const uint16_t* grpOff = &grpData.sfbOffsets[numSwbShort * gr];
          const uint8_t* grpSte = (sfbStereoData == nullptr ? nullptr : &sfbStereoData[numSwbShort * gr]);
          const uint32_t* grpRms = &rms[numSwbShort * gr];
          const uint32_t* refRms = &groupData[firstChannelIndex + nChannels - 1 - ch]->sfbRmsValues[numSwbShort * gr];
          uint32_t* grpStepSizes = &stepSizes[numSwbShort * gr];
          uint32_t  b, grpRmsMin = INT_MAX; // min. RMS value, used for overcoding reduction
          uint64_t  s = (eightShorts ? (nSamplesInFrame * grpData.windowGroupLength[gr]) >> 1 : nSamplesInFrame << 2);

          memset (m_tempSfbValue, UCHAR_MAX, maxSfbInCh * sizeof (uint8_t));
          if (lowRateTuning && (maxSfbInCh > 0) && !eightShorts)
          {
            uint32_t numRedBands = nSamplesInFrame; // final result lies between 1/4 and 1/2

            if ((nChannels == 2) && commonWindow && (grpSte != nullptr))
            {
              for (b = 0; b < maxSfbInCh; b++) if (grpSte[b] == 0) numRedBands += grpOff[b + 1] - grpOff[b];
            }
            b = MAX_NUM_SWB_LONG - ((numRedBands * ((SCHAR_MAX >> 1) + 1 - sfm[ch]) + (1 << 11)) >> 12);

            while ((b < maxSfbInCh) && (grpRms[b] > grpRms[b - 1])) b++; // start after peak

            for (b += ((nChannels == 2) && commonWindow ? b & 1 : 0); b < maxSfbInCh; b++)
            {
              grpStepSizes[b] = __max (grpStepSizes[b], grpRms[b] >= (UINT_MAX >> 1) ? UINT_MAX : (grpRms[b] + 1) << 1);
            }
          }

          // undercoding reduction for case where large number of coefs is quantized to zero
          for (b = 0; b < maxSfbInCh; b++)
          {
            const uint32_t rmsComp = (grpSte != nullptr && grpSte[b] > 0 ? squareMeanRoot (refRms[b], grpRms[b]) : grpRms[b]);
            const uint32_t rmsRef9 = (commonWindow ? refRms[b] >> 9 : rmsComp);
            const uint8_t sfbWidth = grpOff[b + 1] - grpOff[b];

            if (redWeight > 0 && !eightShorts && sfbWidth > (samplingRate >= 18783 ? 8 : 12)) // transient SFBs
            {
              const uint32_t gains = m_tnsPredictor->calcParCorCoeffs (&mdctSpec[ch][grpOff[b]], sfbWidth, MAX_PREDICTION_ORDER, tempCoeffs) >> 24;

              m_tempSfbValue[b] = UCHAR_MAX - uint8_t ((512u + gains * gains * redWeight) >> (10 + (sfbWidth > 16 ? 0 : (20 - sfbWidth) >> 2)));
              if ((b >= 2) && (m_tempSfbValue[b - 1] < m_tempSfbValue[b]) && (m_tempSfbValue[b - 1] < m_tempSfbValue[b - 2]))
              {
                m_tempSfbValue[b - 1] = __min (m_tempSfbValue[b], m_tempSfbValue[b - 2]); // remove local peaks
              }
            }
            if (grpRms[b] < grpRmsMin) grpRmsMin = grpRms[b];

            if (undercodingRed && (rmsComp >= rmsRef9) && (rmsComp < (grpStepSizes[b] >> 1))) // zero-quantized
            {
              s -= (sfbWidth * redFactor * __min (1u << 11, rmsComp) + (1u << 10)) >> 11;
            }
          }

          if ((samplingRate >= 27713) && (b < maxSfbL16k) && !eightShorts) // zeroed HF data
          {
            const uint32_t rmsComp = (grpSte != nullptr && grpSte[b] > 0 ? squareMeanRoot (refRms[b], grpRms[b]) : grpRms[b]);
            const uint32_t rmsRef9 = (commonWindow ? refRms[b] >> 9 : rmsComp);
            const uint8_t sfbWidth = grpOff[maxSfbL16k] - grpOff[b];

            if (undercodingRed && (rmsComp >= rmsRef9)) // check only first SFB above max_sfb as simplification
            {
              s -= (sfbWidth * redFactor * __min (1u << 11, rmsComp) + (1u << 10)) >> 11;
            }
          }
          s = (eightShorts ? s / ((nSamplesInFrame * grpData.windowGroupLength[gr]) >> 8) : s / (nSamplesInFrame >> 5));

          if (redWeight > 0 && !eightShorts) memset (tempCoeffs /*= m_tempSfbValue*/, UCHAR_MAX, MAX_PREDICTION_ORDER * sizeof (short));

          for (b = 0; b < maxSfbInCh; b++) // improve step-sizes by limiting and attenuation
          {
            grpStepSizes[b] = uint32_t ((__max (grpRmsMin, grpStepSizes[b]) * s * (m_tempSfbValue[b] + 1u) + (1u << 14)) >> 15);
            if (grpStepSizes[b] <= (grpRms[b] >> 11)) grpStepSizes[b] = __max (BA_EPS, grpRms[b] >> 11);

            if (lowRateTuning) // clip near-0 SNRs to minimum SNR
            {
              uint32_t lim = uint32_t ((grpRms[b] * (8192u - (uint64_t) sfm[ch] * sfm[ch]) + (1u << 12)) >> 13);

              // ADD: Enforce minimum SNR floor from libaaccrop model
              if (m_sfbMinSnr[ch][b] > 0) {
                  // Calculate the maximum allowed step size to maintain the minimum SNR
                  uint32_t minSnrLim = (uint32_t)((double)grpRms[b] / m_sfbMinSnr[ch][b]);
                  // The step size should not be larger than this limit
                  grpStepSizes[b] = __min(grpStepSizes[b], minSnrLim);
              }

              if ((grpStepSizes[b] > grpRms[b]) && ((grpStepSizes[b] >> 1) <= lim)) grpStepSizes[b] = grpRms[b];
            }
          }
        }
      } // for ch

      return 0; // no error
        
    }

    // Estimate tonality of a scale factor band (0 = noise-like, 1 = tonal)
    float BitAllocator::estimateTonality(const int32_t* mdctSpec,
                                       const uint16_t* sfbOffsets,
                                       const unsigned sfb,
                                       const unsigned sfbWidth)
    {
      if (sfbWidth == 0 || mdctSpec == nullptr || sfbOffsets == nullptr) {
        return 0.5f; // Default neutral tonality
      }
      
      const unsigned startLine = sfbOffsets[sfb];
      double geometricMean = 0.0;
      double arithmeticMean = 0.0;
      unsigned nonZeroCount = 0;
      
      // Calculate geometric and arithmetic means of magnitudes
      for (unsigned i = 0; i < sfbWidth; i++) {
        // Use absolute value of MDCT coefficient
        double magnitude = fabs((double)mdctSpec[startLine + i]);
        
        if (magnitude > BA_EPS) { // BA_EPS is already defined as 1
          // Add small epsilon to avoid log(0)
          geometricMean += log(magnitude + 1e-10);
          arithmeticMean += magnitude;
          nonZeroCount++;
        }
      }
      
      // If all coefficients are zero, it's noise-like
      if (nonZeroCount == 0) {
        return 0.0f;
      }
      
      // Complete the mean calculations
      geometricMean = exp(geometricMean / nonZeroCount);
      arithmeticMean = arithmeticMean / nonZeroCount;
      
      // Spectral Flatness Measure (SFM)
      // SFM = geometric_mean / arithmetic_mean
      // SFM close to 1 = noise-like, SFM close to 0 = tonal
      float sfm = (float)(geometricMean / (arithmeticMean + 1e-10));
      
      // Convert SFM to tonality index (0 = noise, 1 = tonal)
      // Clamp between 0 and 1
      float tonality = 1.0f - sfm;
      tonality = __max(0.0f, __min(1.0f, tonality));
      
      // Apply slight smoothing to avoid extreme values
      // This makes the tonality estimation more robust
      if (tonality > 0.9f) {
        tonality = 0.9f + 0.1f * (tonality - 0.9f);
      } else if (tonality < 0.1f) {
        tonality = 0.1f * tonality;
      }
      
      return tonality;
    }
