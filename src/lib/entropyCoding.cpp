/* entropyCoding.cpp - source file for class with lossless entropy coding capability
 * written by C. R. Helmrich, last modified in 2025 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2025 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "entropyCoding.h"

// ISO/IEC 14496-3, Annex 4.A.1
static const uint32_t huffScf[INDEX_SIZE] = { // upper 3 bytes: value, lower byte: length
  0x03FFE812, 0x03FFE612, 0x03FFE712, 0x03FFE512, 0x07FFF513, 0x07FFF113, 0x07FFED13, 0x07FFF613,
  0x07FFEE13, 0x07FFEF13, 0x07FFF013, 0x07FFFC13, 0x07FFFD13, 0x07FFFF13, 0x07FFFE13, 0x07FFF713,
  0x07FFF813, 0x07FFFB13, 0x07FFF913, 0x03FFE412, 0x07FFFA13, 0x03FFE312, 0x01FFEF11, 0x01FFF011,
  0x00FFF510, 0x01FFEE11, 0x00FFF210, 0x00FFF310, 0x00FFF410, 0x00FFF110, 0x007FF60F, 0x007FF70F,
  0x003FF90E, 0x003FF50E, 0x003FF70E, 0x003FF30E, 0x003FF60E, 0x003FF20E, 0x001FF70D, 0x001FF50D,
  0x000FF90C, 0x000FF70C, 0x000FF60C, 0x0007F90B, 0x000FF40C, 0x0007F80B, 0x0003F90A, 0x0003F70A,
  0x0003F50A, 0x0001F809, 0x0001F709, 0x0000FA08, 0x0000F808, 0x0000F608, 0x00007907, 0x00003A06,
  0x00003806, 0x00001A05, 0x00000B04, 0x00000403, 0x00000001, 0x00000A04, 0x00000C04, 0x00001B05, 0x00003906,
  0x00003B06, 0x00007807, 0x00007A07, 0x0000F708, 0x0000F908, 0x0001F609, 0x0001F909, 0x0003F40A,
  0x0003F60A, 0x0003F80A, 0x0007F50B, 0x0007F40B, 0x0007F60B, 0x0007F70B, 0x000FF50C, 0x000FF80C,
  0x001FF40D, 0x001FF60D, 0x001FF80D, 0x003FF80E, 0x003FF40E, 0x00FFF010, 0x007FF40F, 0x00FFF610,
  0x007FF50F, 0x03FFE212, 0x07FFD913, 0x07FFDA13, 0x07FFDB13, 0x07FFDC13, 0x07FFDD13, 0x07FFDE13,
  0x07FFD813, 0x07FFD213, 0x07FFD313, 0x07FFD413, 0x07FFD513, 0x07FFD613, 0x07FFF213, 0x07FFDF13,
  0x07FFE713, 0x07FFE813, 0x07FFE913, 0x07FFEA13, 0x07FFEB13, 0x07FFE613, 0x07FFE013, 0x07FFE113,
  0x07FFE213, 0x07FFE313, 0x07FFE413, 0x07FFE513, 0x07FFD713, 0x07FFEC13, 0x07FFF413, 0x07FFF313
};

// ISO/IEC 23003-3, Annex C.1
static const uint8_t arithLookupM[ARITH_SIZE] = { // arith_lookup_m
  0x01, 0x34, 0x0D, 0x13, 0x12, 0x25, 0x00, 0x3A, 0x05, 0x00, 0x21, 0x13, 0x1F, 0x1A, 0x1D, 0x36,
  0x24, 0x2B, 0x1B, 0x33, 0x37, 0x29, 0x1D, 0x33, 0x37, 0x33, 0x37, 0x33, 0x37, 0x33, 0x2C, 0x00,
  0x21, 0x13, 0x25, 0x2A, 0x00, 0x21, 0x24, 0x12, 0x2C, 0x1E, 0x37, 0x24, 0x1F, 0x35, 0x37, 0x24,
  0x35, 0x37, 0x35, 0x37, 0x38, 0x2D, 0x21, 0x29, 0x1E, 0x21, 0x13, 0x2D, 0x36, 0x38, 0x29, 0x36,
  0x37, 0x24, 0x36, 0x38, 0x37, 0x38, 0x00, 0x20, 0x23, 0x20, 0x23, 0x36, 0x38, 0x24, 0x3B, 0x24,
  0x26, 0x29, 0x1F, 0x30, 0x2D, 0x0D, 0x12, 0x3F, 0x2D, 0x21, 0x1C, 0x2A, 0x00, 0x21, 0x12, 0x1E,
  0x36, 0x38, 0x36, 0x37, 0x3F, 0x1E, 0x0D, 0x1F, 0x2A, 0x1E, 0x21, 0x24, 0x12, 0x2A, 0x3C, 0x21,
  0x24, 0x1F, 0x3C, 0x21, 0x29, 0x36, 0x38, 0x36, 0x37, 0x38, 0x21, 0x1E, 0x00, 0x3B, 0x25, 0x1E,
  0x20, 0x10, 0x1F, 0x3C, 0x20, 0x23, 0x29, 0x08, 0x23, 0x12, 0x08, 0x23, 0x21, 0x38, 0x00, 0x20,
  0x13, 0x20, 0x3B, 0x1C, 0x20, 0x3B, 0x29, 0x20, 0x23, 0x24, 0x21, 0x24, 0x21, 0x24, 0x3B, 0x13,
  0x23, 0x26, 0x23, 0x13, 0x21, 0x24, 0x26, 0x29, 0x12, 0x22, 0x2B, 0x02, 0x1E, 0x0D, 0x1F, 0x2D,
  0x00, 0x0D, 0x12, 0x00, 0x3C, 0x21, 0x29, 0x3C, 0x21, 0x2A, 0x3C, 0x3B, 0x22, 0x1E, 0x20, 0x10,
  0x1F, 0x3C, 0x0D, 0x29, 0x3C, 0x21, 0x24, 0x08, 0x23, 0x20, 0x38, 0x39, 0x3C, 0x20, 0x13, 0x3C,
  0x00, 0x0D, 0x13, 0x1F, 0x3C, 0x09, 0x26, 0x1F, 0x08, 0x09, 0x26, 0x12, 0x08, 0x23, 0x29, 0x20,
  0x23, 0x21, 0x24, 0x20, 0x13, 0x20, 0x3B, 0x16, 0x20, 0x3B, 0x29, 0x20, 0x3B, 0x29, 0x20, 0x3B,
  0x13, 0x21, 0x24, 0x29, 0x0B, 0x13, 0x09, 0x3B, 0x13, 0x09, 0x3B, 0x13, 0x21, 0x3B, 0x13, 0x0D,
  0x26, 0x29, 0x26, 0x29, 0x3D, 0x12, 0x22, 0x28, 0x2E, 0x04, 0x08, 0x13, 0x3C, 0x3B, 0x3C, 0x20,
  0x10, 0x3C, 0x21, 0x07, 0x08, 0x10, 0x00, 0x08, 0x0D, 0x29, 0x08, 0x0D, 0x29, 0x08, 0x09, 0x13,
  0x20, 0x23, 0x39, 0x08, 0x09, 0x13, 0x08, 0x09, 0x16, 0x08, 0x09, 0x10, 0x12, 0x20, 0x3B, 0x3D,
  0x09, 0x26, 0x20, 0x3B, 0x24, 0x39, 0x09, 0x26, 0x20, 0x0D, 0x13, 0x00, 0x09, 0x13, 0x20, 0x0D,
  0x26, 0x12, 0x20, 0x3B, 0x13, 0x21, 0x26, 0x0B, 0x12, 0x09, 0x3B, 0x16, 0x09, 0x3B, 0x3D, 0x09,
  0x26, 0x0D, 0x13, 0x26, 0x3D, 0x1C, 0x12, 0x1F, 0x28, 0x2E, 0x07, 0x0B, 0x08, 0x09, 0x00, 0x39,
  0x0B, 0x08, 0x26, 0x08, 0x09, 0x13, 0x20, 0x0B, 0x39, 0x10, 0x39, 0x0D, 0x13, 0x20, 0x10, 0x12,
  0x09, 0x13, 0x20, 0x3B, 0x13, 0x09, 0x26, 0x0B, 0x09, 0x3B, 0x1C, 0x09, 0x3B, 0x13, 0x20, 0x3B,
  0x13, 0x09, 0x26, 0x0B, 0x16, 0x0D, 0x13, 0x09, 0x13, 0x09, 0x13, 0x26, 0x3D, 0x1C, 0x1F, 0x28,
  0x2E, 0x07, 0x10, 0x39, 0x0B, 0x39, 0x39, 0x13, 0x39, 0x0B, 0x39, 0x0B, 0x39, 0x26, 0x39, 0x10,
  0x20, 0x3B, 0x16, 0x20, 0x10, 0x09, 0x26, 0x0B, 0x13, 0x09, 0x13, 0x26, 0x1C, 0x0B, 0x3D, 0x1C,
  0x1F, 0x28, 0x2B, 0x07, 0x0C, 0x39, 0x0B, 0x39, 0x0B, 0x0C, 0x0B, 0x26, 0x0B, 0x26, 0x3D, 0x0D,
  0x1C, 0x14, 0x28, 0x2B, 0x39, 0x0B, 0x0C, 0x0E, 0x3D, 0x1C, 0x0D, 0x12, 0x22, 0x2B, 0x07, 0x0C,
  0x0E, 0x3D, 0x1C, 0x10, 0x1F, 0x2B, 0x0C, 0x0E, 0x19, 0x14, 0x10, 0x1F, 0x28, 0x0C, 0x0E, 0x19,
  0x14, 0x26, 0x22, 0x2B, 0x0C, 0x0E, 0x19, 0x14, 0x26, 0x28, 0x0E, 0x19, 0x14, 0x26, 0x28, 0x0E,
  0x19, 0x14, 0x28, 0x0E, 0x19, 0x14, 0x22, 0x28, 0x2B, 0x0E, 0x14, 0x2B, 0x31, 0x00, 0x3A, 0x3A,
  0x05, 0x05, 0x1B, 0x1D, 0x33, 0x06, 0x35, 0x35, 0x20, 0x21, 0x37, 0x21, 0x24, 0x05, 0x1B, 0x2C,
  0x2C, 0x2C, 0x06, 0x34, 0x1E, 0x34, 0x00, 0x08, 0x36, 0x09, 0x21, 0x26, 0x1C, 0x2C, 0x00, 0x02,
  0x02, 0x02, 0x3F, 0x04, 0x04, 0x04, 0x34, 0x39, 0x20, 0x0A, 0x0C, 0x39, 0x0B, 0x0F, 0x07, 0x07,
  0x07, 0x07, 0x34, 0x39, 0x39, 0x0A, 0x0C, 0x39, 0x0C, 0x0F, 0x07, 0x07, 0x07, 0x00, 0x39, 0x39,
  0x0C, 0x0F, 0x07, 0x07, 0x39, 0x0C, 0x0F, 0x07, 0x39, 0x0C, 0x0F, 0x39, 0x39, 0x0C, 0x0F, 0x39,
  0x0C, 0x39, 0x0C, 0x0F, 0x00, 0x11, 0x27, 0x17, 0x2F, 0x27, 0x00, 0x27, 0x17, 0x00, 0x11, 0x17,
  0x00, 0x11, 0x17, 0x11, 0x00, 0x27, 0x15, 0x11, 0x17, 0x01, 0x15, 0x11, 0x15, 0x11, 0x15, 0x15,
  0x17, 0x00, 0x27, 0x01, 0x27, 0x27, 0x15, 0x00, 0x27, 0x11, 0x27, 0x15, 0x15, 0x15, 0x27, 0x15,
  0x15, 0x15, 0x15, 0x17, 0x2F, 0x11, 0x17, 0x27, 0x27, 0x27, 0x11, 0x27, 0x15, 0x27, 0x27, 0x15,
  0x15, 0x27, 0x17, 0x2F, 0x27, 0x17, 0x2F, 0x27, 0x17, 0x2F, 0x27, 0x17, 0x2F, 0x27, 0x17, 0x2F,
  0x27, 0x17, 0x2F, 0x27, 0x17, 0x2F, 0x27, 0x17, 0x2F, 0x27, 0x17, 0x2F, 0x27, 0x17, 0x2F, 0x27,
  0x17, 0x2F, 0x27, 0x17, 0x2F, 0x27, 0x17, 0x2F, 0x17, 0x2F, 0x2B, 0x00, 0x27, 0x00, 0x00, 0x11,
  0x15, 0x00, 0x11, 0x11, 0x27, 0x27, 0x15, 0x17, 0x15, 0x17, 0x15, 0x17, 0x27, 0x17, 0x27, 0x17,
  0x27, 0x17, 0x27, 0x17, 0x27, 0x17, 0x27, 0x17, 0x27, 0x17, 0x27, 0x17, 0x27, 0x17, 0x27, 0x17,
  0x27, 0x15, 0x27, 0x27, 0x15, 0x27
};

// ISO/IEC 23003-3, Annex C.2
static const uint32_t arithHashM[ARITH_SIZE] = { // arith_hash_m
  0x00000104, 0x0000030A, 0x00000510, 0x00000716, 0x00000A1F, 0x00000F2E, 0x00011100, 0x00111103,
  0x00111306, 0x00111436, 0x00111623, 0x00111929, 0x00111F2E, 0x0011221B, 0x00112435, 0x00112621,
  0x00112D12, 0x00113130, 0x0011331D, 0x00113535, 0x00113938, 0x0011411B, 0x00114433, 0x00114635,
  0x00114F29, 0x00116635, 0x00116F24, 0x00117433, 0x0011FF0F, 0x00121102, 0x0012132D, 0x00121436,
  0x00121623, 0x00121912, 0x0012213F, 0x0012232D, 0x00122436, 0x00122638, 0x00122A29, 0x00122F2B,
  0x0012322D, 0x00123436, 0x00123738, 0x00123B29, 0x0012411D, 0x00124536, 0x00124938, 0x00124F12,
  0x00125535, 0x00125F29, 0x00126535, 0x0012B837, 0x0013112A, 0x0013131E, 0x0013163B, 0x0013212D,
  0x0013233C, 0x00132623, 0x00132F2E, 0x0013321E, 0x00133521, 0x00133824, 0x0013411E, 0x00134336,
  0x00134838, 0x00135135, 0x00135537, 0x00135F12, 0x00137637, 0x0013FF29, 0x00140024, 0x00142321,
  0x00143136, 0x00143321, 0x00143F25, 0x00144321, 0x00148638, 0x0014FF29, 0x00154323, 0x0015FF12,
  0x0016F20C, 0x0018A529, 0x00210031, 0x0021122C, 0x00211408, 0x00211713, 0x00211F2E, 0x0021222A,
  0x00212408, 0x00212710, 0x00212F2E, 0x0021331E, 0x00213436, 0x00213824, 0x0021412D, 0x0021431E,
  0x00214536, 0x00214F1F, 0x00216637, 0x00220004, 0x0022122A, 0x00221420, 0x00221829, 0x00221F2E,
  0x0022222D, 0x00222408, 0x00222623, 0x00222929, 0x00222F2B, 0x0022321E, 0x00223408, 0x00223724,
  0x00223A29, 0x0022411E, 0x00224436, 0x00224823, 0x00225134, 0x00225621, 0x00225F12, 0x00226336,
  0x00227637, 0x0022FF29, 0x0023112D, 0x0023133C, 0x00231420, 0x00231916, 0x0023212D, 0x0023233C,
  0x00232509, 0x00232929, 0x0023312D, 0x00233308, 0x00233509, 0x00233724, 0x0023413C, 0x00234421,
  0x00234A13, 0x0023513C, 0x00235421, 0x00235F1F, 0x00236421, 0x0023FF29, 0x00240024, 0x0024153B,
  0x00242108, 0x00242409, 0x00242726, 0x00243108, 0x00243409, 0x00243610, 0x00244136, 0x00244321,
  0x00244523, 0x00244F1F, 0x00245423, 0x0024610A, 0x00246423, 0x0024FF29, 0x00252510, 0x00253121,
  0x0025343B, 0x00254121, 0x00254510, 0x00254F25, 0x00255221, 0x0025FF12, 0x00266513, 0x0027F529,
  0x0029F101, 0x002CF224, 0x00310030, 0x0031122A, 0x00311420, 0x00311816, 0x0031212C, 0x0031231E,
  0x00312408, 0x00312710, 0x0031312A, 0x0031321E, 0x00313408, 0x00313623, 0x0031411E, 0x0031433C,
  0x00320007, 0x0032122D, 0x00321420, 0x00321816, 0x0032212D, 0x0032233C, 0x00322509, 0x00322916,
  0x0032312D, 0x00323420, 0x00323710, 0x00323F2B, 0x00324308, 0x00324623, 0x00324F25, 0x00325421,
  0x00325F1F, 0x00326421, 0x0032FF29, 0x00331107, 0x00331308, 0x0033150D, 0x0033211E, 0x00332308,
  0x00332420, 0x00332610, 0x00332929, 0x0033311E, 0x00333308, 0x0033363B, 0x00333A29, 0x0033413C,
  0x00334320, 0x0033463B, 0x00334A29, 0x0033510A, 0x00335320, 0x00335824, 0x0033610A, 0x00336321,
  0x00336F12, 0x00337623, 0x00341139, 0x0034153B, 0x00342108, 0x00342409, 0x00342610, 0x00343108,
  0x00343409, 0x00343610, 0x00344108, 0x0034440D, 0x00344610, 0x0034510A, 0x00345309, 0x0034553B,
  0x0034610A, 0x00346309, 0x0034F824, 0x00350029, 0x00352510, 0x00353120, 0x0035330D, 0x00353510,
  0x00354120, 0x0035430D, 0x00354510, 0x00354F28, 0x0035530D, 0x00355510, 0x00355F1F, 0x00356410,
  0x00359626, 0x0035FF12, 0x00366426, 0x0036FF12, 0x0037F426, 0x0039D712, 0x003BF612, 0x003DF81F,
  0x00410004, 0x00411207, 0x0041150D, 0x0041212A, 0x00412420, 0x0041311E, 0x00413308, 0x00413509,
  0x00413F2B, 0x00414208, 0x00420007, 0x0042123C, 0x00421409, 0x00422107, 0x0042223C, 0x00422409,
  0x00422610, 0x0042313C, 0x00423409, 0x0042363B, 0x0042413C, 0x00424320, 0x0042463B, 0x00425108,
  0x00425409, 0x0042FF29, 0x00431107, 0x00431320, 0x0043153B, 0x0043213C, 0x00432320, 0x00432610,
  0x0043313C, 0x00433320, 0x0043353B, 0x00433813, 0x00434108, 0x00434409, 0x00434610, 0x00435108,
  0x0043553B, 0x00435F25, 0x00436309, 0x0043753B, 0x0043FF29, 0x00441239, 0x0044143B, 0x00442139,
  0x00442309, 0x0044253B, 0x00443108, 0x00443220, 0x0044353B, 0x0044410A, 0x00444309, 0x0044453B,
  0x00444813, 0x0044510A, 0x00445309, 0x00445510, 0x00445F25, 0x0044630D, 0x00450026, 0x00452713,
  0x00453120, 0x0045330D, 0x00453510, 0x00454120, 0x0045430D, 0x00454510, 0x00455120, 0x0045530D,
  0x00456209, 0x00456410, 0x0045FF12, 0x00466513, 0x0047FF22, 0x0048FF25, 0x0049F43D, 0x004BFB25,
  0x004EF825, 0x004FFF18, 0x00511339, 0x00512107, 0x00513409, 0x00520007, 0x00521107, 0x00521320,
  0x00522107, 0x00522409, 0x0052313C, 0x00523320, 0x0052353B, 0x00524108, 0x00524320, 0x00531139,
  0x00531309, 0x00532139, 0x00532309, 0x0053253B, 0x00533108, 0x0053340D, 0x00533713, 0x00534108,
  0x0053453B, 0x00534F2B, 0x00535309, 0x00535610, 0x00535F25, 0x0053643B, 0x00541139, 0x00542139,
  0x00542309, 0x00542613, 0x00543139, 0x00543309, 0x00543510, 0x00543F2B, 0x00544309, 0x00544510,
  0x00544F28, 0x0054530D, 0x0054FF12, 0x00553613, 0x00553F2B, 0x00554410, 0x0055510A, 0x0055543B,
  0x00555F25, 0x0055633B, 0x0055FF12, 0x00566513, 0x00577413, 0x0059FF28, 0x005CC33D, 0x005EFB28,
  0x005FFF18, 0x00611339, 0x00612107, 0x00613320, 0x0061A724, 0x00621107, 0x0062140B, 0x00622107,
  0x00622320, 0x00623139, 0x00623320, 0x00631139, 0x0063130C, 0x00632139, 0x00632309, 0x00633139,
  0x00633309, 0x00633626, 0x00633F2B, 0x00634309, 0x00634F2B, 0x0063543B, 0x0063FF12, 0x0064343B,
  0x00643F2B, 0x0064443B, 0x00645209, 0x00665513, 0x0066610A, 0x00666526, 0x0067A616, 0x0069843D,
  0x006CF612, 0x006EF326, 0x006FFF18, 0x0071130C, 0x00721107, 0x00722239, 0x0072291C, 0x0072340B,
  0x00731139, 0x00732239, 0x0073630B, 0x0073FF12, 0x0074430B, 0x00755426, 0x00776F28, 0x00777410,
  0x0078843D, 0x007CF416, 0x007EF326, 0x007FFF18, 0x00822239, 0x00831139, 0x0083430B, 0x0084530B,
  0x0087561C, 0x00887F25, 0x00888426, 0x008AF61C, 0x008F0018, 0x008FFF18, 0x00911107, 0x0093230B,
  0x0094530B, 0x0097743D, 0x00998C25, 0x00999616, 0x009EF825, 0x009FFF18, 0x00A3430B, 0x00A4530B,
  0x00A7743D, 0x00AA9F2B, 0x00AAA616, 0x00ABD61F, 0x00AFFF18, 0x00B3330B, 0x00B44426, 0x00B7643D,
  0x00BB971F, 0x00BBB53D, 0x00BEF512, 0x00BFFF18, 0x00C22139, 0x00C5330E, 0x00C7633D, 0x00CCAF2E,
  0x00CCC616, 0x00CFFF18, 0x00D4440E, 0x00D6420E, 0x00DDCF2E, 0x00DDD516, 0x00DFFF18, 0x00E4330E,
  0x00E6841C, 0x00EEE61C, 0x00EFFF18, 0x00F3320E, 0x00F55319, 0x00F8F41C, 0x00FAFF2E, 0x00FF002E,
  0x00FFF10C, 0x00FFF33D, 0x00FFF722, 0x00FFFF18, 0x01000232, 0x0111113E, 0x01112103, 0x0111311A,
  0x0112111A, 0x01122130, 0x01123130, 0x0112411D, 0x01131102, 0x01132102, 0x01133102, 0x01141108,
  0x01142136, 0x01143136, 0x01144135, 0x0115223B, 0x01211103, 0x0121211A, 0x01213130, 0x01221130,
  0x01222130, 0x01223102, 0x01231104, 0x01232104, 0x01233104, 0x01241139, 0x01241220, 0x01242220,
  0x01251109, 0x0125223B, 0x0125810A, 0x01283212, 0x0131111A, 0x01312130, 0x0131222C, 0x0131322A,
  0x0132122A, 0x0132222D, 0x0132322D, 0x01331207, 0x01332234, 0x01333234, 0x01341139, 0x01343134,
  0x01344134, 0x01348134, 0x0135220B, 0x0136110B, 0x01365224, 0x01411102, 0x01412104, 0x01431239,
  0x01432239, 0x0143320A, 0x01435134, 0x01443107, 0x01444134, 0x01446134, 0x0145220E, 0x01455134,
  0x0147110E, 0x01511102, 0x01521239, 0x01531239, 0x01532239, 0x01533107, 0x0155220E, 0x01555134,
  0x0157110E, 0x01611107, 0x01621239, 0x01631239, 0x01661139, 0x01666134, 0x01711107, 0x01721239,
  0x01745107, 0x0177110C, 0x01811107, 0x01821107, 0x0185110C, 0x0188210C, 0x01911107, 0x01933139,
  0x01A11107, 0x01A31139, 0x01F5220E, 0x02000001, 0x02000127, 0x02000427, 0x02000727, 0x02000E2F,
  0x02110000, 0x02111200, 0x02111411, 0x02111827, 0x02111F2F, 0x02112411, 0x02112715, 0x02113200,
  0x02113411, 0x02113715, 0x02114200, 0x02121200, 0x02121301, 0x02121F2F, 0x02122200, 0x02122615,
  0x02122F2F, 0x02123311, 0x02123F2F, 0x02124411, 0x02131211, 0x02132311, 0x02133211, 0x02184415,
  0x02211200, 0x02211311, 0x02211F2F, 0x02212311, 0x02212F2F, 0x02213211, 0x02221201, 0x02221311,
  0x02221F2F, 0x02222311, 0x02222F2F, 0x02223211, 0x02223F2F, 0x02231211, 0x02232211, 0x02232F2F,
  0x02233211, 0x02233F2F, 0x02287515, 0x022DAB17, 0x02311211, 0x02311527, 0x02312211, 0x02321211,
  0x02322211, 0x02322F2F, 0x02323311, 0x02323F2F, 0x02331211, 0x02332211, 0x02332F2F, 0x02333F2F,
  0x0237FF17, 0x02385615, 0x023D9517, 0x02410027, 0x02487827, 0x024E3117, 0x024FFF2F, 0x02598627,
  0x025DFF2F, 0x025FFF2F, 0x02687827, 0x026DFA17, 0x026FFF2F, 0x02796427, 0x027E4217, 0x027FFF2F,
  0x02888727, 0x028EFF2F, 0x028FFF2F, 0x02984327, 0x029F112F, 0x029FFF2F, 0x02A76527, 0x02AEF717,
  0x02AFFF2F, 0x02B7C827, 0x02BEF917, 0x02BFFF2F, 0x02C66527, 0x02CD5517, 0x02CFFF2F, 0x02D63227,
  0x02DDD527, 0x02DFFF2B, 0x02E84717, 0x02EEE327, 0x02EFFF2F, 0x02F54527, 0x02FCF817, 0x02FFEF2B,
  0x02FFFA2F, 0x02FFFE2F, 0x03000127, 0x03000201, 0x03111200, 0x03122115, 0x03123200, 0x03133211,
  0x03211200, 0x03213127, 0x03221200, 0x03345215, 0x04000F17, 0x04122F17, 0x043F6515, 0x043FFF17,
  0x044F5527, 0x044FFF17, 0x045F0017, 0x045FFF17, 0x046F6517, 0x04710027, 0x047F4427, 0x04810027,
  0x048EFA15, 0x048FFF2F, 0x049F4427, 0x049FFF2F, 0x04AEA727, 0x04AFFF2F, 0x04BE9C15, 0x04BFFF2F,
  0x04CE5427, 0x04CFFF2F, 0x04DE3527, 0x04DFFF17, 0x04EE4627, 0x04EFFF17, 0x04FEF327, 0x04FFFF2F,
  0x06000F27, 0x069FFF17, 0x06FFFF17, 0x08110017, 0x08EFFF15, 0xFFFFFF00
};

// ISO/IEC 23003-3, Annex C.3
static const uint16_t arithCumFreqM[64][ARITH_ESCAPE + 1] = { // arith_cf_m
  {  708,  706,  579,  569,  568,  567,  479,  469,  297,  138,   97,   91,   72,   52,   38,   34, 0},
  { 7619, 6917, 6519, 6412, 5514, 5003, 4683, 4563, 3907, 3297, 3125, 3060, 2904, 2718, 2631, 2590, 0},
  { 7263, 4888, 4810, 4803, 1889,  415,  335,  327,  195,   72,   52,   49,   36,   20,   15,   14, 0},
  { 3626, 2197, 2188, 2187,  582,   57,   47,   46,   30,   12,    9,    8,    6,    4,    3,    2, 0},
  { 7806, 5541, 5451, 5441, 2720,  834,  691,  674,  487,  243,  179,  167,  139,   98,   77,   70, 0},
  { 6684, 4101, 4058, 4055, 1748,  426,  368,  364,  322,  257,  235,  232,  228,  222,  217,  215, 0},
  { 9162, 5964, 5831, 5819, 3269,  866,  658,  638,  535,  348,  258,  244,  234,  214,  195,  186, 0},
  {10638, 8491, 8365, 8351, 4418, 2067, 1859, 1834, 1190,  601,  495,  478,  356,  217,  174,  164, 0},
  {13389,10514,10032, 9961, 7166, 3488, 2655, 2524, 2015, 1140,  760,  672,  585,  426,  325,  283, 0},
  {14861,12788,12115,11952, 9987, 6657, 5323, 4984, 4324, 3001, 2205, 1943, 1764, 1394, 1115,  978, 0},
  {12876,10004, 9661, 9610, 7107, 3435, 2711, 2595, 2257, 1508, 1059,  952,  893,  753,  609,  538, 0},
  {15125,13591,13049,12874,11192, 8543, 7406, 7023, 6291, 4922, 4104, 3769, 3465, 2890, 2486, 2275, 0},
  {14574,13106,12731,12638,10453, 7947, 7233, 7037, 6031, 4618, 4081, 3906, 3465, 2802, 2476, 2349, 0},
  {15070,13179,12517,12351,10742, 7657, 6200, 5825, 5264, 3998, 3014, 2662, 2510, 2153, 1799, 1564, 0},
  {15542,14466,14007,13844,12489,10409, 9481, 9132, 8305, 6940, 6193, 5867, 5458, 4743, 4291, 4047, 0},
  {15165,14384,14084,13934,12911,11485,10844,10513,10002, 8993, 8380, 8051, 7711, 7036, 6514, 6233, 0},
  {15642,14279,13625,13393,12348, 9971, 8405, 7858, 7335, 6119, 4918, 4376, 4185, 3719, 3231, 2860, 0},
  {13408,13407,11471,11218,11217,11216, 9473, 9216, 6480, 3689, 2857, 2690, 2256, 1732, 1405, 1302, 0},
  {16098,15584,15191,14931,14514,13578,12703,12103,11830,11172,10475, 9867, 9695, 9281, 8825, 8389, 0},
  {15844,14873,14277,13996,13230,11535,10205, 9543, 9107, 8086, 7085, 6419, 6214, 5713, 5195, 4731, 0},
  {16131,15720,15443,15276,14848,13971,13314,12910,12591,11874,11225,10788,10573,10077, 9585, 9209, 0},
  {16331,16330,12283,11435,11434,11433, 8725, 8049, 6065, 4138, 3187, 2842, 2529, 2171, 1907, 1745, 0},
  {16011,15292,14782,14528,14008,12767,11556,10921,10591, 9759, 8813, 8043, 7855, 7383, 6863, 6282, 0},
  {16380,16379,15159,14610,14609,14608,12859,12111,11046, 9536, 8348, 7713, 7216, 6533, 5964, 5546, 0},
  {16367,16333,16294,16253,16222,16143,16048,15947,15915,15832,15731,15619,15589,15512,15416,15310, 0},
  {15967,15319,14937,14753,14010,12638,11787,11360,10805, 9706, 8934, 8515, 8166, 7456, 6911, 6575, 0},
  { 4906, 3005, 2985, 2984,  875,  102,   83,   81,   47,   17,   12,   11,    8,    5,    4,    3, 0},
  { 7217, 4346, 4269, 4264, 1924,  428,  340,  332,  280,  203,  179,  175,  171,  164,  159,  157, 0},
  {16010,15415,15032,14805,14228,13043,12168,11634,11265,10419, 9645, 9110, 8892, 8378, 7850, 7437, 0},
  { 8573, 5218, 5046, 5032, 2787,  771,  555,  533,  443,  286,  218,  205,  197,  181,  168,  162, 0},
  {11474, 8095, 7822, 7796, 4632, 1443, 1046, 1004,  748,  351,  218,  194,  167,  121,   93,   83, 0},
  {16152,15764,15463,15264,14925,14189,13536,13070,12846,12314,11763,11277,11131,10777,10383,10011, 0},
  {14187,11654,11043,10919, 8498, 4885, 3778, 3552, 2947, 1835, 1283, 1134,  998,  749,  585,  514, 0},
  {14162,11527,10759,10557, 8601, 5417, 4105, 3753, 3286, 2353, 1708, 1473, 1370, 1148,  959,  840, 0},
  {16205,15902,15669,15498,15213,14601,14068,13674,13463,12970,12471,12061,11916,11564,11183,10841, 0},
  {15043,12972,12092,11792,10265, 7446, 5934, 5379, 4883, 3825, 3036, 2647, 2507, 2185, 1901, 1699, 0},
  {15320,13694,12782,12352,11191, 8936, 7433, 6671, 6255, 5366, 4622, 4158, 4020, 3712, 3420, 3198, 0},
  {16255,16020,15768,15600,15416,14963,14440,14006,13875,13534,13137,12697,12602,12364,12084,11781, 0},
  {15627,14503,13906,13622,12557,10527, 9269, 8661, 8117, 6933, 5994, 5474, 5222, 4664, 4166, 3841, 0},
  {16366,16365,14547,14160,14159,14158,11969,11473, 8735, 6147, 4911, 4530, 3865, 3180, 2710, 2473, 0},
  {16257,16038,15871,15754,15536,15071,14673,14390,14230,13842,13452,13136,13021,12745,12434,12154, 0},
  {15855,14971,14338,13939,13239,11782,10585, 9805, 9444, 8623, 7846, 7254, 7079, 6673, 6262, 5923, 0},
  { 9492, 6318, 6197, 6189, 3004,  652,  489,  477,  333,  143,   96,   90,   78,   60,   50,   47, 0},
  {16313,16191,16063,15968,15851,15590,15303,15082,14968,14704,14427,14177,14095,13899,13674,13457, 0},
  { 8485, 5473, 5389, 5383, 2411,  494,  386,  377,  278,  150,  117,  112,  103,   89,   81,   78, 0},
  {10497, 7154, 6959, 6943, 3788, 1004,  734,  709,  517,  238,  152,  138,  120,   90,   72,   66, 0},
  {16317,16226,16127,16040,15955,15762,15547,15345,15277,15111,14922,14723,14671,14546,14396,14239, 0},
  {16382,16381,15858,15540,15539,15538,14704,14168,13768,13092,12452,11925,11683,11268,10841,10460, 0},
  { 5974, 3798, 3758, 3755, 1275,  205,  166,  162,   95,   35,   26,   24,   18,   11,    8,    7, 0},
  { 3532, 2258, 2246, 2244,  731,  135,  118,  115,   87,   45,   36,   34,   29,   21,   17,   16, 0},
  { 7466, 4882, 4821, 4811, 2476,  886,  788,  771,  688,  531,  469,  457,  437,  400,  369,  361, 0},
  { 9580, 5772, 5291, 5216, 3444, 1496, 1025,  928,  806,  578,  433,  384,  366,  331,  296,  273, 0},
  {10692, 7730, 7543, 7521, 4679, 1746, 1391, 1346, 1128,  692,  495,  458,  424,  353,  291,  268, 0},
  {11040, 7132, 6549, 6452, 4377, 1875, 1253, 1130,  958,  631,  431,  370,  346,  296,  253,  227, 0},
  {12687, 9332, 8701, 8585, 6266, 3093, 2182, 2004, 1683, 1072,  712,  608,  559,  458,  373,  323, 0},
  {13429, 9853, 8860, 8584, 6806, 4039, 2862, 2478, 2239, 1764, 1409, 1224, 1178, 1077,  979,  903, 0},
  {14685,12163,11061,10668, 9101, 6345, 4871, 4263, 3908, 3200, 2668, 2368, 2285, 2106, 1942, 1819, 0},
  {13295,11302,10999,10945, 7947, 5036, 4490, 4385, 3391, 2185, 1836, 1757, 1424,  998,  833,  785, 0},
  { 4992, 2993, 2972, 2970, 1269,  575,  552,  549,  530,  505,  497,  495,  493,  489,  486,  485, 0},
  {15419,13862,13104,12819,11429, 8753, 7220, 6651, 6020, 4667, 3663, 3220, 2995, 2511, 2107, 1871, 0},
  {12468, 9263, 8912, 8873, 5758, 2193, 1625, 1556, 1187,  589,  371,  330,  283,  200,  149,  131, 0},
  {15870,15076,14615,14369,13586,12034,10990,10423, 9953, 8908, 8031, 7488, 7233, 6648, 6101, 5712, 0},
  { 1693,  978,  976,  975,  194,   18,   16,   15,   11,    7,    6,    5,    4,    3,    2,    1, 0},
  { 7992, 5218, 5147, 5143, 2152,  366,  282,  276,  173,   59,   38,   35,   27,   16,   11,   10, 0}
};

// ISO/IEC 23003-3, Annex C.4
static const uint16_t arithCumFreqR[3][4] = { // arith_cf_r
  {12571, 10569, 3696, 0},
  {12661,  5700, 3751, 0},
  {10827,  6884, 2929, 0}
};

static const uint8_t arithFastPkIndex[32] = {
  1, 4, 0, 49, 0, 0, 0, 0, 0, 0, 0, 0, 58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 58, 3, 0, 62
};

// static helper functions
static inline uint8_t arithGetPkIndex (const unsigned ctx) // cumul. frequency table index pki = arith_get_pk(c)
{
  if ((ctx & 0xEEEEE) == 0)
  {
    const unsigned tmp = ctx | (ctx >> 6);

    return arithFastPkIndex[(tmp | (tmp >> 9)) & 31];
  }

  int32_t iMax = ARITH_SIZE - 1;
  int32_t iMin = -1;

  do
  {
    const int32_t  i = iMin + ((iMax - iMin) >> 1);
    const uint32_t j = arithHashM[i];
    const uint32_t k = j >> 8;

    if (ctx < k)      iMax = i;
    else if (ctx > k) iMin = i;
    else return  j & UCHAR_MAX;
  }
  while (iMax > iMin + 1);

  return arithLookupM[iMax]; // pki
}

static inline unsigned writeSymbol (OutputStream* const stream, const bool leadingBitIs1, const uint16_t trailingBits)
{
  const uint8_t lowBits = trailingBits & 0x1F;

  stream->write (leadingBitIs1 ? 1 : 0, 1);
  stream->write (leadingBitIs1 ? 0 : UINT_MAX, trailingBits & 0x20);
  stream->write (leadingBitIs1 ? 0 : (1u << lowBits) - 1u, lowBits);

  return trailingBits + 1; // count
}

// private helper functions
unsigned EntropyCoder::arithCodeSymbol (const uint16_t symbol, const uint16_t* table, OutputStream* const stream /*= nullptr*/)
{
  const unsigned  range = m_acHigh + 1 - m_acLow;
  uint16_t high = m_acHigh;
  uint16_t low  = m_acLow;
  unsigned bitCount = 0;

  if (symbol > 0)
  {
    high = low + ((range * table[symbol - 1]) >> 14) - 1;
  }
  low += (range * table[symbol]) >> 14; // NOTE: the spec incorrectly reads "[symbol-1]" here

  if (stream != nullptr) // write-out
  {
    while (true)
    {
      if (high <= SHRT_MAX)
      {
        bitCount += writeSymbol (stream, false, m_acBits);
        m_acBits = 0;
      }
      else if (low > SHRT_MAX)
      {
        bitCount += writeSymbol (stream, true, m_acBits);
        m_acBits = 0;

        high += SHRT_MIN;
        low  += SHRT_MIN;
      }
      else if ((low > (SHRT_MAX >> 1)) && (high < ((-3 * SHRT_MIN) >> 1)))
      {
        m_acBits++;

        high += SHRT_MIN >> 1;
        low  += SHRT_MIN >> 1;
      }
      else break;

      high <<= 1; high |= 1;
      low  <<= 1;
    }
  }
  else // stream == nullptr, counting
  {
    while (true)
    {
      if (high <= SHRT_MAX)
      {
        bitCount += m_acBits + 1;
        m_acBits = 0;
      }
      else if (low > SHRT_MAX)
      {
        bitCount += m_acBits + 1;
        m_acBits = 0;

        high += SHRT_MIN;
        low  += SHRT_MIN;
      }
      else if ((low > (SHRT_MAX >> 1)) && (high < ((-3 * SHRT_MIN) >> 1)))
      {
        m_acBits++;

        high += SHRT_MIN >> 1;
        low  += SHRT_MIN >> 1;
      }
      else break;

      high <<= 1; high |= 1;
      low  <<= 1;
    }
  }
  m_acHigh = high;
  m_acLow  = low;

  return bitCount;
}

unsigned EntropyCoder::arithGetContext (const unsigned ctx, const unsigned idx) // c = arith_get_context(c, i, N)
{
  unsigned c = (ctx & 0xFFFF) >> 4; // NOTE: the "& 0xFFFF" was part of some USAC corrigendum

  c = (c | ((unsigned) m_qcPrev[idx + 1] << 12)) & 0xFFF0;  // add top-left previous neighbor

  if (idx > 0) // lower neighbor(s)
  {
    c |= (unsigned) m_qcCurr[idx - 1];

    if ((idx > 3) && (m_qcCurr[idx - 3] + m_qcCurr[idx - 2] + m_qcCurr[idx - 1] < 5))
    {
      return c | 0x10000;
    }
  }
  return c;  // updated context ctx
}

unsigned EntropyCoder::arithMapContext (const bool arithResetFlag)  // c = arith_map_context(N, arith_reset_flag)
{
  if (arithResetFlag)
  {
    memset (m_qcPrev, 0, m_maxTupleLength * sizeof (uint8_t));
  }
  else if (m_shiftTrafoCurr == m_shiftTrafoPrev)
  {
    memcpy (m_qcPrev, m_qcCurr, m_acSize * sizeof (uint8_t));
  }
  else if (m_shiftTrafoCurr > m_shiftTrafoPrev)
  {
    const int shiftDiff = (int) m_shiftTrafoCurr - m_shiftTrafoPrev;

    for (int i = m_acSize - 1; i >= 0; i--)
    {
      m_qcPrev[i] = m_qcCurr[i << shiftDiff];
    }
  }
  else // (m_shiftTrafoCurr < m_shiftTrafoPrev)
  {
    const int shiftDiff = (int) m_shiftTrafoPrev - m_shiftTrafoCurr;

    for (int i = m_acSize - 1; i >= 0; i--)
    {
      m_qcPrev[i] = m_qcCurr[i >> shiftDiff];
    }
  }
  m_qcPrev[m_acSize] = 0; // for encoding speed

  return (unsigned) m_qcPrev[0] << 12; // initial context ctx with top-left previous neighbor
}

#if EC_TRELLIS_OPT_CODING
void EntropyCoder::arithSetContext (const unsigned newCtxState, const uint16_t sigEnd)
{
  m_csCurr = newCtxState;
  m_acBits = (m_csCurr >> 17) & 31;
  if (sigEnd > 0) m_qcCurr[sigEnd - 1] = (m_csCurr >> 22) & 0xF;
  if (sigEnd > 1) m_qcCurr[sigEnd - 2] = (m_csCurr >> 26) & 0xF;
  if (sigEnd > 2) m_qcCurr[sigEnd - 3] = (m_csCurr >> 30);
}
#endif

// constructor
EntropyCoder::EntropyCoder ()
{
  // initialize all helper buffers
  m_qcCurr = nullptr;
  m_qcPrev = nullptr;

  // initialize encoding variables
  m_acBits = 0;
  m_acHigh = USHRT_MAX;
  m_acLow  = 0;
  m_acSize = 0;
  m_csCurr = 0;
  m_maxTupleLength = 0;
  m_shiftTrafoCurr = 0;
  m_shiftTrafoPrev = 0;
}

// destructor
EntropyCoder::~EntropyCoder ()
{
  // free allocated helper buffers
  MFREE (m_qcCurr);
  MFREE (m_qcPrev);
}

// public functions
unsigned EntropyCoder::arithCodeSigMagn (const uint8_t* const magn, const uint16_t sigOffset, const uint16_t sigLength,
                                         const bool arithFinish /*= false*/, OutputStream* const stream /*= nullptr*/)
{
  const uint8_t* a = &magn[sigOffset    ];
  const uint8_t* b = &magn[sigOffset + 1];
  unsigned c = m_csCurr & 0x1FFFF;
  unsigned bitCount = 0;
  uint16_t r[7];
  uint16_t sigEnd = (sigOffset >> 1) + (sigLength >> 1);

  if (arithFinish && (sigLength > 0)) // try to save bits via signalling of ARITH_STOP symbol
  {
    int i = sigLength - 2;

    while ((i >= 0) && ((a[i] | b[i]) == 0)) i -= 2;
    i = (sigOffset + i + 2) >> 1;

    if (i + 26 < (int) sigEnd) sigEnd = (uint16_t) i;
  }

  for (uint16_t s = sigOffset >> 1; s < sigEnd; s++)
  {
    uint32_t lev = 0;
    uint16_t a1 = *a;
    uint16_t b1 = *b;

    a += 2; b += 2;

    // arith_get_context, cf Scl. 7.4
    c = arithGetContext (c, s);
    // arith_update_context, Scl. 7.4
    m_qcCurr[s] = __min (0xF, a1 + b1 + 1);

    // MSB encoding as in Scl. B.25.3
    while ((a1 > 3) || (b1 > 3))
    {
      // write escaped codeword value
      bitCount += arithCodeSymbol (ARITH_ESCAPE, arithCumFreqM[arithGetPkIndex (c | (lev << 17))], stream);
      // store LSBs in r, right-shift
      r[lev++] = (a1 & 1) | ((b1 & 1) << 1);
      a1 >>= 1; b1 >>= 1;
    }
    // write the m MSB codeword value
    bitCount += arithCodeSymbol (a1 | (b1 << 2), arithCumFreqM[arithGetPkIndex (c | (lev << 17))], stream);

    // LSB encoding, Table 38, B.25.3
    while (lev--)
    {
      const uint16_t rLev = r[lev];

      bitCount += arithCodeSymbol (rLev, arithCumFreqR[a1 == 0 ? 1 : (b1 == 0 ? 0 : 2)], stream);
      a1 = (a1 << 1) | (rLev & 1);
      b1 = (b1 << 1) | ((rLev >> 1) & 1);
    }
  } // for s

  if (arithFinish) // flush last bits
  {
    // NOTE: the spec incorrectly reads "m_acBits" below (bits_to_follow++ missing in B.25.3)
    if (sigLength > 0)
    {
      if (sigEnd < (sigOffset >> 1) + (sigLength >> 1)) // write ARITH_STOP flag to save bits
      {
        c = arithGetContext (c, sigEnd);
        bitCount += arithCodeSymbol (ARITH_ESCAPE, arithCumFreqM[arithGetPkIndex (c)], stream);
        bitCount += arithCodeSymbol (0 /*m=STOP*/, arithCumFreqM[arithGetPkIndex (c | (1 << 17))], stream);
      }
      bitCount += writeSymbol (stream, m_acLow > (SHRT_MAX >> 1), m_acBits + 1);
    }
    m_csCurr = m_acBits = 0;
  }
  else
  {
    m_csCurr = 0;
    if (sigEnd > 0) m_csCurr |= m_qcCurr[sigEnd - 1] << 22;
    if (sigEnd > 1) m_csCurr |= m_qcCurr[sigEnd - 2] << 26;
    if (sigEnd > 2) m_csCurr |= __min (3, m_qcCurr[sigEnd - 3]) << 30;
  }
  m_csCurr |= ((unsigned) __min (31, m_acBits) << 17) | c;

  return bitCount;
}

#if EC_TRELLIS_OPT_CODING
unsigned EntropyCoder::arithCodeSigTest (const uint8_t* const magn, const uint16_t sigOffset, const uint16_t sigLength)
{
  const unsigned inAcBits = m_acBits;
  const uint8_t* a = &magn[sigOffset    ];
  const uint8_t* b = &magn[sigOffset + 1];
  unsigned c = m_csCurr & 0x1FFFF;
  unsigned bitCount = 0;
  uint16_t r[7];
  int16_t s = sigOffset >> 1;

  for (uint16_t sigEnd = (uint16_t) s + (sigLength >> 1); s < sigEnd; s++)
  {
    uint32_t lev = 0;
    uint16_t a1 = *a;
    uint16_t b1 = *b;

    a += 2; b += 2;

    // arith_get_context, cf Scl. 7.4
    c = arithGetContext (c, (unsigned) s);
    // arith_update_context, Scl. 7.4
    m_qcCurr[s] = __min (0xF, a1 + b1 + 1);

    // MSB encoding as in Scl. B.25.3
    while ((a1 > 3) || (b1 > 3))
    {
      // write escaped codeword value
      bitCount += arithCodeSymbol (ARITH_ESCAPE, arithCumFreqM[arithGetPkIndex (c | (lev << 17))]);
      // store LSBs in r, right-shift
      r[lev++] = (a1 & 1) | ((b1 & 1) << 1);
      a1 >>= 1; b1 >>= 1;
    }
    // write the m MSB codeword value
    bitCount += arithCodeSymbol (a1 | (b1 << 2), arithCumFreqM[arithGetPkIndex (c | (lev << 17))]);

    // LSB encoding, Table 38, B.25.3
    while (lev--)
    {
      const uint16_t rLev = r[lev];

      bitCount += arithCodeSymbol (rLev, arithCumFreqR[a1 == 0 ? 1 : (b1 == 0 ? 0 : 2)]);
      a1 = (a1 << 1) | (rLev & 1);
      b1 = (b1 << 1) | ((rLev >> 1) & 1);
    }
  } // for s

  m_csCurr = m_qcCurr[--s] << 22;
  if ((s--) > 0) m_csCurr |= m_qcCurr[s] << 26;
  if ((s--) > 0) m_csCurr |= __min (3, m_qcCurr[s]) << 30;
  m_csCurr |= ((unsigned) __min (31, m_acBits) << 17) | c;
  bitCount += m_acBits;

  return bitCount - __min (inAcBits, bitCount);
}

unsigned EntropyCoder::arithCodeTupTest (const uint8_t* const magn, const uint16_t sigOffset)
{
  const unsigned inAcBits = m_acBits;
  uint16_t a1 = magn[sigOffset    ];
  uint16_t b1 = magn[sigOffset + 1];
  unsigned bitCount = 0;
  uint32_t lev = 0;
  uint16_t r[7];
  int16_t  s = sigOffset >> 1;

  // arith_get_context, cf Scl. 7.4
  const unsigned c = arithGetContext (m_csCurr & 0x1FFFF, (unsigned) s);
  // arith_update_context, Scl. 7.4
  m_qcCurr[s] = __min (0xF, a1 + b1 + 1);

  // MSB encoding as in Scl. B.25.3
  while ((a1 > 3) || (b1 > 3))
  {
    // write escaped codeword value
    bitCount += arithCodeSymbol (ARITH_ESCAPE, arithCumFreqM[arithGetPkIndex (c | (lev << 17))]);
    // store LSBs in r, right-shift
    r[lev++] = (a1 & 1) | ((b1 & 1) << 1);
    a1 >>= 1; b1 >>= 1;
  }
  // write the m MSB codeword value
  bitCount += arithCodeSymbol (a1 | (b1 << 2), arithCumFreqM[arithGetPkIndex (c | (lev << 17))]);

  // LSB encoding, Table 38, B.25.3
  while (lev--)
  {
    const uint16_t rLev = r[lev];

    bitCount += arithCodeSymbol (rLev, arithCumFreqR[a1 == 0 ? 1 : (b1 == 0 ? 0 : 2)]);
    a1 = (a1 << 1) | (rLev & 1);
    b1 = (b1 << 1) | ((rLev >> 1) & 1);
  }

  m_csCurr = m_qcCurr[s] << 22;
  if ((s--) > 0) m_csCurr |= m_qcCurr[s] << 26;
  if ((s--) > 0) m_csCurr |= __min (3, m_qcCurr[s]) << 30;
  m_csCurr |= ((unsigned) __min (31, m_acBits) << 17) | c;
  bitCount += m_acBits;

  return bitCount - __min (inAcBits, bitCount);
}
#endif // EC_TRELLIS_OPT_CODING

unsigned EntropyCoder::arithGetResetBit (const uint8_t* const magn, const uint16_t sigOffset, const uint16_t sigLength)
{
  const uint16_t sigEnd = (sigOffset >> 1) + (sigLength >> 1);
  const uint8_t* a = &magn[sigOffset    ];
  const uint8_t* b = &magn[sigOffset + 1];
  unsigned  qcDiff = 0;

  for (uint16_t s = sigOffset >> 1; s < sigEnd; s++)
  {
    const int qcCurrS = __min (0xF, *a + (int) *b + 1);
    const int qcDiffS = qcCurrS - m_qcPrev[s];

    qcDiff += qcDiffS * qcDiffS;
    a += 2; b += 2;
  }

  return (qcDiff >= sigLength * 4u ? 1 : 0); // reset encoder if difference exceeds threshold
}

unsigned EntropyCoder::indexGetBitCount (const int scaleFactorDelta) const
{
  return huffScf[CLIP_PM (scaleFactorDelta, INDEX_OFFSET) + INDEX_OFFSET] & UCHAR_MAX;
}

unsigned EntropyCoder::indexGetHuffCode (const int scaleFactorDelta) const
{
  return huffScf[CLIP_PM (scaleFactorDelta, INDEX_OFFSET) + INDEX_OFFSET] >> 8;
}

unsigned EntropyCoder::initCodingMemory (const unsigned maxTransfLength)
{
  const unsigned max2TupleLength = maxTransfLength >> 1; // tuple buffer size, maxWinLength/4

  if ((maxTransfLength < 128) || (maxTransfLength > 8192) || (maxTransfLength & 7))
  {
    return 1; // invalid arguments error
  }

  m_maxTupleLength = max2TupleLength;
  MFREE (m_qcCurr);
  MFREE (m_qcPrev);

  if ((m_qcCurr = (uint8_t*) malloc (max2TupleLength * sizeof (uint8_t))) == nullptr ||
      (m_qcPrev = (uint8_t*) malloc ((max2TupleLength + 1) * sizeof (uint8_t))) == nullptr)
  {
    return 2; // memory allocation error
  }

  memset (m_qcCurr, 0, max2TupleLength * sizeof (uint8_t));

  return 0; // no error
}

unsigned EntropyCoder::initWindowCoding (const bool forceArithReset, const uint8_t winLenShift /*= 0*/)
{
  // arith_first_symbol() in Scl. B.25.3
  m_acBits = 0;
  m_acHigh = USHRT_MAX;
  m_acLow  = 0;
  m_acSize = m_maxTupleLength >> __min (3u, winLenShift);

  m_shiftTrafoPrev = m_shiftTrafoCurr;
  m_shiftTrafoCurr = winLenShift;

  m_csCurr = arithMapContext (forceArithReset); // m_qcPrev
  memset (m_qcCurr, 1, m_acSize * sizeof (uint8_t)); // reset m_qcCurr, see also arith_finish

  return 0; // no error
}
