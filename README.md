exhale edit
======

This fork of exhale has custom experimental features do not use for archiving!.

First experimental feature is tonality encoding, dynamic bandwidth controls (by bitrate), and minimum SNR.

exhale, which is an acronym for "Ecodis eXtended High-efficiency And
Low-complexity Encoder", is a lightweight library and application to
encode uncompressed WAVE-format audio files into MPEG-4-format files
complying with the ISO/IEC 23003-3 (MPEG-D) Unified Speech and Audio
Coding (USAC, also known as Extended High-Efficiency AAC) standard.

In addition, exhale writes program peak-level and loudness data into
the generated MPEG-4 files according to the ISO/IEC 23003-4, Dynamic
Range Control (DRC) specification for use by decoders providing DRC.

exhale currently makes use of all frequency-domain (FD) coding tools
in the scalefactor based MDCT processing path. Its objective is high
quality mono, stereo, and multichannel coding at medium and high bit
rates, so the lower-rate USAC coding tools (ACELP, TCX, Enhanced SBR
and MPEG Surround with Unified Stereo coding) won't be integrated.

**Important:** Due to the missing lower-rate coding tools, the audio
quality at the lowest of exhale's bit-rate modes (18 kbit/s mono, 36
kbit/s stereo) doesn't reflect the full capabilities of the Extended
HE-AAC standard. Therefore, use the lowest bit-rate modes *only when
required*. Also, please don't attempt to modify exhale's source code
or to configure the command-line encoder to produce lower bit-rates.
Use only existing presets and input sampling rates of 32...48 kHz.

____________________________________________________________________


Copyright
---------

(c) 2025 Christian R. Helmrich, project ecodis. All rights reserved.


License
-------

exhale is being made available under an open-source license which is
based on the 3-clause BSD license but modified to address particular
aspects dictated by the nature and the output of this application.

The license text and release notes for the current version 1.2.1 can
be found in the `include` subdirectory of the exhale distribution.


Compilation
-----------

This section describes how to compile the exhale source code into an
executable application under Linux and Microsoft Windows. The binary
application files will show up in a newly created `bin` subdirectory
of the exhale distribution directory and/or a subdirectory thereof.

Note that, for advanced use cases, cmake files are provided as well.
See https://gitlab.com/ecodis/exhale/-/merge_requests/2 for details.

### Linux and MacOS (GNU Compiler Collection, gcc):

In a terminal, change to the exhale distribution directory and enter

`
make release
`

to build a release-mode executable with the default (usually 64-bit)
configuration. A 32-bit debug-mode executable can be built by typing

`
make BUILD32=1 debug
`

### Microsoft Windows (Visual Studio 2012 and later):

Doubleclick the exhale_vs2012.sln file to open the project in Visual
Studio. Once it's loaded, rightclick on `exhaleApp` in the "Solution
Explorer" window on the right-hand side, then select `Set as StartUp
Project`. Now simply press `F7` to build the solution in debug mode.

To change the debugging command, rightclick again on `exhaleApp` and
select `Properties`. In the newly opened window click on `Debugging`
under "Configuration Properties" on the left-hand side. Then you can
edit the "Command Arguments" entry on the right-hand side as needed.

For fastest encoding speed, please select `Release` and `x64` before
building the solution. This will create a release-mode 64-bit binary.
If you would like to build a dynamically linked library (DLL) of the
exhale source instead of an application binary, select `Release DLL`
instead of `Release`, rightclick on `exhaleLib`, and select `Build`.


Usage
-----

This section describes how to run the exhale application either from
the command-line or using a third-party software providing WAVE data
to exhale's standard input pipe (stdin), such as foobar2000. See the
Wiki at https://gitlab.com/ecodis/exhale/-/wikis/home for more info.

### Standalone (command-line):

In a terminal, change to exhale's `bin` subdirectory and then enter

`./exhale` (on Linux and MacOS) or `exhale.exe` (on Windows)

to print out usage information. As an example, the following command

`exhale.exe 5 C:\Music\Input.wav C:\Music\Output.m4a`

converts file Input.wav to file Output.m4a at roughly 128 kbit/s (if
the input signal is 2-channel stereo) and in Extended HE-AAC format.

There is also an **expert mode** providing two additional arguments:

`exhale.exe b s 42 C:\Music\Input.wav C:\Music\Output.m4a`

e.g. encodes Input.wav to Output.m4a at roughly 48 kbit/s stereo and
with SBR enabled, seamless operation (`s` forces media time 0 in the
edit list), and an independent frame interval of 42 (range 10...99).

### Third-party stdin (foobar2000):

After downloading from www.foobar2000.org and starting the software,
load the desired input audio files into the playlist. Mark all files
to be converted, rightclick on one of them, and select `Convert` ->
`...`. In the newly opened window click on `Output format`. Once the
window content changed, double-click on entry `AAC (exhale)` and set
up the conversion. If that entry does not exist, click on `Add New`,
select `Custom` under "Encoder" and enter the following information:

- *Encoder file:* exhale.exe (including path to the executable)
- *Extension:* m4a
- *Parameters:* # %d (where # is the bit-rate mode, i.e., 0...9 when
                      SBR is disabled, or a...g when SBR is enabled)
- *Format is:* lossy
- *Highest BPS mode supported:* 24 (or 32, doesn't matter much)
- *Encoder name:* Extended HE-AAC (exhale)
- *Bitrate (kbps):* (depends on bit-rate mode, see Usage above)
- *Settings:* CVBR mode # (where # equals that in *Parameters*)

Then click on `OK` and on `Back` and, in the first "Converter Setup"
window, on `Other` and ensure the "Transfer..." box for the class of
input metadata that you wish to copy to the output files is checked.
Now set the destination settings as desired and click on `Convert`.


Development
-----------

If you are interested in contributing to exhale, please email one of
the developers. Merge requests with fixes and/or speedups are highly
appreciated.



Software License for The Fraunhofer FDK AAC Codec Library for Android

© Copyright 1995 - 2012 Fraunhofer-Gesellschaft zur Förderung der angewandten Forschung e.V.

 All rights reserved.
1. INTRODUCTION The Fraunhofer FDK AAC Codec Library for Android ("FDK AAC Codec") is software that implements the MPEG Advanced Audio Coding ("AAC") encoding and decoding scheme for digital audio. This FDK AAC Codec software is intended to be used on a wide variety of Android devices.

AAC's HE-AAC and HE-AAC v2 versions are regarded as today's most efficient general perceptual audio codecs. AAC-ELD is considered the best-performing full-bandwidth communications codec by independent studies and is widely deployed. AAC has been standardized by ISO and IEC as part of the MPEG specifications.

Patent licenses for necessary patent claims for the FDK AAC Codec (including those of Fraunhofer) may be obtained through Via Licensing (www.vialicensing.com) or through the respective patent owners individually for the purpose of encoding or decoding bit streams in products that are compliant with the ISO/IEC MPEG audio standards. Please note that most manufacturers of Android devices already license these patent claims through Via Licensing or directly from the patent owners, and therefore FDK AAC Codec software may already be covered under those patent licenses when it is used for those licensed purposes only.

Commercially-licensed AAC software libraries, including floating-point versions with enhanced sound quality, are also available from Fraunhofer. Users are encouraged to check the Fraunhofer website for additional applications information and documentation.

2. COPYRIGHT LICENSE

Redistribution and use in source and binary forms, with or without modification, are permitted without payment of copyright license fees provided that you satisfy the following conditions:

You must retain the complete text of this software license in redistributions of the FDK AAC Codec or your modifications thereto in source code form.

You must retain the complete text of this software license in the documentation and/or other materials provided with redistributions of the FDK AAC Codec or your modifications thereto in binary form. You must make available free of charge copies of the complete source code of the FDK AAC Codec and your modifications thereto to recipients of copies in binary form.

The name of Fraunhofer may not be used to endorse or promote products derived from this library without prior written permission.

You may not charge copyright license fees for anyone to use, copy or distribute the FDK AAC Codec software or your modifications thereto.

Your modified versions of the FDK AAC Codec must carry prominent notices stating that you changed the software and the date of any change. For modified versions of the FDK AAC Codec, the term "Fraunhofer FDK AAC Codec Library for Android" must be replaced by the term "Third-Party Modified Version of the Fraunhofer FDK AAC Codec Library for Android."

3. NO PATENT LICENSE

NO EXPRESS OR IMPLIED LICENSES TO ANY PATENT CLAIMS, including without limitation the patents of Fraunhofer, ARE GRANTED BY THIS SOFTWARE LICENSE. Fraunhofer provides no warranty of patent non-infringement with respect to this software.

You may use this FDK AAC Codec software or modifications thereto only for purposes that are authorized by appropriate patent licenses.

4. DISCLAIMER

This FDK AAC Codec software is provided by Fraunhofer on behalf of the copyright holders and contributors "AS IS" and WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES, including but not limited to the implied warranties of merchantability and fitness for a particular purpose. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE for any direct, indirect, incidental, special, exemplary, or consequential damages, including but not limited to procurement of substitute goods or services; loss of use, data, or profits, or business interruption, however caused and on any theory of liability, whether in contract, strict liability, or tort (including negligence), arising in any way out of the use of this software, even if advised of the possibility of such damage.

5. CONTACT INFORMATION

Fraunhofer Institute for Integrated Circuits IIS Attention: Audio and Multimedia Departments - FDK AAC LL Am Wolfsmantel 33 91058 Erlangen, Germany

www.iis.fraunhofer.de/amm amm-info@iis.fraunhofer.de


