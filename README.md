# Albacore Southern California offshore Velocity Model (ALBACORE)

<a href="https://github.com/sceccode/albacore.git"><img src="https://github.com/sceccode/albacore/wiki/images/albacore_logo.png"></a>

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![GitHub repo size](https://img.shields.io/github/repo-size/sceccode/albacore)
[![albacore-ucvm-ci Actions Status](https://github.com/SCECcode/albacore/workflows/albacore-ucvm-ci/badge.svg)](https://github.com/SCECcode/albacore/actions)

Albacore Southern California offshore Velocity Model

ALBACORE is a shear wave velocity model for offshore Southern California that images plate boundary deformation including both thickening and thinning of the crustal and mantle lithosphere at the westernmost edge of the North American continent. The Asthenospheric and Lithospheric Broadband Architecture from the California Offshore Region Experiment (ALBACORE) ocean bottom seismometer array, together with 65 stations of the onshore Southern California Seismic Network, was used to measure ambient noise correlation functions and Rayleigh wave dispersion curves which are inverted for 3-D shear wave velocities. The resulting velocity model defines the transition from continental lithosphere to oceanic, illuminating the complex history and deformation in the region.

Bowden, D. C., Kohler, M. D., Tsai, V. C., & Weeraratne, D. S. (2016). Offshore southern California lithospheric velocity structure from noise cross-correlation functions. Journal of Geophysical Research, 121(5), 3415-3427. doi:10.1002/2016JB012919

## Installation

This package is intended to be installed as part of the UCVM framework,
version 25.7 or higher. 

## Library

The library ./lib/libalbacore.a may be statically linked into any
user application. Also, if your system supports dynamic linking,
you will also have a ./lib/libalbacore.so file that can be used
for dynamic linking. The header file defining the API is located
in ./include/albacore.h.

## Contact the authors

If you would like to contact the authors regarding this software,
please e-mail software@scec.org. Note this e-mail address should
be used for questions regarding the software itself (e.g. how
do I link the library properly?). Questions regarding the model's
science (e.g. on what paper is the ALBACORE based?) should be directed
to the model's authors, located in the AUTHORS file.
