==================
Active Contraction
==================

Active contraction 

Building the example
====================

This example can be configure and built with CMake using the following commands::

  git clone https://github.com/OpenCMISS-Examples/active_contraction.git
  cd active_contraction
  mkdir build
  cd build
  cmake -DOpenCMISS_INSTALL_ROOT=/path/to/opencmiss/install ../.
  make  # cmake --build . will also work here and is much more platform agnostic.

Running the example
===================

To run the example execute the following commands (assumes current directory is the build directory from 'Building the example')::

  cd src/fortran
  ./active_contraction


Prerequisites
=============

None

License
=======

Apache 2.0 License
