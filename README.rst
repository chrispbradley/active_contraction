==================
active_contraction
==================

Active contraction 

Building the example
====================

Instructions on how to configure and build with CMake::

  git clone https://github.com/OpenCMISS-Examples/active_contraction.git
  mkdir build
  cmake -DOpenCMISSLibs_DIR=/path/to/opencmisslib/install ../active_contraction
  make  # cmake --build . will also work here and is much more platform agnostic.

Running the example
===================

Explain how the example is run::

  cd build
  ./src/fortran/active_contraction.F90

or maybe it is a Python only example::

  source /path/to/opencmisslibs/install/virtaul_environments/oclibs_venv_pyXY_release/bin/activate
  python src/python/active_contraction.py

where the XY in the path are the Python major and minor versions respectively.

Prerequisites
=============

None

License
=======

Apache 2.0 License
