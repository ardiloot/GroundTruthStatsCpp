from setuptools import setup, Extension
from Cython.Distutils import build_ext
import numpy

cmdclass = {'build_ext': build_ext}

ext = Extension("GroundTruthStatsCpp._statsCpp",
                sources = ["GroundTruthStatsCpp/src/GroundTruthStatCpp.pyx",
                           "GroundTruthStatsCpp/src/adapter.cpp",
                           "GroundTruthStatsCpp/src/GroundTruthStatistics.cpp",
                           "GroundTruthStatsCpp/src/munkres.cpp",
                           "GroundTruthStatsCpp/src/std2dvectordapter.cpp"],
                include_dirs = [numpy.get_include()],
                language = "c++")

setup(name = "GroundTruthStatsCpp",
      author = "Ardi Loot",
      author_email = "ardi.loot@outlook.com",
      packages = ["GroundTruthStatsCpp"],
      cmdclass = cmdclass,
      ext_modules = [ext],)

