import cython
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector


cdef extern from "GroundTruthStatistics.h" namespace "GroundTruthStatistics":

    cdef cppclass GroundTruthStats:
        GroundTruthStats()
        GroundTruthStats(int _nTotalGts, int * _gtNrs, double * _gtXs, double * _gtYs, double * _gtZs,
            double * _gtPhotons, int _nrGroundTruthMolecules, double _tolXY, double _tolZ) except +
        void AddLocations(int m, int *nrs, double *xs, double *ys, double *zs, double *photons)
        void RemoveLocations(int m, int *nrs, double *xs, double *ys, double *zs, double *photons);
        MatchingResult MatchFrame(int frameNr, int m, int *nrs, double *xs, double *ys, double *zs, double *photons)
        Stats GetStat()
        void Clear()
    
    cdef cppclass Stats:
        int nTP
        int nFP
        int nFN
        double rmsXY
        double rmsZ
        double jac
        double recall
        double precision
        double averageDx
        double averageDy
        double averageDz
        double avgPhotonsScale
        double avgPhotonsDelta
        
        void Print()
        
    cdef cppclass MatchingResult:
        MatchingResult()
        vector[int] idsTest
        vector[int] idsTrue
        vector[double] distsXY
        vector[double] distsZ
    
cdef class PyGroundTruthStats:
    cdef GroundTruthStats _cppStats
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __cinit__(self, np.ndarray[int, ndim = 1, mode = "c"] nrs not None,
                  np.ndarray[double, ndim = 1, mode = "c"] xs not None,
                  np.ndarray[double, ndim = 1, mode = "c"] ys not None,
                  np.ndarray[double, ndim = 1, mode = "c"] zs not None,
                  np.ndarray[double, ndim = 1, mode = "c"] photons not None,
                  int nrGroundTruthMolecules, double tolXY, double tolZ):
        
        cdef int _nTotalGts = nrs.shape[0]
        self._cppStats = GroundTruthStats(_nTotalGts, & nrs[0], & xs[0], \
                                          & ys[0], & zs[0], & photons[0], \
                                          nrGroundTruthMolecules, tolXY, tolZ)
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def AddLocations(self, np.ndarray[int, ndim = 1, mode = "c"] nrs not None,
                  np.ndarray[double, ndim = 1, mode = "c"] xs not None,
                  np.ndarray[double, ndim = 1, mode = "c"] ys not None,
                  np.ndarray[double, ndim = 1, mode = "c"] zs not None,
                  np.ndarray[double, ndim = 1, mode = "c"] photons not None,):
        
        cdef int n = nrs.shape[0]
        self._cppStats.AddLocations(n, &nrs[0], &xs[0], &ys[0], &zs[0], &photons[0])
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def RemoveLocations(self, np.ndarray[int, ndim = 1, mode = "c"] nrs not None,
                  np.ndarray[double, ndim = 1, mode = "c"] xs not None,
                  np.ndarray[double, ndim = 1, mode = "c"] ys not None,
                  np.ndarray[double, ndim = 1, mode = "c"] zs not None,
                  np.ndarray[double, ndim = 1, mode = "c"] photons not None,):
        
        cdef int n = nrs.shape[0]
        self._cppStats.RemoveLocations(n, &nrs[0], &xs[0], &ys[0], &zs[0], &photons[0])

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def MatchFrame(self, int frameNr,
                   np.ndarray[int, ndim = 1, mode = "c"] nrs not None,
                   np.ndarray[double, ndim = 1, mode = "c"] xs not None,
                   np.ndarray[double, ndim = 1, mode = "c"] ys not None,
                   np.ndarray[double, ndim = 1, mode = "c"] zs not None,
                   np.ndarray[double, ndim = 1, mode = "c"] photons not None,):
        
        cdef int n = nrs.shape[0]
        cdef MatchingResult res = self._cppStats.MatchFrame(frameNr, n, &nrs[0], &xs[0], \
                                                            &ys[0], &zs[0], &photons[0])
        return res.idsTest, res.idsTrue, res.distsXY, res.distsZ
        
    def Clear(self):
        return self._cppStats.Clear()
    
    def GetStat(self):
        cdef Stats stats = self._cppStats.GetStat()
        res = {"nFP": stats.nFP,
               "nFN": stats.nFN,
               "rmsXY": stats.rmsXY,
               "rmsZ": stats.rmsZ,
               "jac": stats.jac,
               "recall": stats.recall,
               "precision": stats.precision,
               "averageDx": stats.averageDx,
               "averageDy": stats.averageDy,
               "averageDz": stats.averageDz,
               "avgPhotonsScale": stats.avgPhotonsScale,
               "avgPhotonsDelta": stats.avgPhotonsDelta}
        return res