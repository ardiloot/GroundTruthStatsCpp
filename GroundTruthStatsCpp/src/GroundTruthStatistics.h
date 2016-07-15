#ifndef _GROUND_TRUTH_STATISTICS_H_
#define _GROUND_TRUTH_STATISTICS_H_

#include <iostream>
#include <set>
#include <vector>
#include <map>
#include "munkres.h"
#include "std2dvectordapter.h"

namespace GroundTruthStatistics {

	// Methods
	double sqr(double x);

	// Locations
	class Locations {
	public:
		std::vector<int> nrs, indices;
		std::vector<double> xs, ys, zs, photons;

		Locations();
		void Add(int nr, double x, double y, double z, double photon, int index = 0);
		void Remove(int nr, double x, double y, double z, double photon);
		int IndexOf(int nr, double x, double y, double z, double photon);
		int size();
	private:

	};

	// MatchingResult
	class MatchingResult {
	public:
		int frameNr;
		std::vector<int> idsTest, idsTrue;
		std::vector<double> distsXY, distsZ;

		MatchingResult(int _frameNr = -1);
		void Add(int idTest, int idTrue, double distXY, double distZ);
		int size();

	private:
	};

	// Stats
	class Stats {
	public:
		int nTP;
		int nFP;
		int nFN;
		double rmsXY;
		double rmsZ;
		double jac;
		double recall;
		double precision;
		double averageDx;
		double averageDy;
		double averageDz;
		double avgPhotonsScale;
		double avgPhotonsDelta;

		Stats () {};
		Stats(int nTruePos, int nFalsePos, int nFalseNeg, int nrGroundTruthMolecules, double distXYSqrSum,
			double distZSqrSum, double totalDx, double totalDy, double totalDz, double photonsScale, double photonsDelta);
		void Print();
	};

	//GroundTruthStats
	class GroundTruthStats {
	public:

		GroundTruthStats() {};
		GroundTruthStats(int _nTotalGts, int *_gtNrs, double *_gtXs, double *_gtYs, double *_gtZs, 
			double *_gtPhotons, int _nrGroundTruthMolecules, double _tolXY = 250e-9, double _tolZ = 500e-9);
		void AddLocations(int m, int *nrs, double *xs, double *ys, double *zs, double *photons);
		void RemoveLocations(int m, int *nrs, double *xs, double *ys, double *zs, double *photons);
		MatchingResult MatchFrame(int frameNr, Locations &locs);
		MatchingResult MatchFrame(int frameNr, int m, int *nrs, double *xs, double *ys, double *zs, double *photons);
		void Clear();
		Stats GetStat();

	private:
		double tolXY, tolZ;

		int nTotalGts;
		int *gtNrs;
		double *gtXs, *gtYs, *gtZs, *gtPhotons;
		std::map<int, Locations> gtLocsByFrames;
		int nrGroundTruthMolecules;

		std::map<int, Locations> addedLocsByFrames;
		std::map<int, MatchingResult> framesMatchingResults;
		int nTruePos, nFalsePos, nFalseNeg;
		double distXYSqrSum, distZSqrSum, totalDx, totalDy, totalDz, photonsScale, photonsDelta;

		void AddToStatistics(MatchingResult &r, Locations &addedLocs, double coef = 1.0);
		void RemoveFromStatistics(MatchingResult &r, Locations &addedLocs);
	};

}

#endif
