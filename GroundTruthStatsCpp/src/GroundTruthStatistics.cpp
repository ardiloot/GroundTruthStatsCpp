#include "GroundTruthStatistics.h"

namespace GroundTruthStatistics {

	// Methods

	double sqr(double x) {
		return x * x;
	}

	// Locations

	Locations::Locations() {
		// reserve space
		/*
		int reserveSize = 20;
		nrs.reserve(reserveSize);
		xs.reserve(reserveSize);
		ys.reserve(reserveSize);
		zs.reserve(reserveSize);
		photons.reserve(reserveSize);
		indices.reserve(reserveSize);
		*/
	}

	void Locations::Add(int nr, double x, double y, double z, double photon, int index) {
		nrs.push_back(nr);
		xs.push_back(x);
		ys.push_back(y);
		zs.push_back(z);
		photons.push_back(photon);
		indices.push_back(index);
	}

	void Locations::Remove(int nr, double x, double y, double z, double photon) {
		int index = IndexOf(nr, x, y, z, photon);
		nrs.erase(nrs.begin() + index);
		xs.erase(xs.begin() + index);
		ys.erase(ys.begin() + index);
		zs.erase(zs.begin() + index);
		photons.erase(photons.begin() + index);
		indices.erase(indices.begin() + index);
	}

	int Locations::IndexOf(int nr, double x, double y, double z, double photon) {
		for (int i = 0; i < size(); i++) {
			if (nr == nrs[i] && x == xs[i] && y == ys[i] && z == zs[i] && photon == photons[i]) {
				return i;
			}
		}
		return -1;
	}

	int Locations::size() {
		return nrs.size();
	}

	// MatchingResult
	
	MatchingResult::MatchingResult(int _frameNr) {
		frameNr = _frameNr;

		// reserve space
		int reserveSize = 20;
		idsTest.reserve(reserveSize);
		idsTrue.reserve(reserveSize);
		distsXY.reserve(reserveSize);
		distsZ.reserve(reserveSize);
	}

	void MatchingResult::Add(int idTest, int idTrue, double distXY, double distZ) {
		idsTest.push_back(idTest);
		idsTrue.push_back(idTrue);
		distsXY.push_back(distXY);
		distsZ.push_back(distZ);
	}

	int MatchingResult::size() {
		return idsTest.size();
	}

	// Stats

	Stats::Stats(int nTruePos, int nFalsePos, int nFalseNeg, int nrGroundTruthMolecules, double distXYSqrSum,
		double distZSqrSum, double totalDx, double totalDy, double totalDz, double photonsScale, double photonsDelta) {
		nTP = nTruePos;
		nFP = nFalsePos;
		nFN = nFalseNeg;
		int actMolsProcessed = nFN + nTP;
		int actMolsNotProcessed = nrGroundTruthMolecules - actMolsProcessed;
		nFN += actMolsNotProcessed;

		rmsXY = sqrt(distXYSqrSum / nTP);
		rmsZ = sqrt(distZSqrSum / nTP);
		jac = 100.0 * double(nTP) / (nFN + nFP + nTP);
		recall = double(nTP) / (nFN + nTP);
		precision = double(nTP) / (nFP + nTP);
		averageDx = totalDx / nTP;
		averageDy = totalDy / nTP;
		averageDz = totalDz / nTP;
		avgPhotonsScale = photonsScale / nTP;
		avgPhotonsDelta = photonsDelta / nTP;
	}

	void Stats::Print() {
		std::cout << "nTP: " << nTP << std::endl;
		std::cout << "nFP: " << nFP << std::endl;
		std::cout << "nFN: " << nFN << std::endl;
		std::cout << "rmsXY: " << 1e9 * rmsXY << std::endl;
		std::cout << "rmsZ: " << 1e9 * rmsZ << std::endl;
		std::cout << "jac: " << jac << std::endl;
		std::cout << "recall: " << recall << std::endl;
		std::cout << "precision: " << precision << std::endl;
		std::cout << "averageDx: " << averageDx << std::endl;
		std::cout << "averageDy: " << averageDy << std::endl;
		std::cout << "averageDz: " << averageDz << std::endl;
		std::cout << "avgPhotonsScale: " << avgPhotonsScale << std::endl;
		std::cout << "avgPhotonsDelta: " << avgPhotonsDelta << std::endl;
	}

	// GroundTruthStats

	GroundTruthStats::GroundTruthStats(int _nTotalGts, int *_gtNrs, double *_gtXs, double *_gtYs, double *_gtZs, 
		double *_gtPhotons, int _nrGroundTruthMolecules, double _tolXY, double _tolZ) {
		nTotalGts = _nTotalGts; // Number of ground-truths
		gtNrs = _gtNrs;
		gtXs = _gtXs;
		gtYs = _gtYs;
		gtZs = _gtZs;
		gtPhotons = _gtPhotons;
		nrGroundTruthMolecules = _nrGroundTruthMolecules;
		tolXY = _tolXY;
		tolZ = _tolZ;

		for (int i = 0; i < nTotalGts; i++) {
			gtLocsByFrames[gtNrs[i]].Add(gtNrs[i], gtXs[i], gtYs[i], gtZs[i], gtPhotons[i], i);
		}
		Clear();
	}

	AddLocationsResult GroundTruthStats::AddLocations(int m, int *nrs, double *xs, double *ys, double *zs, double *photons) {
		// Find frames to match
		std::set<int> framesToMatch;
		for (int i = 0; i < m; i++) {
			framesToMatch.insert(nrs[i]);
		}

		//Remove stats from already match frames
		for (std::set<int>::iterator it = framesToMatch.begin(); it != framesToMatch.end(); ++it) {
			int frameNr = *it;
			std::map<int, MatchingResult>::iterator findIt = framesMatchingResults.find(frameNr);
			if (findIt == framesMatchingResults.end()) {
				continue;
			}
			//Frame already added previously, need to remove, match frame again and then read the statistics
			RemoveFromStatistics(findIt->second, addedLocsByFrames[frameNr]);
		}

		// Update locs by frames
		std::vector<std::pair<int, int> > inputIndexToLocsByFramesIndex;
		for (int i = 0; i < m; i++) {
			Locations &tmp = addedLocsByFrames[nrs[i]];
			tmp.Add(nrs[i], xs[i], ys[i], zs[i], photons[i], i);
			inputIndexToLocsByFramesIndex.push_back(std::make_pair(nrs[i], tmp.size() - 1));
		}

		// Add stats for every frame in framesToMatch
		for (std::set<int>::iterator it = framesToMatch.begin(); it != framesToMatch.end(); ++it) {
			int frameNr = *it;
			Locations &locs = addedLocsByFrames[frameNr];
			MatchingResult r = MatchFrame(frameNr, locs);
			framesMatchingResults[frameNr] = r;
			AddToStatistics(r, locs);
			//break;
		}

		// AddLocationsResult
		AddLocationsResult res;
		for (int i = 0; i < inputIndexToLocsByFramesIndex.size(); i++){
			int frameNr = inputIndexToLocsByFramesIndex[i].first;
			int indexInLocsByFrames = inputIndexToLocsByFramesIndex[i].second;
			Locations &tmp = addedLocsByFrames[frameNr];
			res.distsXY.push_back(tmp.distsXY[indexInLocsByFrames]);
			res.distsZ.push_back(tmp.distsZ[indexInLocsByFrames]);
		}
		return res;
	}

	void GroundTruthStats::RemoveLocations(int m, int *nrs, double *xs, double *ys, double *zs, double *photons) {
		// Find frames to match
		std::set<int> framesToMatch;
		for (int i = 0; i < m; i++) {
			framesToMatch.insert(nrs[i]);
		}

		//Remove stats from already match frames
		for (std::set<int>::iterator it = framesToMatch.begin(); it != framesToMatch.end(); ++it) {
			int frameNr = *it;
			std::map<int, MatchingResult>::iterator findIt = framesMatchingResults.find(frameNr);
			if (findIt == framesMatchingResults.end()) {
				continue;
			}
			//Frame already added previously, need to remove, match frame again and then read the statistics
			RemoveFromStatistics(findIt->second, addedLocsByFrames[frameNr]);
		}

		// Update locs by frames
		for (int i = 0; i < m; i++) {
			addedLocsByFrames[nrs[i]].Remove(nrs[i], xs[i], ys[i], zs[i], photons[i]);
		}

		// Add stats for every frame in framesToMatch
		for (std::set<int>::iterator it = framesToMatch.begin(); it != framesToMatch.end(); ++it) {
			int frameNr = *it;
			Locations &locs = addedLocsByFrames[frameNr];
			MatchingResult r = MatchFrame(frameNr, locs);
			framesMatchingResults[frameNr] = r;
			AddToStatistics(r, locs);
		}
	}

	MatchingResult GroundTruthStats::MatchFrame(int frameNr, Locations &locs) {
		// Get ground-truth for this frame
		Locations &gtLocsThisFrame = gtLocsByFrames[frameNr];

		if (locs.size() == 0 || gtLocsThisFrame.size() == 0) {
			return MatchingResult(frameNr);
		}

		// Init matrices
		int maxdim = std::max(locs.size(), gtLocsThisFrame.size());
		std::vector<std::vector<double> > distXY(locs.size(), std::vector<double>(gtLocsThisFrame.size()));
		std::vector<std::vector<double> > distZ(locs.size(), std::vector<double>(gtLocsThisFrame.size()));
		Matrix<double> cost(maxdim, maxdim);

		// Fill cost with ones
		for (int i = 0; i < maxdim; i++) {
			for (int j = 0; j < maxdim; j++) {
				cost(i, j) = 1.0;
			}
		}

		// Build cost and distance matrices
		for (int i = 0; i < locs.size(); i++) {
			for (int j = 0; j < gtLocsThisFrame.size(); j++) {
				distXY[i][j] = sqrt(sqr(locs.xs[i] - gtLocsThisFrame.xs[j]) + sqr(locs.ys[i] - gtLocsThisFrame.ys[j]));
				distZ[i][j] = fabs(locs.zs[i] - gtLocsThisFrame.zs[j]);
				if (distXY[i][j] > tolXY || distZ[i][j] > tolZ) {
					cost(i,j) = 1.0;
				}
				else {
					cost(i,j) = sqrt(sqr(distXY[i][j]) + sqr(distZ[i][j]));
				}
			}
		}

		// Hungarian matching
		Munkres<double> munkres;
		Matrix<bool> munkresRes = munkres.solve(cost);

		// Result
		locs.distsXY.resize(locs.size(), tolXY);
		locs.distsZ.resize(locs.size(), tolZ);
		MatchingResult res(frameNr);
		for (int i = 0; i < locs.size(); i++) {
			locs.distsXY[i] = tolXY;
			locs.distsZ[i] = tolZ;
			for (int j = 0; j < gtLocsThisFrame.size(); j++) {
				if (!munkresRes(i, j) || cost(i, j) >= 1.0) {
					continue;
				}
				res.Add(i, j, distXY[i][j], distZ[i][j]);
				locs.distsXY[i] = distXY[i][j];
				locs.distsZ[i] = distZ[i][j];
			}
		}

		return res;
	}

	MatchingResult GroundTruthStats::MatchFrame(int frameNr, int m, int *nrs, double *xs, double *ys, double *zs, double *photons){
		Locations locs;
		for (int i = 0; i < m; i++){
			locs.Add(nrs[i], xs[i], ys[i], xs[i], photons[i]);
		}
		return MatchFrame(frameNr, locs);
	}

	void GroundTruthStats::Clear() {
		nTruePos = 0;
		nFalsePos = 0;
		nFalseNeg = 0;
		distXYSqrSum = 0.0;
		distZSqrSum = 0.0;
		totalDx = 0.0;
		totalDy = 0.0;
		totalDz = 0.0;
		photonsScale = 0.0;
		photonsDelta = 0.0;
		addedLocsByFrames.clear();
		framesMatchingResults.clear();
	}

	Stats GroundTruthStats::GetStat() {
		Stats res(nTruePos, nFalsePos, nFalseNeg, nrGroundTruthMolecules, distXYSqrSum,
			distZSqrSum, totalDx, totalDy, totalDz, photonsScale, photonsDelta);
		return res;
	}

	void GroundTruthStats::AddToStatistics(MatchingResult &r, Locations &addedLocs, double coef) {
		Locations &gtLocsThisFrame = gtLocsByFrames[r.frameNr];
		for (int i = 0; i < r.size(); i++) {
			int testIndex = r.idsTest[i], trueIndex = r.idsTrue[i];
			distXYSqrSum += coef * sqr(r.distsXY[i]);
			distZSqrSum += coef * sqr(r.distsZ[i]);
			totalDx += coef * (addedLocs.xs[testIndex] - gtLocsThisFrame.xs[trueIndex]);
			totalDy += coef * (addedLocs.ys[testIndex] - gtLocsThisFrame.ys[trueIndex]);
			totalDz += coef * (addedLocs.zs[testIndex] - gtLocsThisFrame.zs[trueIndex]);
			photonsScale += coef * (addedLocs.photons[testIndex] / gtLocsThisFrame.photons[trueIndex]);
			photonsDelta += coef * (addedLocs.photons[testIndex] - gtLocsThisFrame.photons[trueIndex]);
		}

		int nTP = r.size();
		int nFP = addedLocs.size() - nTP;
		int nFN = gtLocsThisFrame.size() - nTP;
		nTruePos += int(coef) * nTP;
		nFalsePos += int(coef) * nFP;
		nFalseNeg += int(coef) * nFN;
	}

	void GroundTruthStats::RemoveFromStatistics(MatchingResult &r, Locations &addedLocs) {
		AddToStatistics(r, addedLocs, -1.0);
	}

}
