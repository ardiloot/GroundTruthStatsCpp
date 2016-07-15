import numpy as np
from GroundTruthStatsCpp import GroundTruthStatsCpp

def LoadLocsFromFile(filename):
    nrs, xs, ys, zs, photons = np.loadtxt(filename, unpack = True)
    nrs = nrs.astype(int)
    return {"nrs": np.ascontiguousarray(nrs), 
            "xs": np.ascontiguousarray(xs),
            "ys": np.ascontiguousarray(ys),
            "zs": np.ascontiguousarray(zs), 
            "photons": np.ascontiguousarray(photons)}
        
if __name__ == "__main__":
    trueLocs = LoadLocsFromFile("ground_truth_locs.txt")
    testLocs = LoadLocsFromFile("postprocessed_locs.txt")
    
    stats = GroundTruthStatsCpp(trueLocs["nrs"],
                                         (trueLocs["xs"], trueLocs["ys"], trueLocs["zs"]),
                                         trueLocs["photons"],
                                         1194)
    
    stats.AddLocations(testLocs["nrs"][:],
                       (testLocs["xs"][:], testLocs["ys"][:], testLocs["zs"][:]),
                       testLocs["photons"][:])
    print stats.GetStat()
    stats.Clear()