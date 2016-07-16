from _statsCpp import PyGroundTruthStats # @UnresolvedImport

class GroundTruthStatsCpp(object):
    
    def __init__(self, frameNrs, coords, photons, nrGroundTruthMolecules, tolXY = 250e-9, tolZ = 500e-9):
        self._stats = PyGroundTruthStats(frameNrs, coords[0], coords[1], \
                                         coords[2], photons, nrGroundTruthMolecules, \
                                         tolXY, tolZ)

    def AddLocations(self, frameNrs, coords, photons):
        # return (distXY, distZ)
        return self._stats.AddLocations(frameNrs, coords[0], coords[1], coords[2], photons)
    
    def RemoveLocations(self, frameNrs, coords, photons):
        self._stats.RemoveLocations(frameNrs, coords[0], coords[1], coords[2], photons)
    
    def MatchFrame(self, frameToMatch, frameNrs, coords, photons):
        return self._stats.MatchFrame(frameToMatch, frameNrs, coords[0], coords[1], coords[2], photons)
        
    def Clear(self):
        return self._stats.Clear()
        
    def GetStats(self):
        return self._stats.GetStat()

if __name__ == "__main__":
    pass