from collections import defaultdict

from Config import ALL_AAS_EXCEPT_E, ALL_AAS_EXCEPT_PE, PSEUDO_COUNTS
from util_functions import sampleFromFreqDict


class PosFreqCounter:
    def __init__(self, attribute, peptides):
        self._attribute = attribute
        self._peptides = peptides
        self.freqs = defaultdict(lambda:defaultdict(int))
        self.count_freqs()


    def count_freqs(self):

        for pep in [p.__dict__[self._attribute] for p in self._peptides]:
            if pep is not None:
                for pos,aa in enumerate(pep):
                    if aa != '.':
                        self.freqs[pos][aa] += 1.0

    def add_pseudo_counts(self, positions, aas, weight):
        for pos in positions:
            for aa in aas:
                self.freqs[pos][aa] += weight

    def __repr__(self):
        return str(self.__dict__)


class PosSampler:
    def __init__(self, freqDict, index):
        assert len(freqDict)>0, index
        self._freqDict = freqDict

    def sample(self):
        return sampleFromFreqDict(self._freqDict)

    def getProb(self, symbol):
        return self._freqDict[symbol] / sum(self._freqDict.values())


class FullSeqModel:
    def __init__(self, freqDictList):
        self._posDistrList = [PosSampler(freq, index) for index, freq in enumerate(freqDictList)]

    def sample(self):
        seq = ""
        for pos,distr in enumerate(self._posDistrList):
            seq += distr.sample()
        return seq

    def getSeqProbs(self, seq):
        print(seq)
        full_prob = 1
        for i,s in enumerate(seq):
            prob = self._posDistrList[i].getProb(s)
            print("{:.3f}".format(prob),end=' ')
            full_prob *= prob
        print('\n',"{:.4f}".format(full_prob*100),'%')


class MixedModel:
    def __init__(self, weightedFullSeqModels):
        self._weightedFullSeqModels = weightedFullSeqModels

    def sample(self):
        modelItems = list(self._weightedFullSeqModels.items())
        freqDict = dict([(index,item[1]) for index,item in enumerate(modelItems)])
        index = sampleFromFreqDict(freqDict)
        return modelItems[index][0].sample()


def setup_models(peptides):
    bf = PosFreqCounter("base_seq", peptides)
    sf = PosFreqCounter("shortDeam", peptides)
    lf = PosFreqCounter("longDeam", peptides)
    bf.add_pseudo_counts(list(range(5)) + list(range(6, 9)), ALL_AAS_EXCEPT_E, PSEUDO_COUNTS)
    sf.add_pseudo_counts([1], ALL_AAS_EXCEPT_PE, PSEUDO_COUNTS)
    sf.add_pseudo_counts([1, 2, 3], ALL_AAS_EXCEPT_PE, PSEUDO_COUNTS)
    fsm1 = FullSeqModel(
        [bf.freqs[i] for i in range(5)] + [sf.freqs[i] for i in range(3)] + [bf.freqs[i] for i in range(8, 9)])
    fsm2 = FullSeqModel(
        [bf.freqs[i] for i in range(3)] + [sf.freqs[i] for i in range(3)] + [bf.freqs[i] for i in range(6, 9)])
    fsm3 = FullSeqModel(
        [bf.freqs[i] for i in range(3)] + [lf.freqs[i] for i in range(2)] + [sf.freqs[i] for i in range(3)] + [
            bf.freqs[i] for i in range(8, 9)])
    mm = MixedModel({fsm1: 6 / 16, fsm2: 8 / 16, fsm3: 2 / 16})
    return mm