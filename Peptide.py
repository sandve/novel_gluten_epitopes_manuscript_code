from re import search, findall


class Peptide:
    def __init__(self, seq):
        self.origseq = seq
        self.base_seq = seq
        self.shortDeam = self._extractPattern("E[^P]P")
        self.longDeam = self._extractPattern("E[^P][^P][FYWMLIV]")

    def __repr__(self):
        return str(self.__dict__)

    def _extractPattern(self, pattern):
        hit = search(pattern, self.origseq)
        if hit is not None:
            assert len(findall(pattern, self.origseq)) == 1
            deam = self.origseq[hit.start(): hit.end()]
            self.base_seq = self.base_seq[0:hit.start()] + '.' * (hit.end() - hit.start()) + self.base_seq[hit.end():]
            return deam
        return None