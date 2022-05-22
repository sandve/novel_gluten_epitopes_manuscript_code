from random import random
from re import finditer
from collections import defaultdict

from Config import PEPTIDES_FN, PSEUDO_COUNTS, DEPTH
from Peptide import Peptide


def getSeqListFromFasta(fn):
    seqList = []
    seqInProgress = ""

    for line in open(fn):
        if line.startswith(">"):
            if seqInProgress != "":
                seqList.append(seqInProgress)
            seqInProgress = ""
        else:
            seqInProgress += line.strip()
    if seqInProgress != "":
        seqList.append(seqInProgress)
    return seqList

def sampleFromFreqDict(freqDict):
    aas = list(freqDict.keys())
    thresholds = []
    threshold = 0
    for aa in aas:
        threshold += freqDict[aa]
        thresholds.append(threshold)
    assert thresholds[-1] == sum(freqDict.values())
    draw = random() * thresholds[-1]
    for index, t in enumerate(thresholds):
        if draw < t:
            return aas[index]
    raise


def getOutputFile(fn):
    outFile = open(fn,'w')
    outFile.write("Pseudo-count: " + f'{PSEUDO_COUNTS:.2f}' +'\n')
    outFile.write("Peptides: " + PEPTIDES_FN + '\n')
    outFile.write("Simulation depth: " + str(DEPTH) + '\n')
    outFile.write("spm means seq per million (estimated probability times a million)\n")
    return outFile


def get_sampled_seq_counts(mm):
    sampledSeqCounts = defaultdict(int)
    for i in range(DEPTH):
        sampledSeqCounts[mm.sample()] += 1
    sortedSampledSeqCounts = sorted(sampledSeqCounts.items(), key=lambda x: x[1], reverse=True)
    return sortedSampledSeqCounts


def modify_to_non_deaminated_version(sortedSampledSeqCounts):
    for i, (s, c) in enumerate(sortedSampledSeqCounts):
        if s[3] == "E":
            sl = list(sortedSampledSeqCounts[i][0])
            sl[3] = "Q"
            sortedSampledSeqCounts[i] = list(sortedSampledSeqCounts[i])
            sortedSampledSeqCounts[i][0] = "".join(sl)
        if s[5] == "E":
            sl = list(sortedSampledSeqCounts[i][0])
            sl[5] = "Q"
            sortedSampledSeqCounts[i] = list(sortedSampledSeqCounts[i])
            sortedSampledSeqCounts[i][0] = "".join(sl)


def read_peptides():
    peptides = []
    deamPositions = defaultdict(int)
    for line in open(PEPTIDES_FN):
        name, pep = line.strip().split()
        peptides.append(Peptide(pep))
        deamPositions[tuple([x.start() for x in finditer("[.]", peptides[-1].base_seq)])] += 1
    return peptides


def extend_seqs(seq_list, db_seqs, flank_size):
    extensions = {}
    for seq in seq_list:
        extensions[seq] = []
        for db_seq in db_seqs:
            for hit in finditer(seq, db_seq):
                start = max(0, hit.start()-flank_size)
                end = min(len(db_seq), hit.end()+flank_size)
                extensions[seq].append(db_seq[start:end])
    return extensions


def check_alpha_gamma(queries, fasta_db_fn):
    entries = {}
    for line in open(fasta_db_fn):
        if line.startswith(">"):
            entries[line.strip()] = ""
            lastKey = line.strip()
        else:
            entries[lastKey] += line.strip()

    for query in queries:
        chains = set()
        for id,seq  in entries.items():
            if query in seq:
                chains.add(id.split(" ")[1])
