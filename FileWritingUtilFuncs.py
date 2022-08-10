from Config import MC_SAMPLING_DEPTH, OUT_FOLDER, NUM_DB_PRESENT_PEPTIDES_TO_KEEP, NUM_NON_FILTERED_PEPTIDES_TO_WRITE
from util_functions import getOutputFile, getSeqListFromFasta

def write_filterered_sorted_by_hit_multiplication(spm_times_hits):
    outFfilteredManyHits = getOutputFile(OUT_FOLDER + "/filtered_samples_multiply_by_hits.csv")
    print("TEMP2: ", len(spm_times_hits))
    spm_items = spm_times_hits.items()
    sorted_spm_items = sorted(spm_items, key=lambda x: x[1], reverse=True)
    for pep, spm in sorted_spm_items:
        outFfilteredManyHits.write(pep + '\t' + "{:.1f}".format(spm) + "spm*dbHits\n")
    outFfilteredManyHits.close()


def count_db_hits_and_write_filtered(sortedSampledSeqCounts):
    outFfiltered = getOutputFile(OUT_FOLDER + "/filtered_samples.csv")
    db_seqs = getSeqListFromFasta("Complete_database.fasta")
    foundSeqs = 0
    spm_times_hits = {}
    for s in sortedSampledSeqCounts:
        pep = s[0]
        # pep = s[0].replace("E", "Q")
        if any([pep in db_seq for db_seq in db_seqs]):
            foundSeqs += 1
            numHits = sum([db_seq.count(pep) for db_seq in db_seqs])
            outFfiltered.write(pep + '\t' + "{:.1f}".format(1e6 * s[1] / MC_SAMPLING_DEPTH) + "spm\n")
            spm_times_hits[pep] = 1e6 * s[1] * numHits / MC_SAMPLING_DEPTH

            if foundSeqs > NUM_DB_PRESENT_PEPTIDES_TO_KEEP:
                print("Stopping at configured sufficient number of db hits: ", NUM_DB_PRESENT_PEPTIDES_TO_KEEP)
                break
    outFfiltered.close()
    print("TEMP1: ", foundSeqs, len(spm_times_hits))
    return spm_times_hits


def write_nondeamed_version(sortedSampledSeqCounts):
    outFnoDeams = getOutputFile(OUT_FOLDER + "/de-deamed_samples.csv")
    # foundSeqs = 0
    for s in sortedSampledSeqCounts:
        # manual rule to exclude five Qs in a row, as considered not viable
        if "QQQQQ" in s[0]:
            continue
        outFnoDeams.write(s[0] + '\t' + "{:.1f}".format(1e6 * s[1] / MC_SAMPLING_DEPTH) + "spm\n")
        # foundSeqs += 1
        # if foundSeqs > NUM_DB_PRESENT_PEPTIDES_TO_KEEP:
            # break
    outFnoDeams.close()


def write_unfiltered(sortedSampledSeqCounts):
    outFnonFiltered = getOutputFile(OUT_FOLDER + "/samples.csv")
    for s in sortedSampledSeqCounts[0:NUM_NON_FILTERED_PEPTIDES_TO_WRITE]:
        outFnonFiltered.write(s[0] + '\t' + "{:.1f}".format(1e6 * s[1] / MC_SAMPLING_DEPTH) + "spm\n")
    outFnonFiltered.close()


def write_extensions(extensions, fn):
    outF = open(fn,'w')
    for seq in extensions:
        outF.write('\t'.join([seq] + [x+f'({extensions[seq].count(x)})' for x in set(extensions[seq])]) + '\n')