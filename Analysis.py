
from FileWritingUtilFuncs import write_filterered_sorted_by_hit_multiplication, count_db_hits_and_write_filtered, \
    write_nondeamed_version, write_unfiltered
from Peptide_models import PosFreqCounter, FullSeqModel, MixedModel, setup_models
from util_functions import sampleFromFreqDict, getSeqListFromFasta, getOutputFile, get_sampled_seq_counts, \
    modify_to_non_deaminated_version, read_peptides, check_alpha_gamma


def main():
    peptides = read_peptides()
    mm = setup_models(peptides)

    sortedSampledSeqCounts = get_sampled_seq_counts(mm)
    write_unfiltered(sortedSampledSeqCounts)

    modify_to_non_deaminated_version(sortedSampledSeqCounts)
    write_nondeamed_version(sortedSampledSeqCounts)

    spm_times_hits = count_db_hits_and_write_filtered(sortedSampledSeqCounts)

    write_filterered_sorted_by_hit_multiplication(spm_times_hits)


main()

#extended = extend_seqs([x.strip() for x in open("seqs_to_extend.txt").readlines()], getSeqListFromFasta("Complete_database.fasta"),2)
#write_extensions(extended, "extended_seqs.txt")
#check_alpha_gamma(queries, "Complete_database.fasta")
