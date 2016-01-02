import construct_overlap_graph
from construct_overlap_graph import construct_overlap_graph
from construct_overlap_graph import make_kmer_dict

from read_genome import read_fastq


def greedy_scs(reads):
    """
    :param reads: sequencing read strings.  In the context of this problem, we assume:
        - reads of length 100
        - no sequencing errors
        - no polyploidy; all reads from the forward strand
    :return: super: a greedily approximated shortest common superstring representing the reconstructed genome
    """

    # Given the number of reads in the genome, we construct an overlap graph
    # connecting only those reads that overlap at least 30 bases
    k = 30
    kmers_to_reads = make_kmer_dict(reads, k)
    overlap_edges, overlap_length_dict, max_length = construct_overlap_graph(kmers_to_reads, reads, k)

    while len(overlap_edges) > 0:       # while there still remain edges of length > 30
        most_overlapping_edge = overlap_length_dict[max_length][0]



    # if assembly is not complete with an overlap lower bound of 30, repeat with smaller lower bound
    # (possibly start with 40)



reads, qualities = read_fastq('ads1_week4_reads.fq')
greedy_scs(reads)



