import overlap_graph
from overlap_graph import build_reads_to_overlap_edges_map
from overlap_graph import make_kmer_dict

from read_genome import read_fastq


def merge_most_overlapping_reads(overlap_length_dict, max_overlap_length):
    reads_to_merge = overlap_length_dict[max_overlap_length].pop()
    return reads_to_merge[0], reads_to_merge[1], r + s[max_overlap_length:]


def update_reads(reads, r, s, merged_read):
    reads.remove(r)
    reads.remove(s)
    reads.add(merged_read)

def update_edge_maps_and_max_overlap(reads_to_edges_map, overlap_length_dict,
                                      max_overlap_length,r,s,merged_read, k):

    

    return max_overlap_length  # if no more edges longer than k, return max_overlap_length of '-1'


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

    reads_to_edges_map, overlap_length_dict, max_overlap_length = \
        build_reads_to_overlap_edges_map(kmers_to_reads, reads, k)

    while max_overlap_length >= k:

        r, s, merged_read = merge_most_overlapping_reads(overlap_length_dict, max_overlap_length)

        update_reads(reads, r, s, merged_read)

        max_overlap_length = \
            update_edge_maps_and_max_overlap(reads_to_edges_map, overlap_length_dict,
                                              max_overlap_length,r,s,merged_read, k)



    # if assembly is not complete with an overlap lower bound of 30, repeat with smaller lower bound
    # (possibly start with 40)



reads, qualities = read_fastq('ads1_week4_reads.fq')
greedy_scs(reads)



