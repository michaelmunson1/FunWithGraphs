import overlap_graph
from overlap_graph import build_reads_to_overlap_edges_map
from overlap_graph import make_kmer_dict

from read_genome import read_fastq


def merge_most_overlapping_reads(overlap_length_dict, max_overlap_length):

    reads_to_merge_dict = overlap_length_dict[max_overlap_length]
    r = reads_to_merge_dict.keys()[0]

    #print 'length of dict is %d' % len(reads_to_merge_dict)

    print 'max_olen is %d' % max_overlap_length
    print 'r is %s' % r

    s = reads_to_merge_dict[r]

    print 's is %s' % s

    return r, s, r + s[max_overlap_length:]


def update_reads(reads, r, s, merged_read):
    reads.remove(r)
    reads.remove(s)
    reads.append(merged_read)


def update_edge_maps_and_max_overlap(reads_to_edges_map, overlap_length_dict,
                                      max_overlap_length,r,s,merged_read, k):

    # remove (r,s) edge
    del reads_to_edges_map[r][0][s]
    del reads_to_edges_map[s][1][r]
    del overlap_length_dict[max_overlap_length][r]

    # update (r, q) and (q',r) edges for q,q' != s
    reads_to_edges_map[merged_read] = [{},{}]

    # could decompose these functions out:
    # update_r_edges(reads_to_edges_map, r, s, merged_read)
    # update_s_edges(reads_to_edges_map, s, r, merged_read)


    for q in reads_to_edges_map[r][0]:              # r's outgoing edges
                                                    # ordering in genome: (r,s,q) - since olap(r,s) > olap(r,q)

        overlap_length = reads_to_edges_map[r][0][q]
        del reads_to_edges_map[q][1][r]
        del overlap_length_dict[overlap_length][r]

        new_overlap_length = reads_to_edges_map[s][0][q]

        reads_to_edges_map[merged_read][0][q] = new_overlap_length
        reads_to_edges_map[q][1][merged_read] = new_overlap_length
        overlap_length_dict[new_overlap_length][merged_read] = q


    for q in reads_to_edges_map[r][1]:              # r's incoming edges:
                                                    # ordering in genome: (q,r,s) - since olap(r,s) > olap(r,q)
                                                    # thus, q may or may not overlap s
        new_overlap_length = reads_to_edges_map[q][0][r]

        del overlap_length_dict[new_overlap_length][q]

        del reads_to_edges_map[q][0][r]
        reads_to_edges_map[q][0][merged_read] = new_overlap_length
        reads_to_edges_map[merged_read][1][q] = new_overlap_length

        if s in reads_to_edges_map[q][0]:
            overlap_length = reads_to_edges_map[q][0][s]
            del overlap_length_dict[overlap_length][q]
            del reads_to_edges_map[q][0][s]

    for q in reads_to_edges_map[s][0]:              # s' outgoing edges
                                                    # genome ordering: (r,s,q), where r and q may or may not overlap
        new_overlap_length = reads_to_edges_map[s][0][q]
        #print 'new_overlap_length is %d' % new_overlap_length
        del overlap_length_dict[new_overlap_length][s]     # there is a unique entry,
                                                           # since redundant reads have been discarded


        reads_to_edges_map[merged_read][1][q] = new_overlap_length
        del reads_to_edges_map[q][1][s]

    del reads_to_edges_map[r]
    del reads_to_edges_map[s]

    if len(overlap_length_dict[max_overlap_length]) > 0:
        pass
    else:
        length = max_overlap_length - 1
        while len(overlap_length_dict[max_overlap_length]) == 0 and length >= k:
            length -= 1

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
        build_reads_to_overlap_edges_map(kmers_to_reads, reads, k)   # stores overlaps of length k to 99, inclusive

    while max_overlap_length >= k:

        r, s, merged_read = merge_most_overlapping_reads(overlap_length_dict, max_overlap_length)

        update_reads(reads, r, s, merged_read)

        max_overlap_length = \
            update_edge_maps_and_max_overlap(reads_to_edges_map, overlap_length_dict,
                                             max_overlap_length,r,s,merged_read, k)

    # TODO: test for completion - how?

    # if assembly is not complete with an overlap lower bound of 30, repeat with smaller lower bound
    # (possibly start with 40)



reads, qualities = read_fastq('ads1_week4_reads.fq')
#reads = list(set(reads))     # remove duplicates - 14 duplicates
greedy_scs(reads)