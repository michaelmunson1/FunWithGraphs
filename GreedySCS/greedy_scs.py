from overlap_graph import build_reads_to_overlap_edges_map
from overlap_graph import make_kmer_dict

from unittest import TestCase

from read_genome import read_fastq

from time import clock


def merge_most_overlapping_reads(overlap_length_dict, max_overlap_length, num_read_uids, uids_to_reads):

    reads_to_merge_dict = overlap_length_dict[max_overlap_length] # find longest overlapping reads

    # TODO: simplify - ignoring identical reads -> lists not required
    last_entry = len(reads_to_merge_dict) - 1
    r_uid = reads_to_merge_dict.keys()[last_entry]
    last_entry = len(reads_to_merge_dict[r_uid]) - 1
    s_uid = reads_to_merge_dict[r_uid][last_entry]

    remove_or_delete_from_overlap_length_dictionary(overlap_length_dict, max_overlap_length, r_uid)

    r = uids_to_reads[r_uid]
    s = uids_to_reads[s_uid]

    merged_read = r + s[max_overlap_length:]
    uids_to_reads[num_read_uids] = merged_read

    return r_uid, s_uid


def remove_merged_reads(uids_to_reads, r_uid, s_uid):

    del uids_to_reads[r_uid]
    del uids_to_reads[s_uid]


def remove_or_delete_from_overlap_length_dictionary(overlap_length_dict, this_length, uid):

    if len(overlap_length_dict[this_length][uid]) == 1:
        del overlap_length_dict[this_length][uid]
    else:
        overlap_length_dict[this_length][uid].pop()

def update_edge_maps_and_max_overlap(reads_to_edges_map, overlap_length_dict,
                                     max_overlap_length, r_uid, s_uid, merged_read_uid, k, uids_to_reads):
    """
    Updates main overlap graph, and max_overlap_length; called upon merging the two
    most overlapping reads (suffix of 'r' = prefix of 's')

    :param reads_to_edges_map:
        Dictionary keyed on read uids. Value: [outgoing edges dict,  incoming edges list],
        where the first element maps destination_nodes to overlap_length,
        the second maps source_nodes to overlap_length,
        and nodes are represented by read uids
    :param overlap_length_dict: uses double-key [r][overlap_length] , value: s
    :return: updated max_overlap_length
    """

    # remove (r,s) edge
    del reads_to_edges_map[r_uid][0][s_uid]
    del reads_to_edges_map[s_uid][1][r_uid]

    # update (r, q) and (q',r) edges for q,q' != s
    reads_to_edges_map[merged_read_uid] = [{}, {}]

    # TODO: refactor, further decompose these functions

    for q_uid in reads_to_edges_map[r_uid][0]:      # r's outgoing edges:
                                                    # ordering in genome: (r,s,q) - since olap(r,s) > olap(r,q)
                                                    # but: must deal with possibility that s & q are the same!

        overlap_length = reads_to_edges_map[r_uid][0][q_uid]

        del reads_to_edges_map[q_uid][1][r_uid]

        remove_or_delete_from_overlap_length_dictionary(overlap_length_dict, overlap_length, r_uid)

        new_overlap_length = reads_to_edges_map[s_uid][0][q_uid]

        reads_to_edges_map[merged_read_uid][0][q_uid] = new_overlap_length
        reads_to_edges_map[q_uid][1][merged_read_uid] = new_overlap_length

        overlap_length_dict[new_overlap_length][merged_read_uid] = [q_uid]

        if new_overlap_length > max_overlap_length:
            max_overlap_length = new_overlap_length



    for q_uid in reads_to_edges_map[r_uid][1]:      # r's incoming edges:
                                                    # ordering in genome: (q,r,s) - since olap(r,s) > olap(r,q)
                                                    # thus, q may or may not overlap s

        new_overlap_length = reads_to_edges_map[q_uid][0][r_uid]

        overlap_length_dict[new_overlap_length][q_uid] = [merged_read_uid]

        remove_or_delete_from_overlap_length_dictionary(overlap_length_dict, new_overlap_length, q_uid)


        overlap_length_dict[new_overlap_length][q_uid] = [merged_read_uid]

        del reads_to_edges_map[q_uid][0][r_uid]
        reads_to_edges_map[q_uid][0][merged_read_uid] = new_overlap_length
        reads_to_edges_map[merged_read_uid][1][q_uid] = new_overlap_length

        if s_uid in reads_to_edges_map[q_uid][0]:
            overlap_length = reads_to_edges_map[q_uid][0][s_uid]
            del overlap_length_dict[overlap_length][q_uid]
            del reads_to_edges_map[q_uid][0][s_uid]

    for q_uid in reads_to_edges_map[s_uid][0]:              # s' outgoing edges
                                                    # genome ordering: (r,s,q), where r and q may or may not overlap
        new_overlap_length = reads_to_edges_map[s_uid][0][q_uid]

        remove_or_delete_from_overlap_length_dictionary(overlap_length_dict, new_overlap_length, s_uid)

        reads_to_edges_map[merged_read_uid][0][q_uid] = new_overlap_length
        reads_to_edges_map[q_uid][1][merged_read_uid] = new_overlap_length
        del reads_to_edges_map[q_uid][1][s_uid]

        overlap_length_dict[new_overlap_length][merged_read_uid] = [q_uid]

    del reads_to_edges_map[r_uid]
    del reads_to_edges_map[s_uid]

    if len(overlap_length_dict[max_overlap_length]) > 0:
        pass
    else:
        length = max_overlap_length - 1

        while len(overlap_length_dict[length]) == 0 and length >= k:
            length -= 1

        max_overlap_length = length

        merged = uids_to_reads[merged_read_uid]


    return max_overlap_length  # if no more edges longer than k, return max_overlap_length of k-1


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
    k = 50
    kmers_to_reads, uids_to_reads, num_read_uids = make_kmer_dict(reads, k)

    # reads_to_edges_map stores overlaps of length k to 99, inclusive
    reads_to_edges_map, overlap_length_dict, max_overlap_length = \
        build_reads_to_overlap_edges_map(kmers_to_reads, uids_to_reads, num_read_uids, k)

    # TODO: implement iterative lowering of an overlap length lower bound (say 50, 40, 30)

    while max_overlap_length >= k:

        r_uid, s_uid = merge_most_overlapping_reads(overlap_length_dict, max_overlap_length,
                                                    num_read_uids, uids_to_reads)

        num_read_uids += 1

        remove_merged_reads(uids_to_reads, r_uid, s_uid)

        merged_read_uid = num_read_uids - 1

        max_overlap_length = \
            update_edge_maps_and_max_overlap(reads_to_edges_map, overlap_length_dict,
                                             max_overlap_length, r_uid, s_uid, merged_read_uid, k, uids_to_reads)


reads, qualities = read_fastq('ads1_week4_reads.fq')
reads = list(set(reads))     # remove duplicates - 14 duplicates

test = TestCase()
start = clock()
test.assertEquals(greedy_scs(reads), 15894)
elapsed = clock() - start
print 'time: %f' % elapsed

