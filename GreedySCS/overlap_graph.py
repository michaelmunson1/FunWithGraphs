from read_genome import read_fastq

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match


# a = 'abcd'
# b = a[1:3]
# print(type(b))


# make a dict of sets keyed on all 30-mers, whose values are the reads which contain that 30-mer
def make_kmer_dict(reads, k):
    kmers_to_reads = {}

    n = 0

    for r in reads:
        for i in xrange(len(r) - k + 1):
            kmer = r[i:i+k]

            # if n < 2 and i > len(r) - k - 5:
            #     print n
            #     print(kmer)

            if kmer not in kmers_to_reads:
                kmers_to_reads[kmer] = set([r])
            else:
                kmers_to_reads[kmer].add(r)
        n += 1

    print(n)

    return kmers_to_reads



# for each read 'a', we find all overlaps containing a suffix of b
def build_reads_to_overlap_edges_map(kmers_to_reads, reads, k):
    """
    Returns reads_to_edges_map, overlap_length_dict, max_overlap_length:
    reads_to_edges_map: a graph of overlapping reads, stored as a map from reads to two lists, of incoming and outgoing edges,
    each stored with the length of the overlap.

    overlap_length_dict: maps overlap lengths to read pairs (edges)

    max_overlap_length: for constant-time identification of next reads to merge

    Only creates edges for overlaps of length at least 30.

    :param kmers_to_reads: dict storing the read strings in which each kmer in 'reads' appears
    :param reads: strings of length 100 bases (A,C,T,G)
    :param k: shortest overlaps to include in graph
    """

    overlap_length_dict = {}    # map from length of overlap list of nodes with that
    max_length = 0

    for i in xrange(1,101):
        overlap_length_dict[i] = []

    reads_to_edges_map = {}     # entry format, for read 'r': ([outgoing edges list],  [incoming edges list]) ,
    # where an entry of the first list is in the form (dest_node, overlap_length),
    # and an entry of the second list is in the form (src_node, overlap_length),
    # where nodes are represented by read strings

    #initialize map
    for r in reads:
        reads_to_edges_map[r] = ([],[])

    #num_edges = 0
    num_nodes_with_outgoing_edge = 0
    for r in reads:
        # found_overlap = False
        r_suffix = r[(-1 * k):]      # ASSUMPTION: (for now) we don't have any reads less than 30 chars long
        #print(len(r_suffix))
        for s in kmers_to_reads[r_suffix]:     #  this suffix should already be present in kmers_to_reads
            if r != s:          # check for read equality, *** DON'T include these ***
                overlap_length = overlap(r,s,k)    # suffix of r (length >= k) occurs as prefix of s

                if overlap_length >= k:
                    #graph.append(((r, s), overlap_length))

                    reads_to_edges_map[r][0].append((s, overlap_length))
                    reads_to_edges_map[s][1].append((r, overlap_length))

                    # num_edges += 1

                    # if found_overlap is False:
                    #     found_overlap = True
                    #     num_nodes_with_outgoing_edge += 1

                    overlap_length_dict[overlap_length].append((r, s))  # map length to edges and store longest overlap

                    if overlap_length > max_overlap_length:
                        max_overlap_length = overlap_length

    #print num_edges
    #print num_nodes_with_outgoing_edge
    #print graph

    return reads_to_edges_map, overlap_length_dict, max_overlap_length


#reads, qualities = read_fastq('ERR266411_1.for_asm.fastq')      # N.B., no two reads are the same in this file
#k = 30
#reads = ['ABCDEFG', 'EFGHIJ', 'HIJABC']

#reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
#k = 5

#kmers_to_reads = make_kmer_dict(reads, k)
#graph= construct_overlap_graph(kmers_to_reads, reads, k)