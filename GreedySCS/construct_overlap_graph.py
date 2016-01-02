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
def construct_overlap_graph(kmers_to_reads, reads, k):

    overlap_length_dict = {}    # map from length of overlap list of nodes with that
    max_length = 0

    graph = []                  # entry (edge) format: ( (src_node, dest_node), overlap_length ),
                                # where nodes are represented by read strings
    num_edges = 0
    num_nodes_with_outgoing_edge = 0
    for r in reads:
        found_overlap = False
        r_suffix = r[(-1 * k):]      # ASSUMPTION: (for now) we don't have any reads less than 30 chars long
        #print(len(r_suffix))
        for s in kmers_to_reads[r_suffix]:     #  this suffix should already be present in kmers_to_reads
            if r != s:          # check for read equality, *** DON'T include these ***
                # overlap_length = overlap(r,s,k)    # suffix of r (length >= k) occurs as prefix of s

                # if overlap_length not in overlap_length_dict:
                #     overlap_length_dict[overlap_length] = [(r, s)]
                # else:
                #     overlap_length_dict[overlap_length].append((r, s))
                #
                # if overlap_length > max_overlap_length:
                #     max_overlap_length = overlap_length

                if overlap_length >= k:
                    if found_overlap is False:
                        found_overlap = True
                        num_nodes_with_outgoing_edge += 1
                    graph.append(((r, s), overlap_length))       # has_been_traversed = False
                    num_edges += 1

    print num_edges
    print num_nodes_with_outgoing_edge
    #print graph

    return graph, overlap_length_dict


#reads, qualities = read_fastq('ERR266411_1.for_asm.fastq')      # N.B., no two reads are the same in this file
#k = 30
#reads = ['ABCDEFG', 'EFGHIJ', 'HIJABC']

#reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
#k = 5

#kmers_to_reads = make_kmer_dict(reads, k)
#graph= construct_overlap_graph(kmers_to_reads, reads, k)