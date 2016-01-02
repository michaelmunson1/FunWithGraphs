# import collections
# import matplotlib.pyplot as plt

# from naive_2mm import naive_2mm
# from naive_with_rc import naive_with_rc

def q_to_phred33(Q):
    """

    :param Q: base quality, rounded to nearest integer
    :return:
    """
    return chr(Q+33)

def phred33_to_q(quality_char):
    """
    :param quality_char: base quality char (ASCII)
    :return: integer base quality
    """
    return ord(quality_char) - 33



def read_genome(filename):
    """
    :rtype: string
    """
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def read_fastq(filename):
    """
    :param: FASTQ file name
    :return: list of read strings, and associated list of quality strings
    """
    sequences = []
    qualities = []

    with open(filename, 'r') as fh:    #file handle
        while True:
            fh.readline()   # tag for read
            seq = fh.readline().rstrip()    # seq and qual lines have \n at end
            fh.readline()
            qual = fh.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


def create_hist(qualities):
    hist = [0] * 50
    for qual in qualities:
        for phred in qual:
            q = phred33_to_q(phred)
            hist[q] += 1
    return hist


def FindGCByPos(reads):

    gc = [0] * 100
    totals = [0] * 100

    for read in reads:

        for i in range(len(read)):
            if read[i] == 'G' or read[i] == 'C':
                gc[i] += 1
            totals[i] += 1

    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])

    return gc


def naive(p,t):
    occurrences = []

    for i in range(len(t) - len(p) + 1):
        is_match = True
        for j in range(len(p)):
            if p[j] != t[i+j]:
                is_match = False
                break
        if is_match:
            occurrences.append(i)

    return occurrences


def reverse_complement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t








#if __name__ == '__main__':




#
# counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
#
# for base in genome:
#     if base not in counts:
#         print("Error! Expected 'A', 'C', 'T', or 'G' ")
#     counts[base] += 1
#
# print(counts)
#
# print(collections.Counter(genome))


def find_low_quality_cycle(quals):

    read_length = len(quals[0])

    qual_scores = [0] * read_length
    num_reads = 0

    for read in quals:
        num_reads += 1
        for i in range(read_length):
            qual_scores[i] += phred33_to_q(read[i])

    min_index = -1
    min_qual_score = 99

    for i in range(read_length):
        qual_scores[i] /= float(num_reads)
        if qual_scores[i] < min_qual_score:
            min_index = i
            min_qual_score = qual_scores[i]

    print(qual_scores)

    print(min(qual_scores))

    print(qual_scores[min_index])


    return min_index


#_, quals = read_fastq('ERR037900_1.first1000.fastq')
#min_index = find_low_quality_cycle(quals)

#print min_index






# print(seqs[:5])
# print(quals[:5])
#
# h = create_hist(quals)
#
# print(h)
#
# plt.bar(range(len(h)), h)
# plt.show()

# gc = FindGCByPos(seqs)
#
# print(gc)
#
# plt.plot(range(len(gc)), gc)
# plt.show()


# count = collections.Counter()
#
# for seq in seqs:
#     count.update(seq)
#
# print(count)

# genome = read_genome('phix.fa')
#
# print(genome)
# print(len(genome))
#
# print(naive('world', 'hello world'))





