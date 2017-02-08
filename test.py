
def consecutive_groups(iterable):
    '''
    Get all possible substrings
    '''
    s = tuple(iterable)
    for size in range(1, len(s)+1):
        for index in range(len(s)+1-size):
            yield iterable[index:index+size]

#print list(consecutive_groups('ATCGTATA'))


def make_kmer_table(seqs, k):
    ''' Given dictionary (e.g. output of parse_fasta) and integer k,
        return a dictionary that maps each k-mer to the set of names
        of reads containing the k-mer. '''
    table = {}  # maps k-mer to set of names of reads containing k-mer
    for name, seq in seqs.iteritems():
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:i+k]
            if kmer not in table:
                table[kmer] = set()
            table[kmer].add(name)
    return table

