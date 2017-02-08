'''
To run: python driver.py

Finds the percentage of contaminants in our sequence reads

'''

from suffix_tree import SuffixTree
from local_align import *
import numpy as np
import time
import urllib
from test import make_kmer_table

def parse_fasta(fh):
    ''' Parse FASTA into a dictionary '''
    fa = {}
    name = None
    # Part 1: compile list of lines for each sequence
    for ln in fh:
        if ln[0] == '>':  # new sequence
            name = ln[1:].split()[0]
            fa[name] = []
        else:
            # append nucleotides to current sequence
            fa[name].append(ln.rstrip())
    # Part 2: join lists into strings
    for name, nuc_list in fa.iteritems():
        fa[name] = ''.join(nuc_list)  # join into one long string
    return fa

def parse_contam_fasta(fh):
    fh = open(fh)
    fa = {}
    dict3 = {}
    name = None

    for ln in fh:
        if ln[0] == '>':
            name = ln[1:].split()[0]
            fa[name] = []
            dict3[name] = []
        else:
            dict3[name].append(len(ln.rstrip()))
            fa[name].append(ln.rstrip())

    for name, nuc_list in fa.iteritems():
        fa[name] = ''.join(nuc_list)

    return fa





####################### Commands Start Here #############################

print "Parsing FASTA file.."
dict1 = parse_contam_fasta('sample_reads1.fa')
total_read_len = len(dict1)

'''
    #DeBugging
print total_read_len
print "# of contams", total_read_len - 5000
print "actual percentage", ((total_read_len - 5000)/float(total_read_len))*100, "%"
quit()
#'''

print "Getting Overlap..."
#get our list of overlap of sequences
overlap_list = bestBuddy(dict1)


###### Alternate Method (kmer overlap) #################################
#'''

print "Getting k-mer weights..."

#get table of k-mers of length 7
tab = make_kmer_table(dict1, 7)
kmer_dict = {}

# k-mer weight (weight) = the number of times each kmer appears in any of the sequences
for i in tab:
    weight = len(tab[i])
    kmer_dict[i] = weight


seq_weights = {}
count = 0
sum_seq_weights = 0
med_list = []

start = time.time()
print "Getting sequence weights (avg runtime: 360 secs)..."
for i in dict1:
    seq = dict1[i]

    #make suffix tree of that seq
    stree = SuffixTree(seq)
    weight_seq = 0

    for j in kmer_dict:
        #if seq has the kmer, then add the kmer weight to its own seq. weight
        #seq weight = the sum of the kmer weights of all kmers that appear in the seq
        if (stree.hasSubstring(j)):
            weight_seq += kmer_dict[j]

    count += 1
    #print "seq weight made...", count
    seq_weights[i] = weight_seq
    sum_seq_weights += weight_seq
    med_list.append(weight_seq)
#print seq_weights
print "Done."
end = time.time()

prev_seq_weights = {}
for i in seq_weights: prev_seq_weights[i] = seq_weights[i]


avg_weight = sum_seq_weights/count
#print "Average seq weight:", avg_weight

#for getting threshold as the weights - average weight
for i in seq_weights:
    seq_weights[i] = seq_weights[i] - avg_weight



print "Getting unitig list..."

#a list of unitigs
unitigs = find_unitigs(overlap_list)
unitig_len_list = {}

############################## Attempt at making the genomes ###############
'''
genome = get_genome(dict1, unitigs)
#dict2 --- key: a unitig, value: list of overlap w. other unitigs
dict2 = {}
#dict3 --- key: a unitig, value: (max overlap value, index of unitig it has max
#                                  overlap with )
dict3 = {}
print "doing genome loop"
for i in genome:
    dict2[i] = []
    for j in genome:
        overlap = suffixPrefixMatch(i, j, 0)
        dict2[i].append(overlap)
    max_overlap = max(dict2[i])
    #print "Length of each dict[i], should be 130", len(dict2[i])
    dict3[i] = (max_overlap, dict2[i].index(max_overlap))


list_unitig_groups = []
new_list = []
rem_list = genome
cpy_rem_list = genome
num = 0
while cpy_rem_list != []:
    i = rem_list[num]
    cpy_rem_list.pop(num)
    new_list.append(i)
    over, idx = dict3[i]
    try:
        next_uni = genome[idx]
    except IndexError:
        continue
    if next_uni in new_list:
        list_unitig_groups.append(new_list)
        new_list = []
        num = idx
    else:
        #new_list.append(next_uni)
        num = idx
'''
#############################################################################
unitig_file = open("unitigs.txt", 'w')
for i in range(len(unitigs)):
    unitig_len_list[unitigs[i][0]] = len(unitigs[i])
    unitig_file.write("START UNITIG " + str(i+1) + " " + unitigs[i][0] + "\n")
    for j in range(1,len(unitigs[i])):
        unitig_file.write(unitigs[i][j] + "\n")
    unitig_file.write("END UNITIG "+ str(i+1) + "\n")



sum_length = 0
for i in unitig_len_list:
    sum_length += unitig_len_list[i]

avg_unitig_len = sum_length/len(unitig_len_list)


print "Contaminant i.d. list: "
contam_list = []

for i in seq_weights:
    if seq_weights[i] < 0:
        possible_contam = i
        if i in unitig_len_list:
            contam_list.append(i)

print contam_list
print "Estimated Percentage"
percent_contam = (len(contam_list)/ float(total_read_len))*100
print str(percent_contam),"%"

con_sum = 0
tot_weight_sum = 0
for i in prev_seq_weights:
    if i in contam_list:
        con_sum += prev_seq_weights[i]
        tot_weight_sum += prev_seq_weights[i]
    else:
        tot_weight_sum += prev_seq_weights[i]




print "Estimated Percentage (seq. weights)"
percent_contam = (con_sum/float(tot_weight_sum))*100
print str(percent_contam), "%"

uni_sum = 0
tot_uni_sum = 0
for i in unitig_len_list:
    if i in contam_list:
        uni_sum += 1
        tot_uni_sum += 1
    else:
        tot_uni_sum += 1

print "Estimated Percentage (unitigs)"
percent_contam = (uni_sum/float(tot_uni_sum))*100
print str(percent_contam), "%"


####################################################

