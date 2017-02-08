
'''
Computational Genomics Project: BaCon (Basic Contamination Detection)
Code: taken from Chris Simmons homework 5 and altered a bit

'''

def suffixPrefixMatch(str1, str2, min_overlap):
    ''' Returns length of longest suffix of str1 that is prefix of
        str2, as long as that suffix is at least as long as min_overlap. '''
    if len(str2) < min_overlap: return 0
    str2_prefix = str2[:min_overlap]
    str1_pos = -1
    while True:
        str1_pos = str1.find(str2_prefix, str1_pos + 1)
        if str1_pos == -1: return 0
        str1_suffix = str1[str1_pos:]
        if str2.startswith(str1_suffix): return len(str1_suffix)


def bestBuddy(fasta):
    '''
    Finds the right side best buddy given
    dictionary of sequence reads
    '''
    ret_list = []
    for id1 in fasta:
        id1_list = []
        str1 = fasta[id1]
        for id2 in fasta:
            if id2 != id1:
                str2 = fasta[id2]
                overlap_len = suffixPrefixMatch(str1,str2,40)
                if overlap_len != 0: #and best_buds != 0:

                    id1_list.append((id1,id2,overlap_len))
                    if len(id1_list) == 2:
                        if id1_list[0][2] > id1_list[1][2]:
                            id1_list.remove(id1_list[1])
                        elif id1_list[0][2] < id1_list[1][2]:
                            id1_list.remove(id1_list[0])
                        elif id1_list[0][2] == id1_list[1][2]:
                            id1_list = []
                        else:
                            continue

        #print "Current list: ", id1_list
        #for k in id1_list: ret_list.append(k)

        if len(id1_list) != 0:
            for k in id1_list: ret_list.append(k)
        #else:
        #    ret_list.append((id1, id1, len(fasta[id1])))
    return ret_list



def find_unitigs(overlap_list):
    '''
    Get list of unitigs
    :param overlap_list: tuple (str1, str2, overlap)
    '''
    left_best_buds = {}
    for (x,y,z) in overlap_list:
        if left_best_buds.get(y) == None or z > left_best_buds[y][1]:
            left_best_buds[y] = (x,z,0)
        elif z == left_best_buds[y][1]:
            left_best_buds[y] = (x,z,1)

    right_best_buds = {}
    for i in list(left_best_buds.items()):
        if (i[1][2] == 0):
            right_best_buds[i[1][0]] = (i[0], i[1][1])
            left_best_buds[i[0]] = (i[1][0], i[1][1])
        else:
            left_best_buds.pop(i[0])

    chars = list(right_best_buds.keys())
    unitig_list = []
    count = 0

    #quit()


    while chars != []:
        srt = chars.pop(0)
        unitig_list.append([srt])
        next = srt
        prev = srt

        while (right_best_buds.get(next) != None):
            if next == right_best_buds[next][0]: break
            next = right_best_buds[next][0]
            unitig_list[count].append(next)
            if right_best_buds.get(next) != None:
                try:
                    chars.remove(next)
                except ValueError:
                    #print "no next"
                    break


        while(left_best_buds.get(prev) != None):
            if prev == left_best_buds[prev][0]:
                unitig_list[count].insert(0,prev)
                break
            prev = left_best_buds[prev][0]
            unitig_list[count].insert(0, prev)
            try:
                chars.remove(prev)
            except ValueError:
                #print "no prev.."
                break

        count += 1

    return unitig_list

def get_genome(seqs, unitigs):
    '''
    Assemble Genome
    '''
    genomes = []
    for i in range(len(unitigs)):
        genome = seqs[unitigs[i][0]]
        for j in range(0, len(unitigs[i]) - 1):
            read = seqs[unitigs[i][j]]
            next = seqs[unitigs[i][j+1]]
            overlap = suffixPrefixMatch(read, next, 40)
            genome = ''.join([genome, next[overlap:]])

        genomes.append(genome)

    return genomes