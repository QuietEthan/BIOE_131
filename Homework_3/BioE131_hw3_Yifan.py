#! /usr/bin/env python3
#####Plz run this script in python 3
##################
#####The Information Content Calculator
#####Yifan Ethan Chen
#####BioE 131 HW3
#####10/15/17 Cal

import sys, math, heapq
from sys import exit
from sys import argv

stockfile_name = sys.argv[1]
stockfile = open(stockfile_name) 
stockcontents = stockfile.readlines()

######function for calculating entropy
def cal_entropy(dist, dist_key):
    entropy = 0
    for key_name in dist[dist_key].keys():
        key_prob = dist[dist_key][key_name]
        if key_prob != 0:
            entropy -= key_prob*math.log2(key_prob)
    return entropy

#####function to find the key
def find_key(input_list, input_dist):
    rep_list = []
    for value in input_list:
        for key in input_dist: 
            if input_dist[key] == value:
                if key not in rep_list:
                    rep_list.append(key)
                    break
    try:
        for rep_ele in rep_list:
            ele1, ele2 = rep_ele
            print (str(ele1)+','+str(ele2))
    except:
        print(*rep_list, sep='\n')
    
#####set up
seq_num = 0
column_seq = {}    #i: [A,T....]
column_dist = {}   #i: A,T,G,C
column_entropy = {}    #i: entropy
colpair_dist = {}    #(i,j): AA ...
colpair_mi = {}    #(i,j): mi

#####deal with the imput file and store them as column
for line in stockcontents:
    if line[0] == '#':
        pass
    elif line[0:2] == '//':
        break
    else:
        line = line.rstrip()
        line = line.replace('.', '-').upper()
        line = line.replace('U', 'T')
        line_list = line.split()
        tg_name = line_list[0]
        tg_seq = line_list[1]
        if seq_num == 0:
            for i in range(len(tg_seq)):
                column_seq.setdefault(i, [tg_seq[i]])
        else:
            for i in range(len(tg_seq)):
                column_seq[i].append(tg_seq[i])
        seq_num += 1

#####calculate the probability and entropy for each column 
for i in column_seq.keys():
    a_num = column_seq[i].count('A')
    t_num = column_seq[i].count('T')
    c_num = column_seq[i].count('C')
    g_num = column_seq[i].count('G')
    dot_num = column_seq[i].count('-')
    total_num = a_num+t_num+c_num+g_num+dot_num 
    column_dist.setdefault(i, {'A':a_num/total_num,'T':t_num/total_num,'G':g_num/total_num,'C':c_num/total_num,'dot':dot_num/total_num})
    entropy_i = cal_entropy(column_dist, i)
    column_entropy.setdefault(i, round(entropy_i,6))

#####calculate the probability p and the mutual information I for every pair of column i and j, i<j
for column_j in range(len(column_entropy.keys())):
    column_i = 0
    col_j_nt_list = column_seq[column_j]
    while column_i < column_j:
        AA_num, AT_num, AG_num, AC_num, AD_num, TA_num, TT_num, TG_num, TC_num, TD_num, GA_num, GT_num, GG_num, GC_num, GD_num, CA_num, CT_num, CG_num, CC_num, CD_num, DA_num, DT_num, DG_num, DC_num, DD_num = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
        col_i_nt_list = column_seq[column_i] 
        for i in range(seq_num):
            if col_i_nt_list[i] == 'A':
                if col_j_nt_list[i] == 'A':
                    AA_num += 1
                elif col_j_nt_list[i] == 'T':
                    AT_num += 1
                elif col_j_nt_list[i] == 'G':
                    AG_num += 1
                elif col_j_nt_list[i] == 'C':
                    AC_num += 1
                elif col_j_nt_list[i] == '-':
                    AD_num += 1
            elif col_i_nt_list[i] == 'T':
                if col_j_nt_list[i] == 'A':
                    TA_num += 1
                elif col_j_nt_list[i] == 'T':
                    TT_num += 1
                elif col_j_nt_list[i] == 'G':
                    TG_num += 1
                elif col_j_nt_list[i] == 'C':
                    TC_num += 1
                elif col_j_nt_list[i] == '-':
                    TD_num += 1
            elif col_i_nt_list[i] == 'G':
                if col_j_nt_list[i] == 'A':
                    GA_num += 1
                elif col_j_nt_list[i] == 'T':
                    GT_num += 1
                elif col_j_nt_list[i] == 'G':
                    GG_num += 1
                elif col_j_nt_list[i] == 'C':
                    GC_num += 1
                elif col_j_nt_list[i] == '-':
                    GD_num += 1
            elif col_i_nt_list[i] == 'C':
                if col_j_nt_list[i] == 'A':
                    CA_num += 1
                elif col_j_nt_list[i] == 'T':
                    CT_num += 1
                elif col_j_nt_list[i] == 'G':
                    CG_num += 1
                elif col_j_nt_list[i] == 'C':
                    CC_num += 1
                elif col_j_nt_list[i] == '-':
                    CD_num += 1
            elif col_i_nt_list[i] == '-':
                if col_j_nt_list[i] == 'A':
                    DA_num += 1
                elif col_j_nt_list[i] == 'T':
                    DT_num += 1
                elif col_j_nt_list[i] == 'G':
                    DG_num += 1
                elif col_j_nt_list[i] == 'C':
                    DC_num += 1
                elif col_j_nt_list[i] == '-':
                    DD_num += 1
        total_seq_num = seq_num+1
        colpair_dist.setdefault((column_i,column_j), {'AA':AA_num/seq_num,'AT':AT_num/seq_num,'AG':AG_num/seq_num,'AC':AC_num/seq_num,'AD':AD_num/seq_num,
                                                      'TA':TA_num/seq_num,'TT':TT_num/seq_num,'TG':TG_num/seq_num,'TC':TC_num/seq_num,'TD':TD_num/seq_num,
                                                      'GA':GA_num/seq_num,'GT':GT_num/seq_num,'GG':GG_num/seq_num,'GC':GC_num/seq_num,'GD':GD_num/seq_num,
                                                      'CA':CA_num/seq_num,'CT':CT_num/seq_num,'CG':CG_num/seq_num,'CC':CC_num/seq_num,'CD':CD_num/seq_num,
                                                      'DA':DA_num/seq_num,'DT':DT_num/seq_num,'DG':DG_num/seq_num,'DC':DC_num/seq_num,'DD':DD_num/seq_num})    
        entropy_ij = cal_entropy(colpair_dist, (column_i,column_j))
        mi_ij = column_entropy[column_i] + column_entropy[column_j] - entropy_ij
        colpair_mi.setdefault((column_i,column_j), round(mi_ij,6))
        column_i += 1

#####find the bottom 10 column 10 columns i with the lowest entropy H and the top 50 columns i,j with the highest mutual information
bottom_10_entropy = heapq.nsmallest(10, column_entropy.values())
top_50_mi = heapq.nlargest(50, colpair_mi.values())
find_key(bottom_10_entropy, column_entropy)
find_key(top_50_mi, colpair_mi)

sys.exit()