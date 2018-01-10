#! /usr/bin/env python
#coding: utf-8

"""
calculate pairwise Haming Distance between insert sequences (argment)
"""

import sys
import numpy as np
from Bio import SeqIO
from Bio import pairwise2

def make_distance_table(inputFastaFilepath):
    def compute_hamming_distance_between_2_seq(seq_arr):
        assert len(seq_arr) == 2
        seq1 = seq_arr[0]
        seq2 = seq_arr[1]
        alignments = pairwise2.align.globalms(seq1, seq2, 1,-1,-1,-1)
        score = alignments[0][2]
        begin = alignments[0][3]
        end = alignments[0][4]

        length = end - begin
        distance = int((length - score) / 2)
        return distance

    taxon_lst = []
    seq_lst = []

    seq = SeqIO.parse(inputFastaFilepath, "fasta")
    for record in seq:
        taxon_lst.append(record.id.split("_"))
        seq_lst.append(record.seq)

    numOfSepc = len(seq_lst)
    seq_mat = np.tile(seq_lst, (numOfSepc, 1))

    stackedSeq_mat = np.stack([seq_mat, seq_mat.T], axis=2)
    dist_mat = np.apply_along_axis(compute_hamming_distance_between_2_seq, 2, stackedSeq_mat)

    def make_isSame_matrix(categoryIndex, taxon_arr):
        numOfSpec = len(taxon_arr)
        categoryProduct_mat = np.ones((numOfSpec, numOfSpec))
        for i in range(categoryIndex + 1):
            category_mat = np.tile(taxon_arr[:,i], (numOfSpec, 1))
            categoryProduct_mat = categoryProduct_mat * (category_mat == category_mat.T)
        return categoryProduct_mat == 1

    taxon_arr = np.array(taxon_lst)
    isSameFamily_mat = make_isSame_matrix(0, taxon_arr)
    isSameGenus_mat = make_isSame_matrix(1, taxon_arr)
    isSameSpec_mat = make_isSame_matrix(2, taxon_arr)
    
    def output_distance_table(dist_mat, isSameCategory_mat, maxDistance):
        for i in range(maxDistance):
            print(i, int(np.sum((dist_mat == i) * ~isSameCategory_mat) / 2))
        print(i, "<", int(np.sum((dist_mat > i) * ~isSameCategory_mat) / 2), "\n")
    print("Species")
    output_distance_table(dist_mat, isSameSpec_mat, 6)
    print("Genus")
    output_distance_table(dist_mat, isSameGenus_mat, 6)
    print("Family")
    output_distance_table(dist_mat, isSameFamily_mat, 6)

if __name__=="__main__":
    argvs = sys.argv
    argc = len(argvs)
    if (argc != 2):
        print('Usage: # python {} inputFilepath'.format(argvs[0]))
        quit()
    
    inputFilepath = argvs[1]
    make_distance_table(inputFilepath)