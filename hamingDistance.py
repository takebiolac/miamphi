#! /usr/bin/env python
#coding: utf-8

"""
calculate Haming Distance between insert sequences (argment)
"""

import sys
import numpy as np
from Bio import SeqIO

def compute_hamming_distance_of_insert(inputFastaFilepath):
    seq_mat = []
    taxon_arr = []

    seq = SeqIO.parse(inputFastaFilepath, "fasta")
    for record in seq:
        taxon_arr.append(record.id.split("_"))
        seq_mat.append(np.array(record.seq))

    taxon_arr = np.array(taxon_arr)
    seq_mat = np.array(seq_mat)
    def compute_HammingDistance(X):
        return (X[:, None, :] != X).sum(2)

    dist_mat = compute_HammingDistance(seq_mat)

    def make_isSame_matrix(categoryIndex, taxon_arr):
        numOfSpec = len(taxon_arr)
        categoryProduct_mat = np.ones((numOfSpec, numOfSpec))
        for i in range(categoryIndex + 1):
            category_mat = np.tile(taxon_arr[:,i], (numOfSpec, 1))
            categoryProduct_mat = categoryProduct_mat * (category_mat == category_mat.T)
        return categoryProduct_mat == 1
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
    compute_hamming_distance_of_insert(inputFilepath)