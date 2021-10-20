#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
from typing import Sequence
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """Read fasta file
    ---
    Read a fasta.gz file and make a generator of sequences in the fasta.gz:\n
    Arguments:
    - amplicon_file (fasta.gz file): take a fasta.gz file in arguments of main.
    - minseqlen (int): minimum lenght of sequence in fasta.gz file, get in arguments of main.\n
    Return:
    - sequences (generator): give a generator of sequences in the fasta.gz file.
    """
    # open the fasta.gz
    with gzip.open(amplicon_file, 'rt') as fasta_gz:
        seq = ""
        # select lines in the fasta.gz
        for line in fasta_gz:
            # check if the line start with a ">"
            if line[0] == ">":
                # check if the sequence have the minimal lenght
                if len(seq) >= minseqlen:
                    # return the sequence
                    yield seq
                # reset the seq and restart the loop
                seq = ""
                continue
            seq += line.strip() # eliminate the carriage return
        # return the last line if have the minimal lenght
        if len(seq) >= minseqlen:
            yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """De-replication of sequences
    ---
    Select sequencs with a minimal lenght and a minimal count and return a list
    ['sequences', 'count']\n
    Arguments:
    - amplicon_file (fasta.gz file): take a fasta.gz file in arguments of main.
    - minseqlen (int): minimum lenght of sequence in fasta.gz file, get in arguments of main.\n
    - mincount (int): minimum count of a sequence to take it.\n
    Return:
    - sequence and count (list): list with sequence (str) and his count (int)
    if the sequence have the minimum lenght and count.
    """
    dict_occurences = {}
    sequences = list(read_fasta(amplicon_file, minseqlen))
    # get every sequences generate by read_fasta
    for sequence in sequences:
        # verify if the number of sequences is superior to the minimal
        if sequences.count(sequence) > mincount:
            # add the sequence in the dictionnary if absent
            if sequence not in dict_occurences.keys():
                dict_occurences[sequence] = 1
            # count one more if the sequence is already in the dictionnary
            else:
                dict_occurences[sequence] += 1
    # change the key(sequence) and the item(count of that sequence) in a list and make the generator
    for key in sorted(dict_occurences, key = dict_occurences.get, reverse= True):
        list_occ = [key, dict_occurences[key]]
        yield list_occ

def common(lst1, lst2):
    """search common elements in two list
    """
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """Get chunks in a sequence
    ---
    Divide the sequence in number of chunk of the same lenght.\n
    Arguments:
    - sequence (str): sequence cut for have chunks.
    - chunk_size (int): length of desire chunks.\n
    Return:
    - chunks (list): list of chunks in the sequence
    """
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size]
              for i in range(0, len_seq, chunk_size)
                if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """Dictionnary in the fasta.gz
    ---
    Take a dictionnary with the k-mer at keys and the index of sequence
    where we can find this k-mer. E.g:\n
    ```Kmer_dict = { “ATTCCCAA”: [0],  “TTCCCAAG”: [0], “TCCCAAGT”: [0],  “CCCAAGTC”: [0,1],
    “CCAAGTCT”: [0, 1], “CAAGTCTT”: [1], “AAGTCTTC”: [1], “AGTCTTCC”: [1]}\n
    Arguments:
    - kmer_dict (dict): dictionnary with unique k-mer and his
    index of sequence where they are finds (can be empty).
    - sequence (str): sequence where we search k-mer.
    - id_seq (int): index of the sequence.
    - kmer_size (int): lenght of k-mer in the sequence.\n
    Return:
    - kmer_dict (dict): dictionnary with new k-mers and new id_seq for current k-mers.
    """
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict:
            kmer_dict.setdefault(kmer, [])
            kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict.setdefault(kmer, [])
            kmer_dict[kmer].append(id_seq)
    return kmer_dict

def search_mates(kmer_dict, sequence, kmer_size):
    """ Search mates of sequence
    ---
    Search the two best mates of a sequence based on the dictionnary of k-mers.\n
    Arguments:
    kmer_dict (dict): dictionnary with k-mers in non-chimerical sequences
    and there affiliations.
    - sequence (str): candidat sequence for know if it's parents.
    - kmer_size (int): lenght of k-mers (need to be the same in the kmer_dict).\n
    Return:
    - list_best_mates (list): list of tuples with the two parents
    for the sequence and their common k-mers.
    """
    list_mates = []
    list_best_mates = []
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict.keys():
            for parent in kmer_dict.get(kmer):
                list_mates.append(parent)
    best_mates = Counter(list_mates).most_common(2)
    for elts in best_mates:
        list_best_mates.append(elts[0])
    return list_best_mates

def detect_chimera(perc_identity_matrix):
    """Detect chimera
    ---
    Detect if a sequence is a chimera between to parents' sequences.\n
    Arguments:
    - perc_identity_matrix (matrix): matrix of identity of each chunk between chunk of parents.\n
    Return:
    - chimera (bool): Boolean of if the sequence is a chimera between the two parents.
    """
    std_identity = []
    sequence_1 = False
    sequence_2 = False
    for elts in perc_identity_matrix:
        std_identity.append(statistics.stdev(elts))
    if statistics.mean(std_identity) > 5:
        for elts in perc_identity_matrix:
            if elts[0] > elts[1]:
                sequence_1 = True
                continue
            elif elts[0] < elts[1]:
                sequence_2 = True
                continue
    if sequence_1 is True and sequence_2 is True:
        return True
    else:
        return False

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Remove chimerical sequences
    ---
    Search for a sequence if it's a chimera of two others sequences
    (first two sequences are defined as non-chimerical by default).\n
    Arguments:
    - amplicon_file (fasta.gz file): take a fasta.gz file in arguments of main.
    - minseqlen (int): minimum lenght of sequence in fasta.gz file, get in arguments of main.\n
    - mincount (int): minimum count of a sequence to take it.
    - chunck_size (int): size of chunck in arguments of main.
    - kmer_size (int): size of k-mer in arguments of main.\n
    Return:
    - generator of sequences non-chimerical (list): list with the sequence
    and his number of occurence.
    """
    list_seq_accept = []
    kmer_dict = {}
    list_chunks = []
    generator_seq = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    seq_temp = next(generator_seq)
    list_seq_accept.append(seq_temp)
    list_chunks.append(get_chunks(seq_temp[0], chunk_size))
    get_unique_kmer(kmer_dict, seq_temp, list_seq_accept.index(seq_temp), kmer_size)
    yield seq_temp
    seq_temp = next(generator_seq)
    list_seq_accept.append(seq_temp)
    list_chunks.append(get_chunks(seq_temp[0], chunk_size))
    get_unique_kmer(kmer_dict, seq_temp, list_seq_accept.index(seq_temp), kmer_size)
    yield seq_temp
    for seq_temp in generator_seq:
        parents = search_mates(kmer_dict, seq_temp[0], kmer_size)
        if len(parents) == 2:
            perc_matrix_identity = []
            list_chunk_seq = get_chunks(seq_temp[0], chunk_size)
            list_chunk_p1 = list_chunks[parents[0][0]]
            list_chunk_p2 = list_chunks[parents[1][0]]
            for index, chunk in enumerate(list_chunk_seq):
                perc_matrix_identity_temp = []
                perc_matrix_identity_temp.append(get_identity(nw.global_align(chunk,
                list_chunk_p1[index], gap_open=-1, gap_extend=-1,
                matrix=os.path.abspath(os.path.join(os.path.dirname(__file__), "MATCH")))))
                perc_matrix_identity_temp.append(get_identity(nw.global_align(chunk,
                list_chunk_p2[index], gap_open=-1, gap_extend=-1,
                matrix=os.path.abspath(os.path.join(os.path.dirname(__file__), "MATCH")))))
                perc_matrix_identity.append(perc_matrix_identity_temp)
            chimera = detect_chimera(perc_matrix_identity)
            if chimera is False:
                list_seq_accept.append(seq_temp)
                list_chunks.append(get_chunks(seq_temp[0], chunk_size))
                get_unique_kmer(kmer_dict, seq_temp, list_seq_accept.index(seq_temp), kmer_size)
                yield seq_temp
        else:
            list_seq_accept.append(seq_temp)
            list_chunks.append(get_chunks(seq_temp[0], chunk_size))
            get_unique_kmer(kmer_dict, seq_temp, list_seq_accept.index(seq_temp), kmer_size)
            yield seq_temp

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Select OTU in fasta.gz file
    ---
    Get OTU based and the abundance of sequences (count of unique sequence)
    and their similarity (<97%).\n
    Arguments:
    - amplicon_file (fasta.gz file): take a fasta.gz file in arguments of main.
    - minseqlen (int): minimum lenght of sequence in fasta.gz file, get in arguments of main.\n
    - mincount (int): minimum count of a sequence to take it.
    - chunck_size (int): size of chunck in arguments of main.
    - kmer_size (int): size of k-mer in arguments of main.\n
    Return:
    - list_otu (list): list of more abundant OTU in the fasta.gz file.
    """
    otu_list = []
    sequences = list(chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size))
    otu_list.append(sequences[0])
    for sequence in sequences:
        for otu in otu_list:
            alignement = nw.global_align(sequence[0], otu[0], gap_open=-1, gap_extend=-1,
            matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            identity = get_identity(alignement)
            if identity < 97.00:
                otu_list.append(sequence)
    return otu_list

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(otu_list, output_file):
    """Write list OTU in fasta
    ---
    Give a fasta file like :\n
    >OTU_1 occurrence:1000\n
    ATGGTCA...
    >OTU_2 occurrence:900\n
    GTAACTA...
    etc...\n
    Arguments:
    - OTU_list (list): list of OTU generate by the abundance_greedy_clustering function.
    - output_file (int): name of the fasta file get in the main programm\n
    Return:
    - output_file.fna (fasta file): fasta with the OTU and their occurrences.
    """
    with open(f"{output_file}", "w", encoding = 'utf-8') as fasta:
        # extract number of contig and contig
        for index,otu in enumerate(otu_list):
            # write a fasta file
            fasta.write(">OTU_" + str(index+1) + " occurrence:" + str(otu[1]) +
            "\n" + fill(otu[0]) + "\n")


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    print("Computing OTU and remove chimera ...")
    otu_list = abundance_greedy_clustering(args.amplicon_file,
    args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
    print("Write OTU fasta")
    write_OTU(otu_list, args.out_file)

if __name__ == '__main__':
    main()
