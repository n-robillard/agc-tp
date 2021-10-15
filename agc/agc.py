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
    Select sequencs with a minimal lenght and a minimal count and return a list ['sequences', 'count']\n
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

def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """"""
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

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Selet OTU in fasta.gz file
    ---
    Get OTU based and the abundance of sequences (count of unique sequence) and their similarity (<97%).\n
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
    sequences = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    otu_list.append(sequences[0])
    for sequence in sequences:
        for otu in otu_list:
            alignement = nw.global_align(sequence[0], otu[0], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            identity = get_identity(alignement)
            if identity < 97.00:
                otu_list.append(sequence)
    return otu_list

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
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
    with open(f"{output_file}", "w") as fasta:
        # extract number of contig and contig
        for index,otu in enumerate(OTU_list):
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
    abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)


if __name__ == '__main__':
    main()
