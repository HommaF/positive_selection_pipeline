#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import csv
import sys
import glob
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight


# Script to identify residues under selection in a protein of interest and look up the alternative residues
# sys.argv[1]: locus name
# sys.argv[2]: path to to mlc_8 file
# sys.argv[3]: path to codon_guided_msa
# sys.argv[4]: path to tmp outfile
# sys.argv[5]: path to output file

name = sys.argv[1]
mlc_file = sys.argv[2]
msa = sys.argv[3]
tmp_outfile =sys.argv[4]
outfile = sys.argv[5]

dic_codons = pd.read_csv('scripts/codons.csv', sep=',', engine='python', index_col=0)
dic_codons = pd.Series(dic_codons.aa.values, index=dic_codons.index).to_dict()

new_file = []

start = 0
new = list()
loci = list()


with open(mlc_file) as f:
    for line in f:

        begin = re.match('Bayes Empirical Bayes \(BEB\)', line)
        end = re.match('The grid', line)

        if begin != None:
            start = 1

        if end != None:
            start = 0

        if start == 1:

            sign = re.match(' +[0-9]+ [A-Z\-]{1} +([0-9\.]{5})', line)
            if sign != None:

                sign_lv = float(sign.group(1))

                if sign_lv >= 0.8:

                    entry = []

                    pos_aa = int(re.match(' +([0-9]+) +[A-Z\-]{1}?', line).group(1))
                    ref_aa = re.match(' +[0-9]+ {1}?([A-Z\-]{1}?)', line).group(1)

                    loci.append(pos_aa-1)

                    entry.append(pos_aa)
                    entry.append(ref_aa)
                    entry.append(sign_lv)

                    new.append(entry)


linear_msa = []
seq = ''


for record in SeqIO.parse(msa, "fasta"):
    linear_msa.append('>{}'.format(record.id))
    linear_msa.append(str(record.seq))


alt_seq = []
for i in range(int(len(loci))):
    alt_seq.append([])

codon_list = []
for i in range(int(len(loci))):
    codon_list.append([])
    
header = []
numbers = []
ref_seq = []

ref_hit = 0

for line in linear_msa:

    if ref_hit < 2:
        ref_hit += 1

    if ref_hit == 2:
        ref_hit = 3

        gap_counter = 0

        for i in range(int(len(line)/3)):

            codon = line[0+3*i:3+3*i]

            if codon  == '---':
                gap_counter += 1

            elif i in loci:
                
                if gap_counter == 0:
                    header.append(str(i+1))
                    numbers.append(i+1)

                elif gap_counter != 0:
                    header.append('{}(-{})'.format(i+1, gap_counter))
                    numbers.append(i+1-gap_counter)


                ref_seq.append(dic_codons[codon.upper()])

    elif (ref_hit == 3) & (line[:1] != '>'):


        counter = 0

        for i in range(int(len(line)/3)):

            codon = line[0+3*i:3+3*i]

            if i in loci:

                if dic_codons[codon.upper()] not in alt_seq[counter]:
                    if dic_codons[codon.upper()] != ref_seq[counter]:
                        alt_seq[counter].append(dic_codons[codon.upper()])

                elif dic_codons[codon.upper()] in alt_seq[counter]:
                    if codon.upper() not in codon_list[counter]:
                        codon_list[counter].append('*')


                if codon.upper() not in codon_list[counter]:
                    if dic_codons[codon.upper()] != ref_seq[counter]:
                        codon_list[counter].append(codon.upper())

                counter += 1

final_one = []
final_two = []

for i in range(int(len(header))):

    temp_one = []
    temp_two = []

    ref_aa_temp = None

    temp_one.append(ref_seq[i])
    temp_one.append(header[i])
    temp_one.append(''.join(alt_seq[i]))


    final_one.append(''.join(temp_one))


numbers = [str(x) for x in numbers]

final_one = ';'.join(final_one)
numbers = ';'.join(numbers)
temp = []

temp.append(name)
temp.append(final_one)
temp.append(numbers)

new_file.append(temp)

data = pd.DataFrame(new_file, columns=['seq_name', 'residues under positive selection', 'pos_under_selection'])


for i in data.index:
    target = data.loc[i, 'seq_name']
    os.system("grep -A 2 'Parameters in M8 (beta&w>1):' {} | grep ' (p1 =' > {}".format(mlc_file, tmp_outfile))

    with open(tmp_outfile) as f:
        for lines in f:
            rho = re.search('p1 \=[ ]+([0-9\.]+)\)', lines)
            omega = re.search('w =[ ]+([0-9\.]+)', lines)

            if omega != None:
                omega = float(omega.group(1))

            if rho != None:
                rho = float(rho.group(1))


        data.loc[i, 'omega'] = omega
        data.loc[i, 'rho'] = rho
           
data.to_csv(outfile, sep='\t', index=None)

os.remove(tmp_outfile)
