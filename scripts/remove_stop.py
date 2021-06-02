#!/usr/bin/env python

import sys
import re
from Bio import SeqIO

new = []

for record in SeqIO.parse(sys.argv[1], "fasta"):

    full = str(record.seq)
    seq = re.search('([RHKDESTNQCUGPAILMFWYV]+?)X', full)
    if seq != None:
        seq = seq.group(1)
    else:
        seq = full
    name = '>{}'.format(record.id)

    new.append(name)
    new.append(seq)

thefile = open(sys.argv[2], 'w')
for item in new:
    thefile.write('{}\n'.format(item))
thefile.close()


