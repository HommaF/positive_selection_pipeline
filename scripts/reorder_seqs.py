#!/usr/bin/env python

import sys
import re
from Bio import SeqIO

seqs = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))

new = []

for i in seqs:
    name = re.search('([A-Za-z]+)[0-9]', i).group(1)

    if name == "Solyc":
        ref_id = '>{}'.format(i)
        ref_seq = seqs[i].seq

    else:
        new.append('>{}'.format(i))
        new.append(seqs[i].seq)

new.insert(0, ref_seq)
new.insert(0, ref_id)


thefile = open(sys.argv[2], 'w')
for item in new:
    thefile.write('{}\n'.format(item))
thefile.close()
