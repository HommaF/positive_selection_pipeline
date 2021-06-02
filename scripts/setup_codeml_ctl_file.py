#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
import glob
import os
import subprocess
import re
from scipy import stats
import string

#####################################################################################
### read in files and run exonerate on them using cds and genomic regions ###########
#####################################################################################

msa = sys.argv[1]
tree = sys.argv[2]
locus = sys.argv[3]

models = (1, 2, 7, 8)

for model in models:
    codeml_ctl = list()

    with open("codeml/codeml_default.ctl", "r") as f:
        for line in f:
            line = line.strip()

            if line[:7] == "seqfile":
                line = "seqfile = {}".format(msa)
                codeml_ctl.append(line)

            elif line[:8] == "treefile":
                line = "treefile = {}".format(tree)
                codeml_ctl.append(line)

            elif line[:7] == "outfile":
                line = "outfile = mlc_files/mlc_{}_{}.txt".format(model,locus)
                codeml_ctl.append(line)

            elif line[:7] == "NSsites":
                line = "NSsites = {}".format(model)
                codeml_ctl.append(line)

            else:
                codeml_ctl.append(line)

    name = 'codeml/{}/ctl_files/{}_{}_codeml.ctl'.format(locus, locus, model)
    thefile = open(name, "w")
    for item in codeml_ctl:
        thefile.write("{}\n".format(item))
    thefile.close()

