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

modeltesting = sys.argv[1]
summary_models = sys.argv[2]
summary_tests = sys.argv[3]

info = list()
new = list()

with open(modeltesting) as f:

    for lines in f:
        print(lines)
        model = re.search('mlc_([0-9])_', lines)
        name = re.search('mlc_[0-9]_(Solyc[0-9]+g[0-9\.]+)\.txt', lines)
        df = re.search('np:.*?([0-9]*)\):', lines)
        ln = re.search('\): +([0-9\.\-]*) ', lines)
        test = re.search(":\t([0-9\.\-]*)\t", lines)


        if name != None:


            info.append(str(name.group(1)))
            info.append(int(model.group(1)))
            info.append(int(df.group(1)))
            info.append(float(ln.group(1)))


for i in range(int(len(info)/4)):
    new.append(info[0+4*i:4+4*i])


df_sum = pd.DataFrame(new, columns=['name', 'model', 'df', 'logL'])
df_sum = df_sum.sort_values(by=['name', 'model'])

comp = list()
new = list()

for i in df_sum.name.unique():

    df_temp = df_sum[df_sum.name == i]

    ln1 = df_temp.logL[df_temp.model == 1].item()
    ln2 = df_temp.logL[df_temp.model == 2].item()
    ln7 = df_temp.logL[df_temp.model == 7].item()
    ln8 = df_temp.logL[df_temp.model == 8].item()

    df1 = df_temp.df[df_temp.model == 1].item()
    df2 = df_temp.df[df_temp.model == 2].item()
    df7 = df_temp.df[df_temp.model == 7].item()
    df8 = df_temp.df[df_temp.model == 8].item()

    dln1_2 = ln2 - ln1
    dln7_8 = ln8 - ln7


    ddf1_2 = df2 - df1
    ddf7_8 = df8 - df7


    pchisq1_2 = stats.chi2.sf(2*dln1_2, ddf1_2)
    pchisq7_8 = stats.chi2.sf(2*dln7_8, ddf7_8)


    comp.append(i)
    comp.append(pchisq1_2)
    comp.append(pchisq7_8)



for i in range(int(len(comp)/3)):
    new.append(comp[0+3*i:3+3*i])

df_ana = pd.DataFrame(new, columns=['name', 'M1-M2', 'M7-M8'])

df_sum.to_csv(summary_models, index=None)
df_ana.to_csv(summary_tests, index=None)


#df_sum.to_csv('summary_selection_models.csv', index=None)
#df_ana.to_csv('summary_model_testing.csv', index=None)


#os.system('mv *codeml.ctl ./model_testing_CODEML/')
#os.system('mv *codonAlignment.txt ./model_testing_CODEML/')
