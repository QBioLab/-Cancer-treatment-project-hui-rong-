# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 23:21:37 2021

@author: chena
"""

import pandas as pd

def to_hcg(prefix):
    count = pd.read_csv('{}.hc.txt'.format(prefix),sep='\t',header=None)
    count = count[count[1].notna()]
    count[0] = count[1]
    count.to_csv('{}.hcg.txt'.format(prefix),sep='\t',header=False,index=False)
    
prefixs = ['Fluc_D1','Fluc_D2','Fluc_D3','Fluc_L1','Fluc_L2','Fluc_L3',
           'mINF_D1','mINF_D3','mINF_D4','mINF_L2','mINF_L3','mINF_L4']
for p in prefixs:
    to_hcg(p)