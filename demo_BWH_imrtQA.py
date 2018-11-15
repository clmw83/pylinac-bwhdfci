#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 23:57:00 2018

@author: christianguthier
"""

import IMRTQA as QA
#
folder = r'C:\Users\cvg4\Dropbox (Partners HealthCare)\QA\IMRT-QA\Gamma\Pat01'
gCalc,rtplan= QA.scanFolder(folder)

for key,gc in gCalc.items():
    print(f'Analyzing Beam: {key}')
    gc.analyzeAll()
    
QA.publish_pdf(gCalc,rtplan,path=folder,filename = 'Report.pdf')    

