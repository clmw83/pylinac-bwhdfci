# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 13:30:19 2018

@author: jt108
"""

from BWHCatPhan import CatPhan504


cbct_folder = "C:\Git\pylinac-bwhdfci\Files\CatPhan\88880181982"

cbct = CatPhan504(cbct_folder)
cbct.analyze()
cbct.plot_analyzed_image()
print(cbct.results())

pdf_fn = "result.pdf"
cbct.publish_pdf(pdf_fn)
