# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 15:05:35 2018

@author: cvg4
"""


import pylinac

from pylinac.core.roi import LowContrastDiskROI, HighContrastDiskROI, DiskROI, bbox_center
import numpy as np
from skimage import feature, filters, morphology, measure



class BWH_LV(pylinac.LasVegas):

    def __init__(self, filepath):

        super().__init__(filepath)
        self.image.array = np.array(self.image.array,dtype=np.float)
        self._processing()
        self._phantom_ski_region = None
        self.image.check_inversion(position=(0.1, 0.1))
        self.cnr_treshold = 2
        
    def _processing(self):
        #
        eps = 10
        
        # processing crops the image to get rid of the jaws. This makes the
        # detection of the phantom more reliable
        thresh = filters.threshold_otsu(self.image.array)
        bw = morphology.closing(self.image.array> thresh, morphology.square(3))
        edges = feature.canny(bw)
        
        labeled = measure.label(edges)
        regions = measure.regionprops(labeled, intensity_image=bw)
        
        area = [(reg.area,reg) for reg in regions]
        region = max(area)[1] 
        
        cropSize = [region.bbox[0]+eps,self.image.shape[0]-region.bbox[2]+eps,region.bbox[1]+eps,self.image.shape[1]-region.bbox[3]+eps]           
        
        self.image.crop(cropSize[0],'left')
        self.image.crop(cropSize[1],'right')
        self.image.crop(cropSize[2],'top')
        self.image.crop(cropSize[3],'bottom')
        
        
    def _determine_low_contrast(self):
        """Sample the detail contrast regions."""
        # create 4 ROIs on each side of the phantom to determine the average background
#        angles = np.array([-10, 80, 170, 260,22.5]) + self.phantom_angle - 4
#        dists = np.array([0.24, 0.24, 0.24, 0.24,0.2]) * self.phantom_radius
        
        angles = np.array([-10, 80, 170, 260]) + self.phantom_angle - 4
        dists = np.array([0.24, 0.24, 0.24, 0.24]) * self.phantom_radius
        
        bg_rois = []
        for dist, angle in zip(dists, angles):
            roi = LowContrastDiskROI(self.image, angle, self.phantom_radius*0.03, dist, self.phantom_center,
                                     0.05)
            bg_rois.append(roi)
        avg_bg = np.mean([roi.pixel_value for roi in bg_rois])

        # create X ROIs to sample the low contrast holes
        angles = np.array([77, 116, 134.5, 0, 13, 77, 142, 153, -21, -29, -107, 182, 174, -37, -55, -105, 206.5, 189.5, -48.1, -67.8]) + self.phantom_angle
        dists = np.array([0.107, 0.141, 0.205, 0.179, 0.095, 0.042, 0.097, 0.178, 0.174, 0.088, 0.024, 0.091, 0.179, 0.189, 0.113, 0.0745, 0.115, 0.191, 0.2085, 0.146]) * self.phantom_radius
        roi_radii = np.array([0.028, 0.028, 0.028, 0.016, 0.016, 0.016, 0.016, 0.016, 0.012, 0.012, 0.012, 0.012, 0.012, 0.007, 0.007, 0.007, 0.007, 0.007, 0.003, 0.003])
        rois = []
        for dist, angle, radius in zip(dists, angles, roi_radii):
            roi = LowContrastDiskROI(self.image, angle, self.phantom_radius*radius, dist, self.phantom_center,
                                     self.threshold, avg_bg)
            rois.append(roi)

        # normalize the threshold
        self.threshold *= max(roi.contrast_constant for roi in rois)
        for roi in rois:
            roi.contrast_threshold = self.threshold
        self.bg_rois = bg_rois
        self.lc_rois = rois
        
        # our results
        self.contrast = rois[3].contrast;
        
        # visible circles
        visCirc = [self.lc_rois[i].contrast_to_noise>self.cnr_treshold for i in [3,8,13,18]];

        try:
            self.smallestVisCirc = min([i for i, x in enumerate(visCirc ) if not x])
        except:
            self.smallestVisCirc = 4
                

        
        
    def excel_lines(self):
        """Generate lines you can cut and paste into the excel workbook "database"
        """
   
        
        out=""
        out+="%.2f\n"%(self.contrast*100)
        out+="%d\n"%self.smallestVisCirc

        for i in [0,2,1,3]:
            out+="%.2f\n"%self.bg_rois[i].pixel_value
            out+="%.2f\n"%self.bg_rois[i].std

        return out

