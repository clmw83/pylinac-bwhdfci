# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 11:24:51 2017

@author: clw34
"""

from pylinac.core import image
from pylinac.core.geometry import Point
from pylinac.core.mask import filled_area_ratio, bounding_box
from pylinac.core.profile import SingleProfile

import numpy as np
from scipy import ndimage

import os,os.path
import matplotlib.pyplot as plt

class IsoCubeSet():

    def __init__(self,imagedir=None):
        self.MV_offsets=None
        self.kV_offsets=None
        
        if imagedir is not None:
            self.scan_dir(os.path.abspath(imagedir))
    
    def scan_dir(self,imagedir,plot=False):
        self.MV_offsets={}
        self.kV_offsets={}
        self.imagedir=os.path.abspath(imagedir)        
        l=os.listdir(imagedir)
        
        for f in l:
            ci = CubeImage(os.path.join(imagedir,f))
            if plot:
                ci.plot()
                
            bb=ci._find_bb()
            
            if plot:
                caxx=(-1*ci.metadata.XRayImageReceptorTranslation[0]*ci.metadata.RTImageSID/1000-ci.metadata.RTImagePosition[0])/(ci.metadata.ImagePlanePixelSpacing[0] * ci.metadata.RTImageOrientation[0])
                caxy=(-1*ci.metadata.XRayImageReceptorTranslation[1]*ci.metadata.RTImageSID/1000-ci.metadata.RTImagePosition[1])/(ci.metadata.ImagePlanePixelSpacing[1] * ci.metadata.RTImageOrientation[4])
                plt.plot(bb.x,bb.y,'x',color='blue')
                plt.plot(ci.cax.x,ci.cax.y,'+',color='red')
                plt.plot(caxx,caxy,'+',color='green')
        
            angle=ci.metadata.ExposureSequence[0].GantryAngle
            around=90*int(round(ci.metadata.ExposureSequence[0].GantryAngle/90))
            if abs(angle-around) > 2:
                print("ERROR!",angle,around)
            around = around % 360  
            energy=ci.metadata.ExposureSequence[0].KVP
    
            dx,dy=self.pix2pos(bb)
            #dx=ci.metadata.XRayImageReceptorTranslation[0]+(ci.metadata.RTImagePosition[0]+bb.x*ci.metadata.ImagePlanePixelSpacing[0] * ci.metadata.RTImageOrientation[0])*1000/ci.metadata.RTImageSID
            #dy=ci.metadata.XRayImageReceptorTranslation[1]+(ci.metadata.RTImagePosition[1]+bb.y*ci.metadata.ImagePlanePixelSpacing[1] * ci.metadata.RTImageOrientation[4])*1000/ci.metadata.RTImageSID
        
            # TODO: need some sanity checking here.
            if energy>1000:
                self.MV_offsets[around]=(dx,dy)
            else:
                self.kV_offsets[around]=(dx,dy)                
    
   
    def excel_lines(self):
        """ print out excel lines for the BWH/DFCI excel workbooks"""
        
        print("kV:")
        out=[None]*8
        for angle,delta in self.kV_offsets.items():
            i=angle//90
            out[i*2]=delta[1]
            out[i*2+1]=delta[0]
        for j in out:
            print(j)
    
        print("\nMV:")
        out=[None]*8
        for angle,delta in self.MV_offsets.items():
            i=angle//90
            out[i*2]=delta[1]
            out[i*2+1]=delta[0]
        for j in out:
                print(j)
    
    

class CubeImage(image.DicomImage):

    @staticmethod
    def _is_symmetric(logical_array):
        """Whether the binary object's dimensions are symmetric, i.e. a perfect circle. Used to find the BB."""
        ymin, ymax, xmin, xmax = bounding_box(logical_array)
        y = abs(ymax - ymin)
        x = abs(xmax - xmin)
        if x > max(y * 1.05, y + 3) or x < min(y * 0.95, y - 3):
            return False
        return True
    
    def _is_modest_size(self,logical_array,shape):
        """Decide whether the ROI is roughly the size of a BB; not noise and not an artifact. Used to find the BB."""
        rad_field_area = (shape[0]*shape[1])
        return rad_field_area * 0.003 < np.sum(logical_array) < rad_field_area * 0.5
    
    @staticmethod
    def _is_round(logical_array):
        """Decide if the ROI is circular in nature by testing the filled area vs bounding box. Used to find the BB."""
        expected_fill_ratio = np.pi / 4
        actual_fill_ratio = filled_area_ratio(logical_array)
        return expected_fill_ratio * 1.2 > actual_fill_ratio > expected_fill_ratio * 0.8
    

    def pix2pos(self,p):
        """ 
        Take in an x/y point in pixels in the image, and return a point with the coordinates in mm based on the RT graticule/coordinate system
        Returns
        -------
        float,float
            The x/y position in 
        """       
        xmm=self.metadata.XRayImageReceptorTranslation[0]+(self.metadata.RTImagePosition[0]+p.x*self.metadata.ImagePlanePixelSpacing[0] * self.metadata.RTImageOrientation[0])*1000/self.metadata.RTImageSID
        ymm=self.metadata.XRayImageReceptorTranslation[1]+(self.metadata.RTImagePosition[1]+p.y*self.metadata.ImagePlanePixelSpacing[1] * self.metadata.RTImageOrientation[4])*1000/self.metadata.RTImageSID
        return (xmm,ymm)
    
    def _find_bb(self):
        """Find the BB within the radiation field. Dervived from pylinac WL test.  Looks at the central 60x60 pixels and finds
        a bb

        Returns
        -------
        Point
            The weighted-pixel value location of the BB.
        """
        span=30
        bbox=[int(self.cax.y-span),int(self.cax.x-span),int(self.cax.y+span),int(self.cax.x+span)]
        subim=np.array(self.array[bbox[0]:bbox[2],bbox[1]:bbox[3]],dtype=np.float)

        
        hmin, hmax = np.percentile(subim, [10, 100.0])
        spread = hmax - hmin
        max_thresh = hmax
        lower_thresh = hmax - spread*.95
        # search for the BB by iteratively lowering the low-pass threshold value until the BB is found.
        found = False
        while not found:
            try:
                binary_arr = np.logical_and((subim <= max_thresh), (subim >= lower_thresh))
                labeled_arr, num_roi = ndimage.measurements.label(binary_arr)
                roi_sizes, bin_edges = np.histogram(labeled_arr, bins=num_roi + 1)
                bw_bb_img = np.where(labeled_arr == np.argsort(roi_sizes)[-2], 1, 0)
                if not self._is_round(bw_bb_img):
                    raise ValueError
                if not self._is_modest_size(bw_bb_img,subim.shape):
                    raise ValueError
                if not self._is_symmetric(bw_bb_img):
                    raise ValueError
            except (IndexError, ValueError):
                lower_thresh += 0.05 * spread
                if lower_thresh > hmax-spread*.5:
                    raise ValueError("Unable to locate the BB. Make sure the field edges do not obscure the BB and that there is no artifacts in the images.")
            else:
                found = True

        # determine the center of mass of the BB
        
        x_arr = np.abs(np.average(subim*bw_bb_img, axis=0))
        x_com = SingleProfile(np.array(x_arr,dtype=np.float)).fwxm_center(interpolate=True)
        y_arr = np.abs(np.average(subim*bw_bb_img, axis=1))
        y_com = SingleProfile(y_arr).fwxm_center(interpolate=True)              
        return Point(x_com+bbox[1], y_com+bbox[0])

    