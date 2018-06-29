# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 15:05:35 2018

@author: cvg4
"""


import pylinac

from pylinac.core.roi import LowContrastDiskROI, HighContrastDiskROI, DiskROI, bbox_center, Point
import numpy as np
from skimage import feature, filters, morphology, measure
from scipy.ndimage.interpolation import rotate
from scipy.spatial import ConvexHull



class BWH_LV(pylinac.LasVegas):

    def __init__(self, filepath):

        super().__init__(filepath)
        self.image.array = np.array(self.image.array,dtype=np.float)
        self._processing()
        self._phantom_ski_region = None
        self.image.check_inversion(position=(0.1, 0.1))
        self.cnr_treshold = 2
        print('Test')
        self._phantom_detection()
      
    @property
    def phantom_center(self):
        return self._phantom_center
    
    @phantom_center.setter
    def phantom_center(self, value):
        self._phantom_center = Point(value[1],value[0],0.)
        
    @property
    def phantom_angle(self):
        return self._phantom_angle
    
    @phantom_angle.setter
    def phantom_angle(self, value):
        self._phantom_angle= value

        
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
        
        cropSize = [region .bbox[1]+eps,self.image.shape[1]-region.bbox[3]+eps,region.bbox[0]+eps,self.image.shape[0]-region.bbox[2]+eps]
        self.image.crop(cropSize[0],'left')
        self.image.crop(cropSize[1],'right')
        self.image.crop(cropSize[2],'top')
        self.image.crop(cropSize[3],'bottom')
        
    
    def _phantom_detection(self):
        # generate phatnom via otsu method and morpological closing
        image = self.image.array
        thresh = filters.threshold_otsu(image)
        bw = morphology.closing(image> thresh, morphology.square(3))
        # get points along the edges
        points = feature.corner_peaks(feature.corner_harris(bw),min_distance = 5)
        # estimation of the rectangular lv phantom via convex hull
        # get the convex hull for the points
        hull_points = points[ConvexHull(points).vertices]
        # calculate edge angles
        edges = np.zeros((len(hull_points)-1, 2))
        edges = hull_points[1:] - hull_points[:-1]
        # find angles between points at the hull
        angles = np.zeros((len(edges)))
        angles = np.arctan2(edges[:, 1], edges[:, 0])
        angles = np.abs(np.mod(angles, np.pi/2))
        angles = np.unique(angles)
        # rotation matrix
        rotations = np.vstack([
        np.cos(angles),
        np.cos(angles-np.pi/2),
        np.cos(angles+np.pi/2),
        np.cos(angles)]).T
    
        rotations = rotations.reshape((-1, 2, 2))
        
        # apply rotations to the hull
        rot_points = np.dot(rotations, hull_points.T)

        # find the bounding points
        min_x = np.nanmin(rot_points[:, 0], axis=1)
        max_x = np.nanmax(rot_points[:, 0], axis=1)
        min_y = np.nanmin(rot_points[:, 1], axis=1)
        max_y = np.nanmax(rot_points[:, 1], axis=1)
        
        # find the box with the best area
        areas = (max_x - min_x) * (max_y - min_y)
        best_idx = np.argmin(areas)
        
        # return the best box
        x1 = max_x[best_idx]
        x2 = min_x[best_idx]
        y1 = max_y[best_idx]
        y2 = min_y[best_idx]
        r = rotations[best_idx]
        
        self.edges = np.zeros((4, 2))
        self.edges[0] = np.dot([x1, y2], r)
        self.edges[1] = np.dot([x2, y2], r)
        self.edges[2] = np.dot([x2, y1], r)
        self.edges[3] = np.dot([x1, y1], r)
        
        self.phantom_center = np.mean(self.edges,axis=0)
        self.phantom_angle  = np.rad2deg(angles[best_idx]);
        
        if self.phantom_angle>45:
            self.phantom_angle -=90
        

        
        
        
        

        
        
    def _determine_low_contrast(self):
        """Sample the detail contrast regions."""
        # create 4 ROIs on each side of the phantom to determine the average background
#        angles = np.array([-10, 80, 170, 260,22.5]) + self.phantom_angle - 4
#        dists = np.array([0.24, 0.24, 0.24, 0.24,0.2]) * self.phantom_radius
        
        angles = np.array([0, 90, 180, 270]) - self.phantom_angle
        dists = np.array([0.24, 0.24, 0.24, 0.24]) * self.phantom_radius
        
        bg_rois = []
        for dist, angle in zip(dists, angles):
            roi = LowContrastDiskROI(self.image, angle, self.phantom_radius*0.03, dist, self.phantom_center,
                                     0.05)
            bg_rois.append(roi)
        avg_bg = np.mean([roi.pixel_value for roi in bg_rois])

        # create X ROIs to sample the low contrast holes
        angles = np.array([77, 116, 134.5, 0, 13, 77, 142, 153, -21, -29, -107, 182, 174, -37, -55, -105, 206.5, 189.5, -48.1, -67.8]) - self.phantom_angle-76
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

