# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 15:05:35 2018

@author: cvg4
"""


import pylinac

from pylinac.core.roi import Point
import numpy as np
from skimage import feature, filters, morphology, measure
from scipy.spatial import ConvexHull

import matplotlib.pyplot as plt

class BWH_LasVegas_Definition():
    def __init__(self):
        self.phantom = (((-25.,0.), 7.5),((-25.,-20.), 7.5),((-25.,-40.), 7.5),
                        ((-10.,40.),5.),((-10.,20.),5.),((-10.,0.),5.),((-10.,-20.),5.),((-10.,-40.),5.),
                        ((5.,40.),3.5),((5.,20.),3.5),((5.,0.),3.5),((5.,-20.),3.5),((5.,-40.),3.5),
                        ((17.,40.),2.),((17.,20.),2.),((17.,0.),2.),((17.,-20.),2.),((17.,-40.),2.),
                        ((27.,40.),1.),((27.,20.),1.),((27.,0.),1.),((27.,-20.),1.),((27.,-40.),1.),
                        ((35.,40.),0.5),((35.,20.),0.5),((35.,0.),0.5),((35.,-20.),0.5),((35.,-40.),0.5))
        
        
    def convert2Polar(self,dpmm=1., frac=1.):
        dist = np.array([np.sqrt(x[0][0]**2+x[0][1]**2)*dpmm for x in self.phantom ])
        angle = np.array([np.rad2deg(np.arctan2(x[0][1],x[0][0])) for x in self.phantom])
        radius = np.array([frac*x[1]*dpmm for x in self.phantom])
        
        return dist, angle,radius
        

class BWH_LowContrastDiskROI(pylinac.roi.LowContrastDiskROI):
    def __init__(self, array, angle, roi_radius, dist_from_center, phantom_center, contrast_threshold=None,
                     background=None, background_noise=None, cnr_threshold=None):
        """
        Parameters
        ----------
        contrast_threshold : float, int
            The threshold for considering a bubble to be "seen".
        """
        super().__init__(array, angle, roi_radius, dist_from_center, phantom_center)
        self.contrast_threshold = contrast_threshold
        self.cnr_threshold = cnr_threshold
        self.background = background
        self.background_noise = background_noise

    @property
    def contrast_to_noise(self):
        """The contrast to noise ratio of the bubble: (Signal - Background)/Stdev."""
        return abs(self.pixel_value - self.background) / self.background_noise
        

class BWH_LV(pylinac.LasVegas):

    def __init__(self, filepath):

        super().__init__(filepath)
        self.image.array = np.array(self.image.array,dtype=np.float)
        self._processing()
        self._phantom_ski_region = None
        self.image.check_inversion(position=(0.1, 0.1))
        self.cnr_treshold = 3

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
        expectedAreaLV = 150.*150.
        
        # processing crops the image to get rid of the jaws. This makes the
        # detection of the phantom more reliable
        thresh = filters.threshold_otsu(self.image.array)
        bw = morphology.closing(self.image.array> thresh, morphology.square(3))
        edges = feature.canny(bw)
        

        
        labeled = measure.label(edges)
        regions = measure.regionprops(labeled, intensity_image=bw)
        
        area = [(reg.bbox_area,reg) for reg in regions]
        maxArea = float(max(area)[0])/(self.image.dpmm**2)
        
        if maxArea/expectedAreaLV >1.05:
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
        points = feature.corner_peaks(feature.corner_harris(bw),min_distance = 20)
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
        
        angle = np.array([np.rad2deg(np.arctan2(x[1]-self.phantom_center.y,x[0]-self.phantom_center.x)) for x in self.edges])
        indSort = np.argsort(angle)
        angle = np.sort(angle)
        self.edges = self.edges[indSort,:]
        self.phantom_angle  = np.mean(angle-(-135, -45, 45, 135))



        
    def _determine_low_contrast(self):
        """Sample the detail contrast regions."""
        # create 4 ROIs on each side of the phantom to determine the average background
        angles_Noise = np.array([0, 90, 180, 270]) - self.phantom_angle
        dists_Noise = np.array([60,60,60,60])*self.image.dpmm
        radius_Noise = 7.5*self.image.dpmm
        bg_rois = []
        for dist, angle in zip(dists_Noise, angles_Noise ):
            roi = BWH_LowContrastDiskROI(self.image, angle,radius_Noise, dist, self.phantom_center)
            bg_rois.append(roi)

        
        avg_bg = np.mean([roi.pixel_value for roi in bg_rois])
        std_bg =np.sqrt(np.sum([roi.std**2 for roi in bg_rois]))/np.sqrt(len(bg_rois))


        # create X ROIs to sample the low contrast holes
        phantom = BWH_LasVegas_Definition()
        dists, angles,roi_radii  = phantom.convert2Polar(self.image.dpmm,0.75)
        angles +=180-self.phantom_angle
        

        rois = []
        for dist, angle, radius in zip(dists, angles, roi_radii):
            roi = BWH_LowContrastDiskROI(self.image, angle, radius, dist, self.phantom_center,
                                     self.threshold, avg_bg,std_bg,self.cnr_treshold)
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
        visCirc = [self.lc_rois[i].contrast_to_noise>self.cnr_treshold for i in [3,8,13,18,23]];

  
        try:
            self.smallestVisCirc = min([i for i, x in enumerate(visCirc ) if not x])
        except:
            self.smallestVisCirc = 5
      

          
    def plot_analyzed_image(self, image=True, low_contrast=True, show=True):
        """Plot the analyzed image, which includes the original image with ROIs marked and low-contrast plots.

        Parameters
        ----------
        image : bool
            Show the image.
        low_contrast : bool
            Show the low contrast values plot.
        show : bool
            Whether to actually show the image when called.
        """
        num_plots = sum((image, low_contrast))
        if num_plots < 1:
            return
        fig, axes = plt.subplots(1, num_plots)
        fig.subplots_adjust(wspace=0.4)
        if num_plots < 2:
            axes = (axes,)
        axes = iter(axes)

        if image:
            img_ax = next(axes)
            # self.image.plot(ax=img_ax, show=False, vmin=self.bg_rois[0].pixel_value*0.92, vmax=self.bg_rois[0].pixel_value*1.08)
            self.image.plot(ax=img_ax, show=False)
            img_ax.plot((self.edges[0,1],self.edges[1,1]),(self.edges[0,0],self.edges[1,0]),'y')
            img_ax.plot((self.edges[1,1],self.edges[2,1]),(self.edges[1,0],self.edges[2,0]),'y')
            img_ax.plot((self.edges[2,1],self.edges[3,1]),(self.edges[2,0],self.edges[3,0]),'y')
            img_ax.plot((self.edges[3,1],self.edges[0,1]),(self.edges[3,0],self.edges[0,0]),'y')
            img_ax.axis('off')
            img_ax.set_title('Las Vegas Phantom Analysis')

            # plot the low contrast ROIs
            for roi in self.lc_rois:
                 if roi.contrast_to_noise>self.cnr_treshold:
                     color = 'b' 
                 else:
                     color = 'r'
                 roi.plot2axes(img_ax, edgecolor=color)
                 
            for roi in self.bg_rois:
                roi.plot2axes(img_ax, edgecolor='g')


        # plot the low contrast values
        if low_contrast:
            lowcon_ax = next(axes)
            self._plot_lowcontrast(lowcon_ax, self.lc_rois, self.threshold)

        if show:
            plt.show()
        
        
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

