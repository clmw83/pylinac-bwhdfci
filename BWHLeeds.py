import pylinac

import copy
from functools import lru_cache
import io
import os.path as osp

import matplotlib.pyplot as plt
import numpy as np
from reportlab.lib.units import cm
from scipy.interpolate.interpolate import interp1d
from skimage import feature, measure

from pylinac.core import image
from pylinac.core.geometry import Point
from pylinac.core.io import get_url, retrieve_demo_file
from pylinac.core.profile import CollapsedCircleProfile
from pylinac.core.roi import LowContrastDiskROI, HighContrastDiskROI, DiskROI, bbox_center,RectangleROI
from pylinac.core import pdf

class BWHLeeds(pylinac.LeedsTOR):
    def analyze(self, low_contrast_threshold=0.005, hi_contrast_threshold=0.4, invert=False,
            angle_offset=0, hc_angle_offset=0, hc_shift=Point(0,0)):
        """Analyze the image.

        Parameters
        ----------
        low_contrast_threshold : float
            The threshold for the low-contrast bubbles to be "seen".
        hi_contrast_threshold : float
            The threshold percentage that the relative MTF must be above to be "seen". Must be between 0 and 1.
        invert : bool
            Whether to force an inversion of the image. Pylinac tries to infer the correct inversion but uneven
            backgrounds can cause this analysis to fail. If the contrasts/MTF ROIs appear correctly located but the
            plots are wonky, try setting this to True.
        angle_offset : int, float
            Some LeedsTOR phantoms have the low contrast regions slightly offset from phantom to phantom. 
            This parameter lets the user correct for any consistent angle effects of the phantom. The offset 
            is in degrees and moves counter-clockwise. Use this if the low contrast ROIs are offset from the real 
            ROIs.
        """
        self.image.check_inversion(box_size=30, position=(0.1, 0.25))
        if invert:
            self.image.invert()
        self.low_contrast_threshold = low_contrast_threshold
        self.hi_contrast_threshold = hi_contrast_threshold

        if not self._is_clockwise():
            self._flip_image_data()
        self.lc_rois, self.lc_ref_rois = self._low_contrast(angle_offset)
        self.hc_rois, self.hc_ref_rois = self._high_contrast(hc_angle_offset,hc_shift)
        
    def _high_contrast(self,angle_offset,shift):
        """Perform high-contrast analysis. This samples disks within the line-pair region and calculates
        relative MTF from the min and max values.

        Returns
        -------
        contrast ROIs : list
            :class:`~pylinac.core.roi.HighContrastDiskROI` instances of the line pairs.
        reference ROIs : list
            :class:`~pylinac.core.roi.HighContrastDiskROI` instances of the solid ROIs that
            determine the normalization value for MTF.
        """
        angle = np.degrees(self.phantom_angle)

        # sample ROIs of the reference areas
        ref_angles = [303, 271]
        ref_dists = [0.3 * self.phantom_radius, 0.25 * self.phantom_radius]
        ref_radius = 0.04 * self.phantom_radius
        rrois = []
        for nominal_angle, dist in zip(ref_angles, ref_dists):
            roi = HighContrastDiskROI(self.image, angle - nominal_angle, ref_radius, dist, self.phantom_center-shift,
                                      self.hi_contrast_threshold)
            rrois.append(roi)
        mtf_norm_val = (rrois[0].pixel_value - rrois[1].pixel_value) / (rrois[0].pixel_value + rrois[1].pixel_value)

        # sample ROIs of each line pair region
        # ordering goes from the "biggest" line pair region downward
        contrast_angles = [-144.8, -115.1, -62.5, -169.7, -153.4, -25, 169.7, 151.6, 27]
        contrast_dists = np.array([0.3, 0.187, 0.187, 0.252, 0.092, 0.094, 0.252, 0.094, 0.0958]) * self.phantom_radius
        contrast_radii = np.array([0.04, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.018, 0.018, 0.015, 0.015, 0.012]) * self.phantom_radius
        crois = []
        for nominal_angle, dist, cradius in zip(contrast_angles, contrast_dists, contrast_radii):
            roi = HighContrastDiskROI(self.image, angle + nominal_angle + 90+angle_offset, cradius, dist, self.phantom_center-shift, self.hi_contrast_threshold, mtf_norm=mtf_norm_val)
            crois.append(roi)

        return crois, rrois
    
    def uniformity(self):
        top=DiskROI(self.image,270,self.phantom_radius*0.07,self.phantom_radius*1.08,self.phantom_center)
        bottom=DiskROI(self.image,90,self.phantom_radius*0.07,self.phantom_radius*1.08,self.phantom_center)
        left=DiskROI(self.image,0,self.phantom_radius*0.07,self.phantom_radius*1.08,self.phantom_center)
        right=DiskROI(self.image,180,self.phantom_radius*0.07,self.phantom_radius*1.08,self.phantom_center)

        
        #top=RectangleROI(self.image,self.phantom_radius*0.1,self.phantom_radius*0.1,270,self.phantom_radius*1.07, self.phantom_center)        
        #bottom=RectangleROI(self.image,self.phantom_radius*0.1,self.phantom_radius*0.1,90,self.phantom_radius*1.07, self.phantom_center)        
        #right=RectangleROI(self.image,self.phantom_radius*0.1,self.phantom_radius*0.1,0,self.phantom_radius*1.07, self.phantom_center)        
        #left=RectangleROI(self.image,self.phantom_radius*0.1,self.phantom_radius*0.1,180,self.phantom_radius*1.07, self.phantom_center)        
        
        for r in (top,bottom,left,right):
            print(r.pixel_value)
            print(r.std)
        
    @property
    def phantom_angle(self):
        """Determine the angle of the phantom.

        This is done by searching for square-like boxes of the canny image. There must be two: one lead and
        one copper. The angle from the lead to the copper ROI centers is determined (modified from original,
        which used the phantom center, and was not as stable)
        Returns
        -------
        angle : float
            The angle in radians.
        """
        if self._phantom_angle is not None:
            return self._phantom_angle
        expected_length = self.phantom_radius * 0.52
        square_rois = [roi for roi in self._blobs if np.isclose(self._regions[roi].major_axis_length, expected_length, rtol=0.2)]
        if len(square_rois) != 2:
            raise ValueError("Could not find the angle of the image.")
        regions = self._regions
        idx_sort= np.argsort([regions[roi].mean_intensity for roi in square_rois])
        lead_idx=idx_sort[1]
        cu_idx = idx_sort[0]
        lead_roi = regions[square_rois[lead_idx]]
        lead_center = bbox_center(lead_roi)

        cu_roi = regions[square_rois[cu_idx]]
        cu_center = bbox_center(cu_roi)

        adjacent = lead_center.x - cu_center.x
        opposite = lead_center.y - cu_center.y
        angle = np.arctan2(opposite, adjacent)
        return angle
    
    def _mtf(self, x=50, lpm=False):
        #norm = max(roi.mtf for roi in self.hc_rois)
        lpms=[.5, .56, .63, .71, .8, .9,1,1.12,1.25]
        norm = 1.0
        ys = [roi.mtf / norm for roi in self.hc_rois]
        if lpm:
            xs=lpms
        else:
            xs = np.arange(len(ys))
        f = interp1d(ys, xs)
        try:
            mtf = f(x / 100)
        except ValueError:
            mtf = min(ys)
        return float(mtf)
