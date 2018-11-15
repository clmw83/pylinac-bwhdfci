# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib
import matplotlib.pyplot as plt
import pydicom as dcm
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import xlrd
import trimesh
import os
import copy
import time
from pylinac.core import pdf
from reportlab.lib.units import cm
import io,urllib
import base64

inch2cm = 2.54

def scanFolder(folder):
    dcmDose  = {}
    sdData   = {}
    beamData = {}
    rtplan = None
    for file in os.listdir(folder):
        if file.endswith('.dcm'):
            data = dcm.read_file(os.path.join(folder,file))
            if data.Modality=='RTDOSE':
                temp = DicomDose(data)
                dcmDose[temp.get_BeamNumber] = temp
            elif data.Modality=='RTPLAN':
                rtplan = DicomPlan(data)
                beamData = rtplan.get_BeamSeq
                for bs in data.BeamSequence:
                    beamData[bs.BeamNumber] = bs.BeamName
                    
        elif file.endswith('.xlsx'):
            wb = xlrd.open_workbook(os.path.join(folder,file)) 
            for sn in wb.sheet_names():
                sheet = wb.sheet_by_name(sn)
                temp,fieldname = processExcelSheet(sheet)
                sdData[fieldname] = temp


            

    temp = None
    for key,dose in dcmDose.items(): 
        if temp is None:
            temp = copy.deepcopy(dose)
        else:
            temp.addDose(dose)
    dcmDose['composite'] = temp 
    
    gCalc = {}
    for key,dose in dcmDose.items():
        try:
            gCalc[beamData[key]] = GammaCalc(sdData[beamData[key]],dose)
        except:
            gCalc[key] = GammaCalc(sdData[key],dose)
            
    
    return gCalc,rtplan

def generateGammaReport(gCalc,rtplan):
    gantry = rtplan.get_Gantry
    energy = rtplan.get_Energy
    columns=['Fields','Gantry','Energy','Norm Dose','Dose Dev','DTA','Gamma','Median Dose Dev']
    data = [] 
    for key,gc in gCalc.items():
        col = []
        col.append(key)
        try:
            col.append(gantry[key])
        except:
            col.append('')
        try:
            col.append(energy[key]+' MV')
        except:
            col.append('')
            
        col.append('{:1.3f} Gy'.format(gc.reference.maxDose))
            
        v,passing = gc.getDoseDevPassingRate()
        col.append('{:3.1f}%'.format(passing))
        
        v,passing = gc.getDTAPassingRate()
        col.append('{:3.1f}%'.format(passing))

        v,passing = gc.getGammaPassingRate()
        col.append('{:3.1f}%'.format(passing))
        
        data.append(col)
            
    

    plt.figure()
    plt.axis('off')
    plt.axis('tight')
    table = plt.table(cellText=data,colLabels=columns,rowLabels=None,loc='center')
    table.set_fontsize(20)
    table.scale(1,2)
    return table
    

class DicomPlan:
    
    def __init__(self,dcmData):
        self.data = dcmData
        
    @property
    def get_ID(self):
        return self.data.PatientID
    @property
    def get_Name(self):
        return str(self.data.PatientName).replace('^',', ')
    @property
    def get_BeamSeq(self):
        beamData = {}
        for bs in self.data.BeamSequence:
                    beamData[bs.BeamNumber] = bs.BeamName
        return beamData
    @property
    def get_PlanName(self):
        return self.data.RTPlanLabel
    @property
    def get_Gantry(self):
        beamData = {}
        for bs in self.data.BeamSequence:
            start = bs.ControlPointSequence[0].GantryAngle
            stop  = bs.ControlPointSequence[int(bs.NumberOfControlPoints)-1].GantryAngle
            beamData[bs.BeamName] = "%s to %s"%(start,stop)
        return beamData
    @property
    def get_Energy(self):
        beamData = {}
        for bs in self.data.BeamSequence:
            beamData[bs.BeamName] = str(bs.ControlPointSequence[0].NominalBeamEnergy)
            
        return beamData
    @property
    def get_Device(self):
        return self.data.BeamSequence[0].TreatmentMachineName
        
        
    
        
    

class DicomDose():
    
    def __init__(self,dcmData):
        
        self.data     = dcmData
        self.doseCube = np.array(self.data.pixel_array)*self.data.DoseGridScaling
        self.doseFun  = RegularGridInterpolator((self.get_Z(),self.get_Y(),self.get_X()),self.doseCube,bounds_error=False,method='linear',fill_value=None)
        self.maxDose  = np.max(self.doseCube)


        print(f'Modality:\t\t{self.data.Modality}')
        print(f'No Colums:\t\t{self.data.Columns}')
        print(f'No Rows:\t\t{self.data.Rows}')
        print(f'Number of Frames:\t{self.data.NumberOfFrames}')
        print(f'Pixel Spacing:\t\t{self.data.PixelSpacing}')
        print(f'Image Orientation:\t{self.data.ImageOrientationPatient}')
        print(f'Image Position:\t\t{self.data.ImagePositionPatient}')
        print(f'ReferencedBeamNumber: {self.data.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence[0].ReferencedBeamSequence[0].ReferencedBeamNumber}')
    
     
    def addDose(self,addDCM):

        if np.all(self.get_X()==addDCM.get_X()) and np.all(self.get_Y()==addDCM.get_Y())  and np.all(self.get_Z()==addDCM.get_Z()):
            self.doseCube += addDCM.doseCube
            self.doseFun  = RegularGridInterpolator((self.get_Z(),self.get_Y(),self.get_X()),self.doseCube,bounds_error=False,method='linear',fill_value=None)
            self.maxDose  = np.max(self.doseCube)
        else:
            raise NameError('Cube dimension missmatch')

        
        
    def get_X(self):
        return self.data.ImagePositionPatient[0]+np.arange(self.data.Rows)*self.data.PixelSpacing[0]
    
    def get_Y(self):
        return self.data.ImagePositionPatient[1]+np.arange(self.data.Rows)*self.data.PixelSpacing[1]
    
    def get_Z(self):
        return self.data.ImagePositionPatient[2]+np.array(self.data.GridFrameOffsetVector)  
    
    def get_dZ(self):
        return self.data.GridFrameOffsetVector[1]-self.data.GridFrameOffsetVector[0]
   
    @property
    def get_BeamNumber(self):
        return self.data.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence[0].ReferencedBeamSequence[0].ReferencedBeamNumber
    
    @property
    def get_MaxDose(self):
        return self.maxDose
    
    def get_Gradient(self,x,y,z,epsilon=0.01):
        dose0 = self.interpolateDose(x,y,z)
        return [(self.interpolateDose(x+epsilon,y,z)-dose0)/epsilon,(self.interpolateDose(x,y+epsilon,z)-dose0)/epsilon,(self.interpolateDose(x,y,z+epsilon)-dose0)/epsilon]/self.get_MaxDose

        
    
        
    def interpolateDose(self,x,y,z):
        try: 
            return np.array([self.doseFun((zz,yy,xx)) for xx,yy,zz in zip(x,y,z) ]).squeeze()
        except:
            return np.float64(self.doseFun((z,y,x)))
        
    def get_IdX(self,x):
        return int((x-self.data.ImagePositionPatient[0])/self.data.PixelSpacing[0])
    
    def get_IdY(self,y):
        return int((y-self.data.ImagePositionPatient[1])/self.data.PixelSpacing[1])
    
    def get_IdZ(self,z):
        return np.argmin(np.abs((z-self.data.ImagePositionPatient[2])-np.array(self.data.GridFrameOffsetVector)))

class Diode:
    def __init__(self,x,y,z,r,value):
        self.x = 10.*x
        self.y = -10.*y  # definition by scandidose is not same as in eclipse
        self.z = 10.*z
        self.r = 10.*r
        self.value = value/100.
        
 # gamma contains dictionary [doseDev,SpatialDev]      
class GammaDiode(Diode):
    def __init__(self,x,y,z,r,value):
        
        super().__init__(x,y,z,r,value)
        self.doseDev = None
        self.gradient = None
        self.dta = None
        self.radiusDoseMat = None
        self.gamma = None
        
        icosphere = trimesh.creation.icosphere(0,1)
        self.directionsDTA = icosphere.vertices
        

    def calcAll(self,doseCube,directionsGamma,dist=[20.,9],doseRange = [0,5],gradThreshold = 0.01):
        #        
        Dcalc = doseCube.interpolateDose(self.x,self.y,self.z)
        if doseRange[0]<=Dcalc and Dcalc<=doseRange[1]:
            self.doseDev = (self.value-Dcalc)/doseCube.get_MaxDose
        else:
            self.doseDev =None  
        
#        
        radius = np.linspace(0,dist[0],dist[1]) 

        start = time.time()
        
        if self.gradient>gradThreshold:
            for dir in directionsGamma:
                x = dir[0]*radius
                y = dir[1]*radius
                z = dir[2]*radius
                Dcalc = doseCube.interpolateDose(x,y,z)
            
                deltaD = self.value-Dcalc
                print(np.sign(deltaD))
        
        
        end = time.time()
        print(end - start)
      
    def setGradient(self,gradient):
        self.gradient = gradient
        
    def getDoseDev(self,doseCube,doseRange = [0,5],norm='max'):
        Dcalc = doseCube.interpolateDose(self.x,self.y,self.z)
        if doseRange[0]<=Dcalc and Dcalc<=doseRange[1]:
            DMeas = self.value
            if norm=='max':
                self.doseDev = (DMeas-Dcalc)/doseCube.get_MaxDose
        else:
            self.doseDev =None
    
    def getDTA(self,doseCube,dtaRadius=[20.,9],gradThreshold = 0.01,doseThreshold=0.005):
        
        radius = np.linspace(0,dtaRadius[0],dtaRadius[1])
        dRadius = radius[1]-radius[0]
              

        if self.gradient>gradThreshold:
            self.dta = dtaRadius[0]
            if self.doseDev is not None:
                if self.doseDev<doseThreshold:
                    self.dta = 0
                        
            for p in self.directionsDTA:
                oldR = None
                oldD = None
                for r in radius:
                    
                    if r>self.dta:
                        break
                    
                    dir = p*r
                    tempR = dir+[self.x,self.y,self.z] 
                    tempD = doseCube.interpolateDose(tempR[0],tempR[1],tempR[2])

                    if oldD is not None:
                        if np.sign(oldD-self.value)!=np.sign(tempD-self.value):
                            if tempD!=oldD and tempD!=self.value:
                                alpha = -(oldD-self.value)/(tempD-oldD)
                                temp = dRadius*alpha+oldR
                                if temp <self.dta:
                                    self.dta = temp
                                break
            
                    oldD = tempD
                    oldR = r

    @staticmethod
    def minDistLineOrigin(x0,y0,x1,y1): 
        m = (y1-y0)/(x1-x0)
        b = y0-m*x0
        xx = (-m*b)/(m**2+1.)
        
        if xx<x0:
            dist = np.sqrt(x0**2+y0**2)
        elif x1<xx:
            dist = np.sqrt(x1**2+y1**2)
        else:
            dist = np.abs(b)/(np.sqrt(m**2+1.))
        
        return dist
        
    def calcGamma(self,doseCube,directionsGamma,rad=[10.,5],includeDTA=[0.1,5.0],doseDevRange=None,spatialDevRange=None):
        # dta and gamma
        radius = np.linspace(0.,rad[0],rad[1])

        if includeDTA[0]*doseCube.get_MaxDose<=self.value and self.value<includeDTA[1]*doseCube.get_MaxDose:
            self.radiusDoseMat = []
            
            self.gamma = {}
            for ddr in doseDevRange:
                for sdr in spatialDevRange:     
                    self.gamma[ddr,sdr] = np.inf
            
            
            for p in directionsGamma:
                oldR = None
                oldD = None
                
                for r in radius:
                    
                    dir = p*r
                    tempR = dir+[self.x,self.y,self.z] 
                    tempD = (self.value-doseCube.interpolateDose(tempR[0],tempR[1],tempR[2]))/doseCube.get_MaxDose
                    
                    if oldD is not None:
                        for ddr in doseDevRange:
                            for sdr in spatialDevRange:
                                tempGamma = self.minDistLineOrigin(oldR/sdr,oldD/ddr,r/sdr,tempD/ddr)
                                if tempGamma<self.gamma[ddr,sdr]:
                                    self.gamma[ddr,sdr] = tempGamma
                                else:
                                    break
                    oldR = r
                    oldD = tempD
                
                                    
                                
                                 

    @property 
    def validDoseDev(self):
        return (self.doseDev is not None)

    def passedDoseDev(self,doseThreshold):
        return np.abs(self.doseDev)<=doseThreshold

    def plotColorDoseDev(self,doseThreshold):
        return 'green' if self.passedDoseDev(doseThreshold) else 'red'
        
    
    def passedDTA(self,spatialThreshold):
        return np.abs(self.dta)<=spatialThreshold

    def plotColorDTA(self,doseThreshold):
        return 'green' if self.passedDTA(doseThreshold) else 'red'   
    
    def passedGamma(self,doseDev,spatialDev):
        return np.abs(self.gamma[doseDev,spatialDev])<=1.0
    
    def plotColorGamma(self,doseDev,spatialDev):
        return 'green' if self.passedGamma(doseDev,spatialDev) else 'red' 
    
        
class DetectorBoard():
    def __init__(self):
        self.diodes = []
        self.X = None
        self.Y = None
        self.Z = None
        self.R = None
        self.rows = None
        self.columns = None
        self.angle = None
        
    def addDiode(self,diode):
        self.diodes.append(diode)
        
    def ready(self):
        self.getDimension()
        self.angle = np.round(np.rad2deg(np.mean([np.arctan2(d.y,d.x) for d in self.diodes if (np.sign(d.r)==1)])))
    
    def getDetectorPlane(self,res=2.5):
        
        NR = int(np.round((max(self.R)-min(self.R)+2*res)/res))+1
        r = np.linspace(min(self.R)-res,max(self.R)+res,NR)
    
        NZ = int(np.round((max(self.Z)-min(self.Z)+2*res)/res))+1
        z = np.linspace(min(self.Z)-res,max(self.Z)+res,NZ)
        
        alpha = np.deg2rad(self.angle)
        points = np.array([np.round([np.cos(alpha)*rr, np.sin(alpha)*rr,zz],2) for zz in z for rr in r])
        
        return points, z, r
    
    def plotDiodes2D(self,color='r'):
        for d in self.diodes:
            self.doseThreshold
            plt.plot(d.r,d.z,'o',markersize = 2.5,color = color)  
        
        
        
                
    
    def getDimension(self):
        self.X = []
        self.Y = []
        self.Z = []
        self.R = []
        for d in self.diodes:
            self.X.append(d.x)
            self.Y.append(d.y)
            self.Z.append(d.z)
            self.R.append(d.r)
        self.X = sorted(list(set(self.X)))
        self.Y = sorted(list(set(self.Y)))
        self.Z = sorted(list(set(self.Z)))
        self.R = sorted(list(set(self.R)))
        self.rows    = len(self.Z)
        self.columns = len(self.R)
            
    
    def count(self):
        return len(self.diodes)

class GammaSettings:
    def __init__(self):
        self.doseDev = 0.03
        self.doseDevRange = [0.,5.]
        self.detGradientThreshold = 0.01
 
class GammaCalc(GammaSettings):
    
    def __init__(self,phantom,reference):
        super().__init__()
        self.phantom = phantom
        self.reference = reference
        
        self.doseDeviation()
        self.calcGradient()
        
        self.doseDevRange = np.round(np.linspace(0.01,0.1,10),3)
        self.spatialDevRange = np.round(np.linspace(1.0,10.0,10),1)
        
    @property
    def get_NormDose(self):
        self.reference.maxDose
    

    
        
    def getProfile(self,type,pIndex,bIndex,doseDev=3.,spatialDev=3.,subplot=False):
        
        if not(subplot): 
            plt.figure()
            plt.ylabel('Dose [% of max Dose]')
        else:
            plt.ylabel('Dose [%]')
            
        plt.title('Board {:d}'.format(bIndex))
        plt.grid(True)
        
        
        x = []
        y = []
        z = []
        r = []
        for d in self.phantom[bIndex].diodes:
            if ((d.r==self.phantom[bIndex].R[pIndex]) and (type == 'R')) or ((d.z==self.phantom[bIndex].Z[pIndex]) and (type == 'Z')):
                x.append(d.x)
                y.append(d.y)
                z.append(d.z)
                r.append(d.r)
                if d.gamma is None:
                    if (type == 'R'):
                        plt.plot(d.z,100*d.value/self.reference.get_MaxDose,'xb')
                    elif (type == 'Z'):
                        plt.plot(d.r,100*d.value/self.reference.get_MaxDose,'xb')
                else:
                    color = d.plotColorGamma(doseDev,spatialDev)
                    if (type == 'R'):
                        plt.plot(d.z,100*d.value/self.reference.get_MaxDose,'o',color = color)
                    elif (type == 'Z'):
                        plt.plot(d.r,100*d.value/self.reference.get_MaxDose,'o',color = color)

        if (type == 'R'):
             plt.xlabel('z [mm]')
             plt.title('R = {:2.1f}mm @ Board {:d}'.format(self.phantom[bIndex].R[pIndex],bIndex))
        elif (type == 'Z'):
            plt.xlabel('r [mm]') 
            plt.title('Z = {:2.1f}mm @ Board {:d}'.format(self.phantom[bIndex].Z[pIndex],bIndex))
        if self.reference is not None:
            refDose = 100*self.reference.interpolateDose(np.array(x),np.array(y),np.array(z))/self.reference.get_MaxDose
            if (type == 'R'):
                plt.plot(z,refDose,'-k')
         
            elif (type == 'Z'):
                plt.plot(r,refDose,'-k')
        
            

        
    def doseDeviation(self,norm='max'):
        for key,det in self.phantom.items():
            for d in det.diodes:
                d.getDoseDev(self.reference)
    
    def distance2Agreement(self):
        for key,det in self.phantom.items():
            for d in det.diodes:
                d.getDTA(self.reference)
                
    def gamma(self):
        directionsGamma=trimesh.creation.icosphere(1,1).vertices
        
 
                
        for key,det in self.phantom.items():
            for d in det.diodes:
                d.calcGamma(self.reference,directionsGamma,[10.,5],[0.2,5.0],self.doseDevRange,self.spatialDevRange)
    
    def getGamma(self,frac,deviation):
        
        gammaVec = []
        for key,det in self.phantom.items():
            for d in det.diodes:
                if d.gamma is not None:
                    gammaVec.append(d.gamma[frac,deviation])
        return gammaVec         

    def calcGradient(self,threshold=0.01):
        for key,det in self.phantom.items():
            for d in det.diodes:
                grad = self.reference.get_Gradient(d.x,d.y,d.z)
                grad = np.sqrt(grad[0]**2+grad[1]**2+grad[2]**2)
                d.setGradient(grad)

    
    def getGammaPassingRate(self,doseDev=0.03,spatialDev=3.):
        v = []
        passed = 0
        for key,det in self.phantom.items():
            for d in det.diodes:
                if d.gamma is not None:
                    val = d.gamma[doseDev,spatialDev]
                    v.append(val)
                    if val<=1:
                        passed += 1
        try:  
            passed = float(100.*passed/len(v))
        except:
              passed = 0. 
        return v,passed
        
    
    def plotGamma(self,doseDev=0.03,spatialDev=3.,range = [0,2],subplot=False):
        
        v,passed = self.getGammaPassingRate(doseDev,spatialDev)
        if not(subplot): 
            plt.figure()
        plt.grid(True)
        plt.hist(np.array(v),bins=32, range=range,weights=np.ones_like(v) * 100. / len(v),edgecolor='black')
        plt.xlabel('Gamma Index')
        plt.ylabel('Frequency [% of {:d}]'.format(len(v)))
        plt.title('{:3.1f}% passed with index<=1'.format(passed))
        if len(v)>0:
            plt.text(plt.xlim()[1]*0.5,plt.ylim()[1]*0.9,'mean {:3.2f}'.format(np.mean(v)))
            plt.text(plt.xlim()[1]*0.5,plt.ylim()[1]*0.8,'max {:3.2f}'.format(np.max(v)))


    def getDoseDevPassingRate(self,doseDev = 0.03):
        v = []
        passed = 0
        for key,det in self.phantom.items():
            for d in det.diodes:
                if d.doseDev is not None:
                    val = d.doseDev
                    v.append(100.*val)
                    if np.abs(val)<=doseDev:
                        passed +=1
        passed = float(100.*passed/len(v))  
        return v,passed
        
    
    def plotDoseDev(self,doseDev = 0.03,range=[-10,10],subplot=False):

        v,passed = self.getDoseDevPassingRate(doseDev)
        if not(subplot):                
            plt.figure()
        plt.grid(True)
        plt.hist(np.array(v),bins=32, range=range,weights=np.ones_like(v) * 100. / len(v),edgecolor='black')
        plt.xlabel('Dose Deviation [%]')
        plt.ylabel('Frequency [% of {:d}]'.format(len(v)))
        plt.title('{:3.1f}% within {:3.1f}% deviation'.format(passed,100.*doseDev)) 
        plt.text(plt.xlim()[0]+0.05*(plt.xlim()[1]-plt.xlim()[0]),0.9*plt.ylim()[1],'median {:3.2f}%'.format(np.median(v)))

    def getDTAPassingRate(self,spatialDev = 3.):
        v = []
        passed = 0
        for key,det in self.phantom.items():
            for d in det.diodes:
                if d.dta is not None:
                    val = d.dta
                    v.append(val)
                    if np.abs(val)<=spatialDev:
                        passed +=1
        passed = float(100.*passed/len(v))  
        return v,passed        
  
    def plotDTA(self,spatialDev = 3.,range=[0,10],subplot=False):

        v,passed = self.getDTAPassingRate(spatialDev)
        if not(subplot):
            plt.figure()   
        plt.grid(True)            
        plt.hist(np.array(v),bins=32, range=range,weights=np.ones_like(v) * 100. / len(v),edgecolor='black')
        plt.xlabel('Distance to Agreement [mm]')
        plt.ylabel('Frequency [% of {:d}]'.format(len(v)))
        plt.title('{:3.1f}% with DTA<={:3.1f}mm'.format(passed,spatialDev))

    def analyzeAll(self):
        self.doseDeviation()
        self.distance2Agreement()
        self.gamma()

    def plotAllBoards(self,type='gamma',threshold=[0.03,3.]):
        plt.figure()
        plt.subplot(1, 2, 1)
        self.plotDetectorBoard(0,type='gamma',threshold=[0.03,3.],subplot=True)
        plt.subplot(1, 2, 2)
        self.plotDetectorBoard(1,type='gamma',threshold=[0.03,3.],subplot=True)
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        
        
    def plotAllInOne(self):
        plt.figure()
        plt.subplot(1, 3, 1)
        self.plotDoseDev(doseDev = 0.03,range=[-10,10],subplot=True)
        plt.subplot(1, 3, 2)
        self.plotDTA(spatialDev = 3.,range=[0,10],subplot=True)
        plt.subplot(1, 3, 3)
        self.plotGamma(doseDev=0.03,spatialDev=3.,range = [0,2],subplot=True)
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    
    def plotAllProfiles(self,board=0,N=3,profile='R'):
        
        if profile=='R':
            M = np.ceil(self.phantom[board].columns/N)
            plt.figure()
        
            for i in range(self.phantom[0].columns):
                plt.subplot(M, N, i+1)
                self.getProfile('R',i,board,0.03,3.,subplot=True)   
        elif profile=='Z':
            M = np.ceil(self.phantom[board].columns/N)
            plt.figure()
        
            for i in range(self.phantom[0].rows):
                plt.subplot(M, N, i+1)
                self.getProfile('Z',i,board,0.03,3.,subplot=True)    
           
    def plotDetectorBoard(self,bIndex,type='doseDev',threshold=3.,subplot=False):
        points,z,r = self.phantom[bIndex].getDetectorPlane() 
        
        im = self.reference.interpolateDose(points[:,0],points[:,1],points[:,2]).reshape(len(z),len(r))
        extent=[r[0],r[-1],z[-1],z[0]]  #scalars (left, right, bottom, top)
        
        if not(subplot):
            plt.figure()  
            
        plt.imshow(im,extent=extent,cmap='bone')
        
        for d in self.phantom[bIndex].diodes:
            if d.validDoseDev and type is 'doseDev':
                color = d.plotColorDoseDev(threshold/100.)
                plt.plot(d.r,d.z,'o',markersize = 2.5,color = color)
            elif d.dta and type is 'dta':
                color = d.plotColorDTA(threshold)
                plt.plot(d.r,d.z,'o',markersize = 2.5,color = color)
            elif d.gamma and type is 'gamma':
                color = d.plotColorGamma(threshold[0],threshold[1])
                plt.plot(d.r,d.z,'o',markersize = 2.5,color = color)
            else:
                plt.plot(d.r,d.z,'o',markersize = 1.0,color = 'b')
            
        plt.xlabel('r [mm]')
        plt.ylabel('z [mm]')
        plt.title("Board: {:d}".format(bIndex))

  
    def generateGammaTable(self,passIf = 95.):
        
        
        columns = []
        for ddr in self.doseDevRange:
            columns.append("{:3.1f} %".format(100*ddr))
            
        rows =[]
        for sdr in self.spatialDevRange:
            rows.append("{:3.1f} mm".format(sdr))

        data = []
        color = []

        for sdr in self.spatialDevRange:
            rowData= []
            rowColor = []
            for ddr in self.doseDevRange:
                v,passed = self.getGammaPassingRate(ddr,sdr)
                rowData.append("{:3.1f}%".format(passed))
                
                if passed>passIf:
                    rowColor.append('g')
                else:
                    rowColor.append('r')
            
            color.append(rowColor)   
            data.append(rowData)

        plt.figure()
        plt.axis('off')
        plt.axis('tight')
        return plt.table(cellText=data,rowLabels=rows,colLabels=columns,loc='center',cellColours=color)
        
        
def processExcelSheet(sheet):
    fieldName = 'composite'
    for ir in range(sheet.nrows):
        if (sheet.cell_value(ir,0)=='Distance'):
            indDist = ir
        if (sheet.cell_value(ir,0)=='X (dicom-left)'):
            indX = ir
        if (sheet.cell_value(ir,0)=='Y (dicom-down)'):
            indY = ir
        if (sheet.cell_value(ir,0)=='Z (dicom-head)'):
            indZ = ir
        if (sheet.cell_type(ir,0)==1):
            if 'Beam' in sheet.cell_value(ir,0):
                fieldName = (sheet.cell_value(ir,0)[6:])
            
                
    detector = {0:DetectorBoard(),1:DetectorBoard()}
    for ir in range(indZ+1,sheet.nrows):
        iBoard = 0
        for ic in range(1,sheet.ncols):
            if sheet.cell_type(ir,ic)==2:
                x = sheet.cell_value(indX,ic)
                y = sheet.cell_value(indY,ic)
                z = sheet.cell_value(ir,0)
                r = sheet.cell_value(indDist,ic)
                value = sheet.cell_value(ir,ic)
                detector[iBoard].addDiode(GammaDiode(x,y,z,r,value))
            
            if sheet.cell_type(ir,ic)==1:
                detector[iBoard].ready()
                iBoard = 1
    detector[iBoard].ready()      
    return detector,fieldName

def publish_pdf(gamma,rtplan,path=None, filename=None, author=None, unit=None, notes=None, open_file=False):


    filename = os.path.join(path,filename)
    # header (name and logo) 
    canvas = pdf.create_pylinac_page_template(filename, analysis_title='IMRT QA Report',author=author, unit=unit, file_name=None,file_created=None)
    # header (patient info)
    text = ['{:s}'.format(rtplan.get_Name),'{:s}'.format(rtplan.get_ID),'{:s}'.format(rtplan.get_PlanName)]
    pdf.draw_text(canvas, x=15* cm, y=28.25 * cm, text=text,fontsize=13)    

    # start content    
    text = ['Summary']
    actHeight = 25.5*cm
    pdf.draw_text(canvas, x=1*cm, y=actHeight, text=text,fontsize=13)   
    actHeight = actHeight-2
    # summary        
    table = generateGammaReport(gamma,rtplan)
    plt.gcf().canvas.draw()
    points = table.get_window_extent(plt.gcf()._cachedRenderer).get_points()
    # add 10 pixel spacing
    points[0,:] -= 10; points[1,:] += 10
    #
    nbbox = matplotlib.transforms.Bbox.from_extents(points/plt.gcf().dpi)    
    
    bx = (nbbox.x1-nbbox.x0)*inch2cm
    by  =(nbbox.y1-nbbox.y0)*inch2cm
    m = 1/(by/(len(gamma)+1))  
    buf = io.BytesIO()
    matplotlib.pyplot.savefig(buf, format='png', bbox_inches=nbbox,dpi=900)
    matplotlib.pyplot.close('all')
    
    actHeight -= cm*(m*by)
  
    img = pdf.create_stream_image(buf)
    canvas.drawImage(img ,cm*(21-m*bx)/2-0.5, actHeight, width=m*bx* cm, height=m*by*cm, preserveAspectRatio=True)    
    actHeight -=1*cm
    #
    for key,gc in gamma.items():
        # Field Info
        actHeightL = actHeight
        actHeight -=1*cm
        text = ['Results for {:s}'.format(key)]
        actHeight0 = actHeight
        actHeight -=1*cm
        
        #  All Plots     
        gc.plotAllInOne()
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(10,3)
        fig.tight_layout()
        buf = io.BytesIO()
        matplotlib.pyplot.savefig(buf, format='png')
        matplotlib.pyplot.close('all')
        buf.seek(0)
        img1 = pdf.create_stream_image(buf)
        m = 19/(10*inch2cm)
        sizeImg1 = [19,3*inch2cm*m]
  
        actHeight -= sizeImg1[1]*cm
        actHeight1 = actHeight
        actHeight -=1*cm
        
        # Gamma Table
        table = gc.generateGammaTable()
        plt.gcf().canvas.draw()
        points = table.get_window_extent(plt.gcf()._cachedRenderer).get_points()
        # add 10 pixel spacing
        points[0,:] -= 10; points[1,:] += 10
        #
        nbbox = matplotlib.transforms.Bbox.from_extents(points/plt.gcf().dpi)    
    
        bx = (nbbox.x1-nbbox.x0)*inch2cm
        by  =(nbbox.y1-nbbox.y0)*inch2cm
        m = 18/bx  
        sizeImg2 = [m*bx,m*by]
        buf = io.BytesIO()
        matplotlib.pyplot.savefig(buf, format='png', bbox_inches=nbbox,dpi=900)
        matplotlib.pyplot.close('all')
        img2 = pdf.create_stream_image(buf)
    
        actHeight -= cm*(m*by)
        actHeight2 = actHeight
        
        if actHeight2<0:
            dh = 25.5*cm-actHeight0
            actHeight0 +=dh
            actHeight1 +=dh
            actHeight2 +=dh
            actHeight = actHeight2
            canvas.showPage()
            pdf.add_pylinac_page_template(canvas, analysis_title='IMRT QA Report')
            textHeader = ['{:s}'.format(rtplan.get_Name),'{:s}'.format(rtplan.get_ID),'{:s}'.format(rtplan.get_PlanName)]
            pdf.draw_text(canvas, x=15* cm, y=28.25 * cm, text=textHeader ,fontsize=13)  

            
        else:
            canvas.line(1 * cm, actHeightL, 20 * cm, actHeightL)   
             
        pdf.draw_text(canvas, x=1*cm, y=actHeight0, text=text,fontsize=13) 
        canvas.drawImage(img1 ,1*cm, actHeight1, width=sizeImg1[0]*cm, height=sizeImg1[1]*cm, preserveAspectRatio=True) 
        canvas.drawImage(img2 ,cm*(21-sizeImg2[0])/2, actHeight2, width=sizeImg2[0]*cm, height=sizeImg2[1]*cm, preserveAspectRatio=True) 
        actHeight -=1*cm

    pdf.finish(canvas, open_file=open_file, filename=filename)
    
    





