'''
Written for Python 3.7.6

The purpose of this script it to be able to parse through gravity wave strain
data detections reported by the Gravitational Wave Open Science Center. With the
goal being to derive our own conclusions about the nature of two different merger
events. For exsample, we would like to find the mass of each object present
in the merger event, the chirp mass of the event, the peak freq of the measured
gravity waves, the freq of rotation of the two bodies throughout the merger and
the relative velocity each body has with the other.
'''

####Imports on pip Wheel####
import os
import scipy
import warnings
import numpy as np
import os.path as pat
import matplotlib.pyplot as plt

from collections import namedtuple

####Prodject Specific Imports####
#tools to manipulate the data
import pwcbc

#used to gather information about some event id
import gwosc.datasets as datasets

#used to download the strain data
from gwpy.timeseries import TimeSeries




####Constant Perameters####

#file location to the file that cotains the perameters used by the script
procPeramLoc = pat.join(os.getcwd(),'processingParameters.txt')
#text file containg the id codes for merger events and the coresponding detector
#codes the measurments were made on
mergerEventListLoc = pat.join(os.getcwd(),'listOfMergerEvents.txt')




class eventHandler:
    '''
    -have is read a perameter file here
    -going to haft ot wright a methode for reading the file
    '''
    def __init__(self,event_id,detector_ids,timeInterval):
        self.event_id = event_id
        self.detector_ids = detector_ids
        self.a,self.b = timeInterval

        #the gps time stamp of the merger event
        self.gps = int(datasets.event_gps(self.event_id))

        if self.a > self.b:
            raise Exception('timeInterval[0]:{} must be ' + \
                        'less then timeInterval[1]:{}'.format(self.a,self.b))

    def detectorID_checker(self,detector_ids):
        #making sure that detector_ids is well defined
        if type(detector_ids) != type(None):
            if type(detector_ids) == str:
                detector_ids = [detector_ids]
            #checking that the given ids can be used
            for detector_id in detector_ids:
                #check = False when none of the ids defined in init match the
                #given id
                check = max([detector_id == c for c in self.detector_ids])
                if check == False:
                    fmt = 'detector {} is not in the list of avalible ones: {}'
                    raise Exception(fmt.format(detector_id,self.detector_ids))
        else:
            #if none are given gather the ones defined in init
            detector_ids = self.detector_ids

        return detector_ids


    def getStrainData(self,detector_ids = None):
        '''gathers data in the interval from a to b centered around the GPS time
        stamp for the merger event.

        a:float
            The starting point of the time interval
        b:float
            The ending point of the itnerval
        detector_id:list(string) [None]
            The id of the detector the strain data should be pulled from. If no
            id is given, then the data for all avalible detectors will be given
        '''

        #checking that the ids check out
        detector_ids = self.detectorID_checker(detector_ids)
        #the inteval of time to select the strain data for
        interval = (self.gps + self.a,self.gps + self.b)
        #dict init
        strainData = {}
        for detect_id in detector_ids:
            #loading the strain data into a dictionary keyed to the id
            #coresponding to the detector it was measured from
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                strainData[detect_id] = \
                                        TimeSeries.fetch_open_data(detect_id,\
                                                           *interval,\
                                                           cache = True,\
                                                           verbose = True)
        return strainData

    def getQtransformData(self,fRange,qRange,detector_ids = None):
        '''This methode is the one that actully returns the qtransform
        data for the event.
        '''
        #checking that the ids check out
        detector_ids = self.detectorID_checker(detector_ids)
        #if no strain data is present
        strainData = self.getStrainData(detector_ids)
        qData = {}
        for detect_id in detector_ids:
            qData[detect_id] = strainData[detect_id].q_transform(frange=fRange,\
                                                                 qrange=qRange)

        return qData

    def qPlotter(self,qdata,title = None):
        '''Tool for plotting the q-transform data for the event
        '''
        plot = qdata.plot()
        ax = plot.gca()
        ax.set_epoch(self.gps)
        ax.set_ylim(30, 500)
        ax.set_yscale('log')
        ax.colorbar(label="Normalised energy")

        if title != None:
            ax.set_title(title)


    def getFFT(self,detector_id,frange):
        #getting the raw strain data
        strainData = self.getStrainData(detector_ids)
        #windowing the strain data
        strainData = *scipy.signal('hann',strainData.size)

        return strainData.fft(])

    def mathched_filtering(self,detector_id):
        #loading in the raw data
        strainData = self.getStrainData(detector_id)
        ####Conditioning the data####
        strain = pwcbc.filter.highpass(strain, 15.0)
        strain = pwcbc.filter.resample_to_delta_t(strain, 1.0/2048)
        #removing discontinuities errors that form at the end due to
        #resampling
        strain = strain.crop(2,2) #crops off the first and last two seconds




'''
We are looking at GW150914 and GW170814 events on the L1 and H1 detectors
'''

def main():
    #init constants
    # global procPeramLoc,\
    #        mergerEventListLoc
    #reading the list of events that should be looked at from the text file
    # events = eventListReader(mergerEventListLoc)/

    #The detectors that we wish to look at for this lab
    detectors = ['L1','H1']
    event_1 = eventHandler('GW150914',detectors)
    event_2 = eventHandler('GW170814',detectors)










if __name__ == '__main__':
    main()
