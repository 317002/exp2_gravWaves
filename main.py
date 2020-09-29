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
import pylab
import warnings
import numpy as np
import os.path as pat
import matplotlib.pyplot as plt
import matplotlib as mpl

from collections import namedtuple
from math import sqrt

####Prodject Specific Imports####
#tools to manipulate the data
import pycbc
import pycbc.filter as filter
import pycbc.catalog as catalog

#used to gather information about some event id
import gwosc.datasets as datasets

#used to download the strain data
from gwpy.timeseries import TimeSeries
#used to generate matched filtering model vector
from pycbc.waveform import get_td_waveform
#methode to load in strain data from pycbc librabry
from pycbc.catalog import Merger
#another filtering methode
from pycbc.filter import sigma


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
    def __init__(self,event_id,detector_id,timeInterval):
        self.event_id     = event_id
        self.detector_id  = detector_id
        #gets the total time duration of the strain data
            #its needed for the mergerStrain type data
        self.eventDuration = sum([abs(c) for c in timeInterval])
        #gets the gps time stamp of the merger event
        self.gps = int(datasets.event_gps(self.event_id))
        #generating the time interval centered around the gps time
        self.timeInterval = [c + self.gps for c in timeInterval]

        self.strain = None
        #strain data handler from a different library
        self.mergerStrain = None


    def getStrainData(self):
        '''
        gathers data in the interval from a to b centered around the GPS time
        stamp for the merger event
        '''
        #sets the number of cores that the strain data loader can use
        args = {'nproc':8}
        #checking that the ids check out
        with warnings.catch_warnings():
            #ignores user warnings
            warnings.simplefilter('ignore')
            #downloading and storing in an obj the grav wave strain data
            strainData = TimeSeries.fetch_open_data(self.detector_id,\
                                                    *self.timeInterval,\
                                                    cache = True,\
                                                    verbose = True,\
                                                    **args)
        return strainData

    def getMergerStrain(self):
        '''
        Uses the pycbc library to load in the strain data for the purposes
        of using the matched filtering methode. It dosent go through the gwpy
        library that would normally use for the fft and spectrogram stuff
        '''
        merger = Merger(self.event_id)
        return merger.strain(self.detector_id)



    def mathched_filtering(self,m1,m2,f_highPass = 15,\
                                      fft_crop = 2,\
                                      psd_interval = 4,\
                                      genWave_f_lowerbound = 20,\
                                      snrCrop = 4):
        #done to avoid loading the data every time when used in a loop
        if self.mergerStrain == None:
            #this methode takes in a duration instead of a time interval
                #This automatically pulls strain data centered around the
                #gps time stamp instead of you specifing it yourself.
            self.mergerStrain = self.getMergerStrain()

        merger = Merger(self.event_id)


        '''
        There is an issue for how the strain data is read using this methode
        when being used with the filter.highpass methode

        Need to find a conversion so that a custome time interval can be used
        when selecting a data set
        '''

        #changing from the class wide strain array to a local one at the same
            #time of performing the highpass filtering.
        strain = filter.highpass(self.mergerStrain, f_highPass)
        strain = filter.resample_to_delta_t(strain, 1.0/2048)

        #removing discontinuities errors that form at the end due to resampling
        conditioned = strain.crop(fft_crop,fft_crop)
                                                     #crops off the first
                                                     #and last two seconds
        #generating the psd, thats used in the matched filtering methode
            #the psd is used to weight "the frequency components of the
            #potential signal and data by the noise amplitude"
        psd = conditioned.psd(psd_interval)
        #this matches the psd to our conditioned strain data
        psd = pycbc.psd.interpolate(psd, conditioned.delta_f)
        #this generated a 1/psd that is used to further filter the data
        psd = pycbc.psd.inverse_spectrum_truncation(psd,\
                                        psd_interval*conditioned.sample_rate,\
                                        low_frequency_cutoff=f_highPass)

        #Generating matched filtering waveform
        hp, hc = get_td_waveform(approximant="SEOBNRv4_opt",
                     mass1=m1,
                     mass2=m2,
                     delta_t=conditioned.delta_t,
                     f_lower=genWave_f_lowerbound)

        #Resizing the matched filtering wave form to the size of the our data
        hp.resize(len(conditioned))
        #shifting the moldeled wave form to the aproximant location of the
            #merger event
        template = hp.cyclic_time_shift(hp.start_time)
        #generating the signal to noise ratio data set
        snr = filter.matched_filter(template, conditioned,
                     psd=psd, low_frequency_cutoff=genWave_f_lowerbound)

        #cropping out the problamatic data points. There are discontinuitie
            #errors at the ends of the interval
        snr          = snr.crop(snrCrop + psd_interval,snrCrop)
        snrPeakIndex = abs(snr).numpy().argmax()
        snrPeak      = abs(snr)[snrPeakIndex]
        snrPeakTime  = snr.sample_times[snrPeakIndex]

        # # Shift the template to the peak time
        # dt = snrPeakTime - conditioned.start_time
        # aligned = template.cyclic_time_shift(dt)
        #
        # # scale the template so that it would have SNR 1 in this data
        # aligned /= sigma(aligned, psd=psd, low_frequency_cutoff=20.0)
        #
        # # Scale the template amplitude and phase to the peak value
        # aligned = (aligned.to_frequencyseries() * snrPeak).to_timeseries()
        # aligned.start_time = conditioned.start_time

        # We do it this way so that we can whiten both the template and the data
        white_data = (conditioned.to_frequencyseries() / psd**0.5).to_timeseries()
        white_data = white_data.highpass_fir(f_highPass*2, 512).lowpass_fir(230, 512)

        # Select the time around the merger
        white_data = white_data.time_slice(merger.time-.1, merger.time+.1)


        outputFormater = namedtuple('signal_to_noise_ratio_data',\
            ['snr','snrPeakIndex','snrPeak','snrPeakTime','white_data'])
        #returning the signal to noise ratio
        return outputFormater(snr,snrPeakIndex,snrPeak,snrPeakTime,white_data)




def gravWaveFittingFunc(x,A1,A2,w1,w2,n1,n2,a1,a2,z1,z2,center):
    base = lambda A,w,n,a,z,x:A*np.cos(w*np.abs(x)**n)*np.exp(-a*np.abs(x)**z)

    center = .001

    first = x[np.where(x <= center)]
    last  = x[np.where(x >= center)]
    # last = last[1:]

    first = base(A1,w1,n1,a1,z1,first)
    last  = base(A2,w2,n2,a2,z2,last)
    result= np.append(first,last)


    return result



'''
We are looking at GW150914 and GW170814 events on the L1 and H1 detectors
'''

def main():

    #initializing the handler for the same event with two different detectors
    l1_event = eventHandler('GW150914','L1',(-1,1))
    h1_event = eventHandler('GW170104','L1',(-1,1))


    # snr,snrPeakIndex,snrPeak,snrPeakTime,white_data = \
    #                                         l1_event.mathched_filtering(30,30)
    snr,snrPeakIndex,snrPeak,snrPeakTime,white_data = \
                                            h1_event.mathched_filtering(31.2,19.4,f_highPass=30)


    peak_shift = white_data.sample_times[white_data.numpy().argmax()]
    outInterval = white_data.sample_times - peak_shift

    # plt.plot(white_data.sample_times,white_data)

    fit_param = scipy.optimize.curve_fit(gravWaveFittingFunc,\
                                         outInterval,
                                         white_data)[0]

    print(fit_param)

    z = np.linspace(outInterval[0],\
                    outInterval[-1],\
                    100000)
    # print(z)
    f = gravWaveFittingFunc(z,*fit_param)

    plt.plot(outInterval,white_data)
    plt.plot(z,f,c = 'r',alpha = .5)


    plt.show()





if __name__ == '__main__':
    main()
