# -*- coding: utf-8 -*-
"""
Created on Wed May 24 09:55:59 2017

@author: BPPTK_COMP03
"""

from obspy.core import read, Trace, Stream, UTCDateTime
from obspy.signal.invsim import paz_to_freq_resp
from os import listdir, path, system
import numpy as np

# Seismic Data file path parameters
stations      = ['grab','jrob','klab','pasb','ijob','lab']
components    = ['-n','-e']
yearFolder   = '2016\\'
monthFolder  = '01\\'
msdExt = '\\*.msd'
gcfExt = '\\*.gcf'

# Instrument Calibration Parameter
poles = [-0.01178+0.01178j, -0.01178-0.01178j, -160+0j, -80+0j, -180+0j]
zeros = [0.0+0j, 0.0+0j]

# Find frequency response of the seismometer
h, f = paz_to_freq_resp(poles,zeros,0.4,0.005,16384,freq = True)

# Poles, Zeros, Gain, Sensitivity --> from Guralp Calibration Sheet
paz_cmg40T = {'poles': poles,
              'zeros': zeros,
              'gain': 2304000.0,
              'sensitivity': 2*2977.0*1e6}

def Gcf2msd(thePath):  
    print('Convert gcf to msd in '+ thePath)
    system("gcf2msd.exe \"" + thePath + "\"")
    
    #print('remove gcf in' + thePath)
    #system("del \"" + thePath)

def BatchMsdConverting():
    for station in stations:
        print ('CONVERTING --------> ' + station)
        for component in components:
            print ('Component ' + component)
            for i in range(30):
                if(i<9):    
                    filePath = yearFolder + monthFolder + '0' + str(i+1) + '\\' + station + component + gcfExt # Define filepath
                else:
                    filePath = yearFolder + monthFolder + str(i+1) + '\\' + station + component + gcfExt # Define filepath
                #filePath = yearFolder + monthFolder + '01\\' + station + component # Define filepath
                Gcf2msd(filePath) # Begin Converting
        
def ReadMonthy(station):
    resultStream = Stream()
    for component in components:
        stream = Stream()
        for i in range(30):
            if(i<9):    
                filePath = yearFolder + monthFolder + '0' + str(i+1) + '\\' + station + component + msdExt # Define filepath
            else:
                filePath = yearFolder + monthFolder + str(i+1) + '\\' + station + component + msdExt # Define filepath
            
            if(path.exists(yearFolder + monthFolder + str(i+1) + '\\' + station + component)):
                print('Reading... ' + filePath)
                newStream = read(filePath)
                newStream.resample(5)
                stream += newStream
            else:
                print(filePath + 'Not Found')
        
        print('MERGING--------' + component)
        stream.merge(fill_value = 'interpolate') # merge seismic data
        resultStream += stream
    return resultStream

def GetTilt(stream):
    # Clone Stream
    cstream = stream.copy()
    
    # Perform instrument correction
    cstream.simulate(paz_remove=paz_cmg40T,remove_sensitivity = True, taper = True, zero_mean=True,simulate_sensitivity=True)
    
    # Perform filering
    cstream.filter('highpass', freq = 0.00000005, corners=2, zerophase=True)
    cstream.filter('lowpass', freq = 0.005, corners=2, zerophase=True)
    
    # Perform diferentiation
    cstream.differentiate(method = 'gradient')
    
    tilt = Stream()
    for traceComp in cstream:
        # Calculate Tilt Value
        tiltVal = traceComp.data/-9.81
        trace = Trace(data=tiltVal,header= traceComp.stats)
        tilt += Stream(trace)
    
    fiveDaysInSec = 432000
    startCut = tilt[0].stats.starttime + fiveDaysInSec
    endCut = tilt[0].stats.endtime - fiveDaysInSec
    tilt = tilt.slice(starttime = UTCDateTime(startCut), endtime = UTCDateTime(endCut))
    
    return tilt
    
    
    
    