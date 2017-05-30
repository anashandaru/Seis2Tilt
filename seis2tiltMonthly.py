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
#monthFolder  = '02\\'
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

def RemoveGcf(thePath):
    print('remove gcf in' + thePath)
    system("del \"" + thePath)
    
def addZeroF(number):
    if(number<10):
        result = '0' + str(number) + '\\'
    else:
        result = str(number) + '\\'
    return result

def BatchMsdConverting(month):
    monthFolder = addZeroF(month)
    
    for station in stations:
        print ('CONVERTING --------> ' + station)
        for component in components:
            print ('Component ' + component)
            for i in range(30):                
                filePath = yearFolder + monthFolder + addZeroF(i+1) + station + component + gcfExt # Define filepath
                Gcf2msd(filePath) # Begin Converting
                
def BatchRemoveGcf(month):
    monthFolder = addZeroF(month)
    
    for station in stations:
        print ('Removing GCF --------> ' + station)
        for component in components:
            print ('Component ' + component)
            for i in range(30):
                filePath = yearFolder + monthFolder + addZeroF(i+1) + station + component + gcfExt # Define filepath
                RemoveGcf(filePath) # Begin Converting
        
def ReadMonthly(month,station):
    monthFolder = addZeroF(month)
    
    print monthFolder + station
    resultStream = Stream()
    for component in components:
        stream = Stream()
        for i in range(30):   
            filePath = yearFolder + monthFolder + addZeroF(i+1) + station + component + msdExt # Define filepath
            if(path.exists(yearFolder + monthFolder + addZeroF(i+1) + station + component)):
                #print('Reading... ' + filePath)
                print'#',
                newStream = read(filePath)
                newStream.resample(1)
                stream += newStream
            else:
                #print(filePath + 'Not Found')
                print 'x',                
        
        print('MERGING--------' + component)
        stream.merge(fill_value = 'interpolate') # merge seismic data
        resultStream += stream
        for trace in resultStream:
            trace.stats.station = station
    return resultStream

def GetTilt(stream):
    print('Clone Stream')
    #cstream = stream.copy()
    
    print('Perform instrument correction')
    stream.simulate(paz_remove=paz_cmg40T,remove_sensitivity = True, taper = True, zero_mean=True,simulate_sensitivity=True)
    
    print('Perform filering')
    stream.filter('highpass', freq = 0.00000005, corners=2, zerophase=True)
    stream.filter('lowpass', freq = 0.005, corners=2, zerophase=True)
    
    print('Perform diferentiation')
    stream.differentiate(method = 'gradient')
    
    print('Calculate Tilt')
    tilt = Stream()
    for traceComp in stream:
        # Calculate Tilt Value
        tiltVal = traceComp.data/-9.81
        trace = Trace(data=tiltVal,header= traceComp.stats)
        tilt += Stream(trace)
    
    print('Cut Tilt Edge')
    fiveDaysInSec = 1*1*10*60
    startCut = tilt[0].stats.starttime + fiveDaysInSec
    endCut = tilt[0].stats.endtime - fiveDaysInSec
    tilt.trim(starttime = UTCDateTime(startCut), endtime = UTCDateTime(endCut))
    return tilt
    
def TiltThisMonth(month):
    for station in stations:
        stream = ReadMonthly(month,station)
        if(len(stream)>0):
            tilt = GetTilt(stream)
            tilt.plot()
    
def TiltMonths(maximumMonth, station):
    stream = Stream()
    for i in range(maximumMonth):
        stream += ReadMonthly(i+1,station)
        stream.merge(fill_value = 'interpolate')
    if(len(stream)>0):
            tilt = GetTilt(stream)
            tilt.plot()
            
def TiltMonthsM(maximumMonth, station): #failed
    tilt = Stream()
    for i in range(maximumMonth):
        stream = ReadMonthly(i+1,station)
        tiltNew = GetTilt(stream)
        if(i>0):
            for i in range(len(tilt)):                
                tiltTail = tilt[i].data[len(tilt[i].data)-1]
                tiltNewHead = tiltNew[i].data[len(tiltNew[i].data)-1]
                offset = tiltTail - tiltNewHead
                tiltNew[i].data += offset
        tilt += tiltNew
    tilt.merge(fill_value = 'interpolate')
    tilt.plot()

        
        
    
    
    
    
    
    
    