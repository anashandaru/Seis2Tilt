# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 09:28:07 2017

@author: anas
"""

from obspy.core import read, Trace, Stream, UTCDateTime
from obspy.signal.invsim import paz_to_freq_resp
import numpy as np
import csv
import matplotlib.pyplot as plt

ijozPath = r'/home/anas/Documents/BPPTKG/Seismometer2Tilt/IJO/ijob-z/*.msd'
ijonPath = r'/home/anas/Documents/BPPTKG/Seismometer2Tilt/IJO/ijob-n/*.msd'
ijoePath = r'/home/anas/Documents/BPPTKG/Seismometer2Tilt/IJO/ijob-e/*.msd'
ijoTiltPath = r'/home/anas/Documents/BPPTKG/Seismometer2Tilt/IJO/gunungijo.csv'

ijoz = read(ijozPath)
ijon = read(ijonPath)
ijoe = read(ijoePath)
ijotRaw = open(ijoTiltPath,'r')
ijotRead = csv.reader(ijotRaw, delimiter=',')
ijotList = []
ijotTime = []
for row in ijotRead:
    if(float(row[1]) == 0):
        ijotList.append(0.10940)
    else:
        ijotList.append(float(row[1]))
    ijotTime.append(UTCDateTime(row[0]))
    
stats = {'starttime':ijotTime[0],'delta':'60'}
print ijotRead
ijotArray = np.array(ijotList)
ijotTrace = Trace(data=ijotArray,header=stats)
#jotTrace.stats.delta = 60
ijot = Stream(ijotTrace)
ijot.plot()
#print ijotRead[0][0]
#ijotStats.starttime = ijotList[0].strip().split()[0]


ijoz.merge()
ijon.merge()
ijoe.merge()

print ijon[0].stats
print ijot[0].stats
#ijoz.plot()
#ijon.plot()
#ijoe.plot()

#poles = [-0.149+0.149j,-0.149-0.149j,-503+0j,-1010+0j,-1130+0j]
poles = [-0.01178+0.01178j, -0.01178-0.01178j, -160+0j, -80+0j, -180+0j]
zeros = [0.0+0j, 0.0+0j]

# Find frequency response of the seismometer
h, f = paz_to_freq_resp(poles,zeros,0.4,0.005,16384,freq = True)

# Plot frequency response of the seismometer
#plt.figure()
#plt.loglog(f, abs(h))
#plt.xlabel('Frequency [Hz]')
#plt.ylabel('Amplitude')

# Poles, Zeros, Gain, Sensitivity --> from Guralp Calibration Sheet
paz_cmg40T = {'poles': poles,
              'zeros': zeros,
              'gain': 2304000.0,
              'sensitivity': 2*2977.0*1e6}

# Apply Instrument Correction to each channel
ijoz.simulate(paz_remove=paz_cmg40T,remove_sensitivity = True, taper = True, zero_mean=True,simulate_sensitivity=True)
ijon.simulate(paz_remove=paz_cmg40T,remove_sensitivity = True, taper = True, zero_mean=True,simulate_sensitivity=True)
ijoe.simulate(paz_remove=paz_cmg40T,remove_sensitivity = True, taper = True, zero_mean=True,simulate_sensitivity=True)

ijon.filter('highpass', freq = 0.00005, corners=2, zerophase=True)
ijon.filter('lowpass', freq = 0.001, corners=2, zerophase=True)

#ijoz.plot()
ijon.plot()
#ijoe.plot()

# Create copy of each channel to get acceleration value
#ijozA = ijoz.copy()
ijonA = ijon.copy()
#ijoeA = ijoe.copy()

# Calculate Acceleration value
#ijozA.differentiate(method = 'gradient')
ijonA.differentiate(method = 'gradient')
#ijoeA.differentiate(method = 'gradient')

#tilt = ijonA[0].data/9.81
#trace = Trace(tilt)
#stream = Stream(trace)

#ijozA.plot()
ijonA.plot()
#ijoeA.plot()

#stream.plot()