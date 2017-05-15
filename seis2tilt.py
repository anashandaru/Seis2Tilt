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

# Time cut for labuhan
#startCut = "20160220T000000.0"
#endCut = "20160311T235959.9"

# Time cut for grawah
#startCut = "20151005T000000.0"
#endCut = "20151025T235959.9"

# Time cut for klatakan
#startCut = "20160927T000000.0"
#endCut = "20161018T235959.9"

# Time cut for ijo
startCut = "20161205T000000.0"
endCut = "20161225T235959.9"

# Set seismic and tilt data file path respectively
ijonPath = r'/home/anas/Documents/BPPTKG/Seismometer2Tilt/Data/ijob-n/*.msd'
ijoTiltPath = r'/home/anas/Documents/BPPTKG/Seismometer2Tilt/Data/gunungijo.csv'

# Read seismic data and return to ijon
ijon = read(ijonPath)
ijon.resample(7); # resample seismic raw data
ijon.merge() # merge seismic data

#ijoz = read(ijozPath)
#ijoz.resample(7); # resample seismic raw data
#ijoz.merge() # merge seismic data

# Read tilt data and return to ijot
ijotRaw = open(ijoTiltPath,'r')
ijotRead = csv.reader(ijotRaw, delimiter=',')
ijotList = []
ijotTime = []
for row in ijotRead:
    if((float(row[2]) > 0.110) | (float(row[1]) == 0)):
        ijotList.append(0.106)
    else:
        ijotList.append(float(row[1]))
    ijotTime.append(UTCDateTime(row[0]))
    
stats = {'starttime':ijotTime[0],'delta':'60'}
print ijotRead
ijotArray = np.array(ijotList)
ijotTrace = Trace(data=ijotArray,header=stats)
#jotTrace.stats.delta = 60
ijot = Stream(ijotTrace)

#poles = [-0.149+0.149j,-0.149-0.149j,-503+0j,-1010+0j,-1130+0j]
poles = [-0.01178+0.01178j, -0.01178-0.01178j, -160+0j, -80+0j, -180+0j]
zeros = [0.0+0j, 0.0+0j]

# Find frequency response of the seismometer
h, f = paz_to_freq_resp(poles,zeros,0.4,0.005,16384,freq = True)

# Poles, Zeros, Gain, Sensitivity --> from Guralp Calibration Sheet
paz_cmg40T = {'poles': poles,
              'zeros': zeros,
              'gain': 2304000.0,
              'sensitivity': 2*2977.0*1e6}

# Create copy of seismic raw data
ijoN = ijon.copy();
#ijoZ = ijoz.copy();

# Apply Detrend
#ijon.detrend(type = 'polynomial',order=2);
#ijoz.detrend(type = 'polynomial',order=2);


# Apply Instrument Correction
ijon.simulate(paz_remove=paz_cmg40T,remove_sensitivity = True, taper = True, zero_mean=True,simulate_sensitivity=True)
ijon.filter('highpass', freq = 0.00000005, corners=2, zerophase=True)
ijon.filter('lowpass', freq = 0.0005, corners=2, zerophase=True)

#ijoz.simulate(paz_remove=paz_cmg40T,remove_sensitivity = True, taper = True, zero_mean=True,simulate_sensitivity=True)
#ijoz.filter('highpass', freq = 0.000055, corners=2, zerophase=True)
#ijoz.filter('lowpass', freq = 0.001, corners=2, zerophase=True)


# Create copy of each channel to get acceleration value
ijonA = ijon.copy()
#ijozA = ijoz.copy()


# Calculate Acceleration value
ijonA.differentiate(method = 'gradient')
#ijozA.differentiate(method = 'gradient')


# Calculate Tilt Value
# tilt =  np.sqrt(np.square(ijonA[0].data)+np.square(ijonA[0].data))/9.81
#tilt = ijonA[0].data/9.81
tilt = ijonA[0].data/-9.81
trace = Trace(data=tilt,header= ijonA[0].stats)
ijonTs = Stream(trace)

# cut all the data to specific time range
ijon = ijon.slice(starttime = UTCDateTime(startCut), endtime = UTCDateTime(endCut))
#ijoz = ijoz.slice(starttime = UTCDateTime(startCut), endtime = UTCDateTime(endCut))
ijonA = ijonA.slice(starttime = UTCDateTime(startCut), endtime = UTCDateTime(endCut))
ijonTs = ijonTs.slice(starttime = UTCDateTime(startCut), endtime = UTCDateTime(endCut))
ijot = ijot.slice(starttime = UTCDateTime(startCut), endtime = UTCDateTime(endCut))

# Plot Seismic Data
ijoN.plot()
print("\tVelocity Stream")

# Plot frequency response of the seismometer
plt.figure()
plt.loglog(f, abs(h))
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
print("\tFrequency response of the seismometer")

ijon.plot()
print("\tVelocity Instrument Corrected and band pass filtered")
ijonA.plot()
print("\tAcceleration Data")

ijonTs.plot()
print("\tTilt derived from velocity of N component")

# Plot tilt data
ijot.plot()
print("\tTilt from tiltmeter")
