import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import find_peaks

# wavelengths of x-rays
kcualpha = 154e-12
kcubeta = 138e-12

# %% NaCl - 5 seconds
NaClfilename = 'XRD-Data-csv/20-50-5s-NaCl.csv'
NaCltheta, NaClcounts = np.genfromtxt(NaClfilename, usecols=(0, 1), skip_header=1, unpack=True)

# detector measured 2*theta, so account for that
NaCltheta /= 2

NaClpeaks, _ = find_peaks(NaClcounts, prominence=15)
NaClpeaktheta = NaCltheta[NaClpeaks]

NaCld = np.zeros(2)
NaCld[0] = kcubeta/(2*np.sin(NaClpeaktheta[0] * np.pi/180))
NaCld[1] = kcubeta/(2*np.sin(NaClpeaktheta[1] * np.pi/180))
NaCld = np.average(NaCld)

print(f'The atomic spacing in NaCl is {NaCld*1e9:.3f} nm')
NaCllit = 5.64e-10 #m
NaClpercenterror = abs(NaCld*2 - NaCllit) / NaCllit * 100
print(f'The percent error for NaCl is {NaClpercenterror:0.3}%')

# %% LiF - 5 seconds
LiFfilename = 'XRD-Data-csv/20-50-5s-LiF.csv'
LiFtheta, LiFcounts = np.genfromtxt(LiFfilename, usecols=(0, 1), skip_header=1, unpack=True)

# detector measured 2*theta, so account for that
LiFtheta /= 2

LiFpeaks, _ = find_peaks(LiFcounts, prominence=15)
LiFpeaktheta = LiFtheta[LiFpeaks]

LiFd = np.zeros(2)
LiFd[0] = kcubeta/(2*np.sin(LiFpeaktheta[0] * np.pi/180))
LiFd[1] = kcubeta/(2*np.sin(LiFpeaktheta[1] * np.pi/180))
LiFd = np.average(LiFd)

print(f'The atomic spacing in LiF is {LiFd*1e9:.3f} nm')
LiFlit = 4.03e-10 #m
LiFpercenterror = abs(LiFd*2 - LiFlit) / LiFlit * 100
print(f'The percent error for LiF is {LiFpercenterror:0.3}%')

# %% RbCl - 5 seconds
RbClfilename = 'XRD-Data-csv/20-50-5s-RbCl.csv'
RbCltheta, RbClcounts = np.genfromtxt(RbClfilename, usecols=(0, 1), skip_header=1, unpack=True)

# detector measured 2*theta, so account for that
RbCltheta /= 2

RbClpeaks, _ = find_peaks(RbClcounts, prominence=15)
RbClpeaktheta = RbCltheta[RbClpeaks]

RbCld = np.zeros(2)
RbCld[0] = kcubeta/(2*np.sin(RbClpeaktheta[0] * np.pi/180))
RbCld[1] = kcubeta/(2*np.sin(RbClpeaktheta[1] * np.pi/180))
RbCld = np.average(RbCld)

print(f'The atomic spacing in RbCl is {RbCld*1e9:.3f} nm')
RbCllit = 6.59e-10 #m
RbClpercenterror = abs(RbCld*2 - RbCllit) / RbCllit * 100
print(f'The percent error for RbCl is {RbClpercenterror:0.3}%')

# %% Plotting
plt.figure(figsize=(7, 5), dpi=300)

#NaCl plot
plt.plot(NaCltheta, NaClcounts, linestyle='solid', color='red', label='NaCl')
plt.plot(NaCltheta[NaClpeaks], NaClcounts[NaClpeaks], 'x', color='black')

#LiF plot
plt.plot(LiFtheta, LiFcounts, linestyle='dotted', color='blue', label='LiF')
plt.plot(LiFtheta[LiFpeaks], LiFcounts[LiFpeaks], 'x', color='black')

#RbCl plot
plt.plot(RbCltheta, RbClcounts, linestyle='dashed', color='green', label='RbCl')
plt.plot(RbCltheta[RbClpeaks], RbClcounts[RbClpeaks], 'x', color='black')

plt.xlabel('Theta (degrees)')
plt.ylabel('Intensity (counts)')
plt.legend()
plt.title('XRD-Analysis')
plt.show()
