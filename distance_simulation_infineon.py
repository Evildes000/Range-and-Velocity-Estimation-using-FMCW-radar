import numpy as np 
import matplotlib.pyplot as plt

#Radar parameters setting
maxR = 200 #Maximum range
rangeRes = 1 #Range resolution 
maxV = 70 #Maximum speed
fc = 77e9 #Carrier frequency
c = 3e8 #Speed of light

r0 = 100 #Target distance
v0 = 70 #Target speed

B = c/(2*rangeRes) #Bandwidth
print(f"bandwidth is: {B}")
Tchirp = 5.5*2*maxR/c #Chirp time
print(f"Tchirp is: {Tchirp}")
endle_time = 6.3e-6 
slope = B/Tchirp #Chirp slope
print(f"slope is: {slope}")

f_IFmax = (slope*2*maxR)/c #Maximum IF frequency
f_IF = (slope*2*r0)/c #Current IF frequency

Nd = 128 #Number of chirp
Nr = 1024 #Numnber ADC sampling points
vres = (c/fc)/(2*Nd*(Tchirp+endle_time)) #Speed resolution 
Fs = Nr/Tchirp #Sampling rate
# Fs is lower than B(against Nyquist?)

print(f"Fs is: {Fs}")

# generate Tx 
# t = np.linspace(0,Nd*Tchirp,Nr*Nd) #Time of Tx and Rx
t = np.linspace(0,Nd*Tchirp,Nr*Nd)
print(f"t is: {t}")
angle_freq = fc*t+(slope*t*t)/2 #Tx signal angle speed
freq = fc + slope*t #Tx frequency
Tx = np.cos(2*np.pi*angle_freq) #Waveform of Tx

"""
plt.figure(0)
plt.plot(t[0:1024], Tx[0:1024])
plt.title("Tx")
plt.show()
"""

plt.figure(0)
plt.plot(t, freq)
plt.show()
# generate Rx
td = 2*r0/c
# freqRx = fc + slope*(t)
Rx = np.cos(2*np.pi*(fc*(t+td) + (slope*(t+td)*(t+td))/2))
print(f"t-td is: {t+td}")
plt.figure(0)
plt.plot((t+td)[0:1024], Rx[0:1024])
plt.title("Rx")
plt.show()

# calculate FFT of IF signal
IF_angle_freq = fc*t+(slope*t*t)/2 - ((fc*(t+td) + (slope*(t+td)*(t+td))/2))
freqIF = slope*td
print(f"freqIF is: {freqIF}")
# IPx is the signal after mixer
# IFx = np.cos(-(2*np.pi*(fc*(t-td) + (slope*(t-td)*(t-td))/2))+(2*np.pi*angle_freq))
IFx = np.cos(2*np.pi*IF_angle_freq)
IF_fft = abs(np.fft.fft(IFx))
IF_freq = (np.fft.fftfreq(len(IF_fft), d=t[1]-t[0]) / slope) * c / 2

plt.figure(1)
plt.plot((t+td)[0:1024], IFx[0:1024])
plt.xlabel("t/s")
plt.ylabel("amp")
plt.title("IFx")
plt.grid()
plt.show()

plt.figure(2)
plt.plot(IF_freq, IF_fft)
plt.xlabel("range")
plt.ylabel("amp.")
plt.grid()
plt.show()


# 2D plot
plt.figure(3)
plt.specgram(IFx, 1024, Fs)
plt.xlabel("time")
plt.ylabel("freq")
plt.title("Specgram")
plt.show()


####################### velocity estimation ###########################

# take one sample from each chirp
# we have 128 chirps, so there will be 128 samples
chirpamp = []
chirpnum = 1

# dt = t[1] - t[0]
# start_index = round(td / dt) 
# print(f"start index is: {start_index}")
while(chirpnum<=Nd):
    # strat = start_index + (chirpnum-1)*1024
    start = (chirpnum-1)*1024
    end = chirpnum*1024
    # why do we need one sample from each chirp?   
    chirpamp.append(IFx[(chirpnum-1)*1024]) 
    chirpnum = chirpnum + 1

print(f"length of chirpamp: {len(chirpamp)}")
print(f"the number is :{(c / fc) / (4*np.pi*Tchirp)  }")

velocities = 3 * np.arange(0, Nd) / 5

# Doppler FFT
# velocities = c * np.fft.fftfreq(Nd, d=1/Fs)/ (2 * fc * Tchirp)
amp = np.fft.fft(chirpamp)

plt.figure(4)
plt.plot(chirpamp)
plt.title("velocity")
plt.xlabel("v m/s ")
plt.ylabel("amp")
plt.grid()
plt.show()
