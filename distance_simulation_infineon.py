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

Nd = 64 #Number of chirp
Nr = 1024 #Numnber ADC sampling points
vres = (c/fc)/(2*Nd*(Tchirp+endle_time)) #Speed resolution 
Fs = Nr/Tchirp #Sampling rate
# Fs is lower than B(against Nyquist?)

print(f"Fs is: {Fs}")

#################################################################################
# generate Tx 
# t = np.linspace(0,Nd*Tchirp,Nr*Nd) #Time of Tx and Rx
t = np.arange(0, Nd*Nr)/Fs

# time period of one chrip
one_Tchirp = t[0:Nr]
# angle_freq = fc*t+(slope*t*t)/2 #Tx signal angle speed
angle_freq_chirp = fc*one_Tchirp + (slope*one_Tchirp*one_Tchirp)/2 # angle frequency of one chirp
# freq = fc + slope*t #Tx frequency
freq_chirp = fc + slope*one_Tchirp
one_chirp = np.cos(2*np.pi*angle_freq_chirp)
# Tx = np.cos(2*np.pi*angle_freq) #Waveform of Tx

Tx = np.array([]) # Tx contains all chirps
freq = np.array([]) # freq is the relation of frequency and time of all chirps 
angle_freq = np.array([]) # angle_freq is angle over periods of all chirps

for i in np.arange(0, Nd):
    # Tx.append(one_chirp)
    Tx = np.concatenate((Tx, one_chirp))
    # freq.append(freq_chirp)
    freq = np.concatenate((freq, freq_chirp))
    angle_freq = np.concatenate((angle_freq, angle_freq_chirp))

# plot chirp and change of frequency over one chirp period
print(f"length of Tx and freq are: {len(Tx)}, {len(freq)}")

plt.figure(0)
plt.subplot(2,1,1)
plt.plot(t[0:Nr], Tx[0:Nr])
plt.xlabel("t/s")
plt.ylabel("amp")
plt.title("Chirp")
plt.grid()

plt.subplot(2,1,2)
plt.plot(t[0:Nr], freq[0:Nr])
plt.xlabel("t/s")
plt.ylabel("freq")
plt.title("Freq vs Time")
plt.grid()

plt.show()

# compute Rx
td = 2*r0 / c
dt = 1/Fs
delay_index = round(td/dt) - 1
print(f"delay_index is: {delay_index}")
IF = slope*t[delay_index] 
print(f"IF is: {IF}")
Rx = Tx

plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t[0:Nr], Tx[0:Nr])
plt.plot((t+td)[0:Nr], Rx[0:Nr])
plt.legend(["Tx", "Rx"])
plt.xlabel("t/s")
plt.ylabel("amp")
plt.grid()

plt.subplot(2,1,2)
plt.plot(t[0:Nr], freq[0:Nr])
plt.plot((t+td)[0:Nr], freq[0:Nr])
plt.legend(["Tx", "Rx"])
plt.xlabel("t/s")
plt.ylabel("amp")
plt.grid()

plt.show()

# IF signal
IF_angle = angle_freq[delay_index:] - angle_freq[:-delay_index]
IF_signal = np.cos(2*np.pi*IF_angle)
IF_t = t[delay_index:]
plt.figure(2)
plt.plot(IF_t[0:1024], IF_signal[0:1024])
plt.title("IF_signal")
plt.xlabel("t/s")
plt.ylabel("amp")
plt.grid()
plt.show()

# FFT of IF signal
IF_fft = np.abs(np.fft.fft(IF_signal))
IF_fft_freq = c * np.fft.fftfreq(len(IF_fft), d=1/Fs) / (2*slope)
plt.figure(3)
plt.plot(IF_fft_freq, IF_fft)
plt.xlabel("range/m")
plt.ylabel("amp")
plt.title("FFT of IF Signal")
plt.grid()
plt.show()

# estimate the velocity




"""
plt.figure(0)
plt.plot(t[0:1024], Tx[0:1024])
plt.title("Tx")
plt.show()


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
"
"""