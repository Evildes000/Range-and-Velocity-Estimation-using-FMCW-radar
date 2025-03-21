import numpy as np 
import matplotlib.pyplot as plt

#Radar parameters setting
maxR = 200 #Maximum range
rangeRes = 1 #Range resolution 
maxV = 70 #Maximum speed
fc = 77e9 #Carrier frequency
c = 3e8 #Speed of light
lamda = c / fc 

r0 = 200 #Target distance
v0 = 50 #Target speed

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


# plot chirp and changing rate of frequecy
t_chirp = np.arange(0, Nr) / Fs 
angle_phase0 = 2 * np.pi * ( fc*t_chirp + slope * t_chirp * t_chirp/2)
Tx = np.cos(angle_phase0)
freq = fc + slope * t_chirp

plt.figure(0)
plt.subplot(2,1,1)
plt.plot(t_chirp, Tx)
plt.xlabel("t/s")
plt.ylabel("amp")
plt.title("Chirp")
plt.grid()

plt.subplot(2,1,2)
plt.plot(t_chirp, freq)
plt.xlabel("t/s")
plt.ylabel("freq")
plt.title("Freq vs Time")
plt.grid()
plt.show()


# estimate range
delay = 2 * r0 / c
dt = 2*v0*Tchirp/c
print(f"dt is: {dt}")
angle_phase1 = 2 * np.pi * ( fc*(t_chirp - delay) + slope * (t_chirp-delay) * (t_chirp-delay)/2)
# Rx = np,cos(phase)
angle_phase_diff = angle_phase1 - angle_phase0
IFx = np.cos(angle_phase_diff)

plt.figure(1)
plt.plot(t_chirp, IFx)
plt.xlabel("t/s")
plt.ylabel("amp")
plt.title("IFx")
plt.grid()
plt.show()


# plot fft of the first IFx
IFx_fft = np.abs(np.fft.fft(IFx))
threshold = np.max(IFx_fft)/10000
for i in np.arange(len(IFx_fft)):
    if IFx_fft[i] < threshold:
        IFx_fft[i] = 0

range =c * np.fft.fftfreq(len(IFx_fft), d=1/Fs) / (2 * slope)
plt.figure(2)
plt.plot(range, IFx_fft)
plt.xlabel("range m/s")
plt.ylabel("amp")
plt.title("range estimation")
plt.grid()
plt.show()


angle_phase_mat = np.zeros(shape=(Nd, Nr))
angle_phase_mat[0] = angle_phase1
for i in np.arange(1, Nd):
    increment_delay = delay + i * dt
    angle_phase_mat[i] = 2 * np.pi * ( fc*(t_chirp - increment_delay) + slope * (t_chirp-increment_delay) * (t_chirp-increment_delay)/2) 


angle_diff_mat = np.zeros(shape=(Nd, Nr)) # used to store angle phase difference 
angle_diff_mat[0] = angle_phase_diff
for i in np.arange(1, Nd):
    # angle_diff_mat[i] = angle_phase_mat[i] - angle_diff_mat[i-1]
    angle_diff_mat[i] = angle_phase_mat[i] - angle_phase0

IFx_mat = np.zeros(shape=(Nd, Nr)) # used to store IF signal
IFx_mat[0] = IFx
for i in np.arange(1,  Nd):
    IFx_mat[i] = np.cos(angle_diff_mat[i])


# compute range FFT
doppler_fft = np.abs(np.fft.fft2(IFx_mat))
factor =  (c / fc) / (4*np.pi*Tchirp)
print(f"the factor is: {factor}")
plt.figure(0)
plt.imshow(doppler_fft[0:64, 0:512])
plt.show()
