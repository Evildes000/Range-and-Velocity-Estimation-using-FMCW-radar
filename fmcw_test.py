import numpy as np 
import matplotlib.pyplot as plt

#Radar parameters setting
maxR = 200 #Maximum range
rangeRes = 1 #Range resolution 
maxV = 70 #Maximum speed
fc = 77e9 #Carrier frequency
c = 3e8 #Speed of light
lamda = c / fc 

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
angle_phase1 = 2 * np.pi * ( fc*(t_chirp + delay) + slope * (t_chirp+delay) * (t_chirp+delay)/2)
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
# incresment = delay + dt
for i in np.arange(1, Nd):
    # for each chirp, dt is a constant number
    incresment = delay + i*dt
    angle_phase_mat[i] = 2 * np.pi * ( fc*(t_chirp + incresment) + slope * (t_chirp+incresment) * (t_chirp+incresment)/2)
    
# calculate FFT for each chirp and its reflected chirp
IF_mat = np.zeros(shape=(Nd,Nr))
IF_signal_mat = np.zeros(shape=(Nd,Nr))
FFT_signal_mat = np.zeros(shape=(Nd,Nr), dtype=np.complex64)
for i in np.arange(0, Nd):
    IF_mat[i] = angle_phase0 - angle_phase_mat[i]
    IF_signal_mat[i] = np.cos(IF_mat[i])
    FFT_signal_mat[i] = np.fft.fft(IF_signal_mat[i])



# check if the range is the same when using other chirps
"""
IFx_fft = np.abs(FFT_signal_mat[20])
range =c * np.fft.fftfreq(len(IFx_fft), d=1/Fs) / (2 * slope)
plt.figure(3)
plt.plot(range, IFx_fft)
plt.xlabel("range m/s")
plt.ylabel("amp")
plt.title("range estimation")
plt.grid()
plt.show()
"""
# compute velocity
doppler_fft = np.abs(np.fft.fft2(IF_signal_mat))
plt.figure(3)
plt.imshow(doppler_fft)
plt.show()
# compute FFT for rach row

"""
angle_phase_array = np.zeros(shape=(Nd,Nr))
angle_phase_array[0] = angle_phase0 
for i in np.arange(1, Nd):
    incresment = i * Tchirp
    angle_phase_array[i] = 2 * np.pi * ( fc*(t_chirp + incresment) + slope * (t_chirp+incresment) * (t_chirp+incresment)/2)"

# angle phase array of Rx
delay_array = np.array([delay])

for i in np.arange(1,Nd):
    delay = delay + i*dt
    delay_array = np.concatenate((delay_array, t_chirp+delay)) 

# calculate phase for each delay
phase_mat = np.zeros(shape=(Nd, Nr)) # mat is used to store phase of each 
for i in np.arange(0, Nd):
    phase_buffer = np.array(2 * np.pi * ( fc*(t_chirp + delay_array[i]) + slope * (t_chirp+delay_array[i]) * (t_chirp+delay_array[i])/2))
    print(f"type of phase_buffer is: {type(phase_buffer)}")
    # phase_buffer = np.transpose(phase_buffer)
    # print(f"shape of phase_buffer: {phase_buffer.shape}")
    # phase_mat = np.concatenate((phase_mat, phase_buffer), axis=0)
    # numpy的axis, along rows 叠加每个相位
    # phase_mat = np.vstack((phase_mat, phase_buffer))
    phase_mat[i] = phase_buffer

print(f"shape of phase_mat: {phase_mat.shape}")

# calculate phase differnce
phase_diff_mat = np.zeros(shape=(Nd, Nr))
# Fx_mat = np.zeros(shape=(Nd, Nr))
for i in np.arange(0, Nd):
    # phase_diff_buffer = phase_buffer[i] - phase0
    # phase_diff_mat = np.concatenate((phase_diff_mat, phase_diff_buffer), axis=0)
    phase_diff_mat[i] = phase_buffer[i] - phase0
    # IFx_mat[i] = np.cos(phase_diff_mat[i])


# take one column data to estimate velocity
first = np.abs(np.fft.fft(np.cos(phase_diff_mat[:,1])))
# velocty = np.arange(0, Nd) * lamda / (4*np.pi*Tchirp)
velocty = phase_diff_mat[1,1] * (lamda * np.fft.fftfreq(len(first), d = Tchirp)) / (4*np.pi*Tchirp)
# Z_fft2 = abs(np.fft.fft2())
plt.figure(3)
plt.plot(velocty, first)
# plt.imshow(Z_fft2[0:64, 0:512])
plt.xlabel("m/s")
plt.ylabel("amp")
plt.title("velocity")
plt.grid()
plt.show()


"""