import numpy as np 
import matplotlib.pyplot as plt

"""
Info:
    construct cos across rows
Paras:
    mat: phase matrix
Return:
    a matrix with cosinus across the rows
"""
def cos_accross_rows(mat:np.ndarray)->np.ndarray:
    num_rows = len(mat[:,1])
    num_columns = len(mat[1,:])
    cos_mat = np.zeros(shape=(num_rows, num_columns))
    for i in np.arange(0, num_columns):
        cos_mat[:,i] = np.cos(mat[:,i])
    
    return cos_mat


"""
Info:
    Compute angle difference between rach Tx and Rx chirp
Paras:
    mat: angle phase matrix
Return:
    angle_diff_mat = a matix with phase difference
"""
def compute_phase_diff(mat:np.ndarray)->np.ndarray:
    num_rows = len(mat[:,1])
    num_columns = len(mat[1,:])
    angle_diff_mat = np.zeros(shape=(num_rows, num_columns))
    first_rows = mat[0]
    for i in np.arange(0, num_rows):
        angle_diff_mat[i] = mat[i] - first_rows 
    
    return angle_diff_mat


"""
Info:
    transfer result of fft into its phase
Paras:
    mat: range FFT
Return:
    phase_mat: matrix of phase
"""
def compute_phase(mat:np.ndarray)->np.ndarray:
    num_rows = len(mat[:,1])
    num_columns = len(mat[1,:])
    phase_mat = np.zeros(shape=(num_rows, num_columns))
    for i in np.arange(0, num_rows):
        phase_mat[i] = np.angle(mat[i])
    
    return phase_mat





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
delay = 2 * r0 / c # delay when the target is static
dt = 2*v0*Tchirp/c # addtional delay when the target moves within every chirp period
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


# copmpute fft of the first IFx
IFx_fft = np.fft.fft(IFx)
IFx_fft_abs = np.abs(IFx_fft)
range =c * np.fft.fftfreq(len(IFx_fft_abs), d=1/Fs) / (2 * slope)

# find the index of the peak value, which corresponds to range of the target
max_index = np.argmax(IFx_fft_abs)
print(f"index of max value is: {max_index}")
plt.figure(2)
plt.plot(range, IFx_fft_abs)
plt.xlabel("range/m")
plt.ylabel("amp")
plt.title("range estimation")
plt.grid()
plt.show()


angle_phase_mat = np.zeros(shape=(Nd, Nr)) # used to store the angle phase of Rx
angle_phase_mat[0] = angle_phase1
for i in np.arange(1, Nd):
    increment_delay = delay + i * dt
    angle_phase_mat[i] = 2 * np.pi * ( fc*(t_chirp + increment_delay) + slope * (t_chirp + increment_delay) * (t_chirp + increment_delay)/2) 


after_mixer_mat = np.zeros(shape=(Nd, Nr)) # used to store angle phase difference 
after_mixer_mat[0] = angle_phase_diff
for i in np.arange(1, Nd):
    # angle_diff_mat[i] = angle_phase_mat[i] - angle_diff_mat[i-1]
    after_mixer_mat[i] = angle_phase_mat[i] - angle_phase0


IFx_mat = np.zeros(shape=(Nd, Nr)) # used to store IF signal
IFx_mat[0] = IFx
for i in np.arange(1,  Nd):
    IFx_mat[i] = np.cos(after_mixer_mat[i])


# compute range FFT
IFx_fft_mat = np.fft.fft(IFx_mat, axis=1) # apply fft along columns of IF signal matrix
print(f"Tfx_fft_mat: {IFx_fft_mat}")
phase_mat = compute_phase(IFx_fft_mat)
print(f"phase_mat 200th column: {phase_mat}")


angle_diff_mat = compute_phase_diff(phase_mat)
# construct cosinus function along rows
cos_accross_rows_mat = cos_accross_rows(angle_diff_mat)


# apply fft across rows
plt.figure(3)
plt.plot(cos_accross_rows_mat[:,200])
plt.title("phase evaluation")
plt.xlabel("chirps")
plt.ylabel("phase")
plt.show()


# doppler FFT(apply FFT across the rows of the 200th column)
doppler_fft = np.abs(np.fft.fft(cos_accross_rows_mat[:,200])) / Nd
doppler_freq = lamda * np.fft.fftfreq(Nd, d=1) / (2*Tchirp)
# doppler_freq = 2*np.arange(0, Nd)

plt.figure(4)
plt.plot(doppler_freq[0:int(Nd/2)], doppler_fft[0:int(Nd/2)])
plt.title("doppler fft")
plt.xlabel("velocity")
plt.ylabel("amp")
plt.show()
