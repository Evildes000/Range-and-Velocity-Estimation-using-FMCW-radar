import numpy as np
import matplotlib.pyplot as plt

range = 100 # meters
velocity = 100 # m/s
c = 3e8 # m/s
range_resolution = 1 # m


target_range = 5000 # m
target_velocity = 3 # m/s
t_sweep = 0.5 # s
slope = 10
B = slope*t_sweep

fc = 1e3 # carrier frequency
fs = 3 * 1e3
N_chirps = 1
N_samples_per_chirp = 1024

t = np.arange(0, t_sweep, 1/fs, dtype=np.float32)
print(t)
# signal generate
Tx = np.cos(np.pi * slope * t**2)
f = slope*t

plt.figure(0)
plt.subplot(211)
plt.plot(t,Tx)
plt.title("modulated cosinus")
plt.xlabel("t/s")
plt.ylabel("amp.")
plt.grid()

plt.subplot(212)
plt.plot(t, f)
plt.title("freq. vs time")
plt.xlabel("t/s")
plt.ylabel("freq.")
plt.grid()

plt.show()

# in spectrum
N = 1024 # number of DFT points
Tx_fft = np.abs(np.fft.fft(Tx))**2

freq = np.arange(0, len(Tx_fft)) * fs

plt.figure(1)
plt.plot(freq, Tx_fft)
plt.xlabel("frequency/hz")
plt.ylabel("power")
plt.grid()
plt.show()

# Rx signal
td = 2*target_range / c
Rx = np.cos(np.pi*slope*(t-td)**2)
print(t-td)
# Tx and Rx signal
plt.figure(2)
plt.plot(t-td, Rx, label = "Rx")
plt.grid()
# plt.legend(["Tx", "Rx"], loc = "upper right")
plt.title("Tx and Rx")
plt.show()

# IF 
IF = slope * td 
print(f"IF is: {IF}")
IF_signal = np.cos(2*np.pi*IF*(t-td))

# using frequency to find out range
IF_fft = abs(np.fft.fft(IF_signal)) ** 2
IF_freq = np.fft.fftshift(np.fft.fftfreq(len(IF_fft), d=1/fs)) / slope
plt.figure(3)
plt.plot(IF_freq, IF_fft)
plt.title("FFT of received Signal")
plt.grid()
plt.xlabel("range/m")
plt.ylabel("power")
plt.show()














