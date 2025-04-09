import numpy as np
import matplotlib.pyplot as plt


f1 = 500
f2 = 200

# sampling frequency 
fs = 1200

t = np.arange(0, 5, 1/fs)
print(f"length of t is: {len(t)}")
# signal generate
signal1 = np.cos(2*np.pi*f1*t + np.pi/3)
signal2 = np.cos(2*np.pi*f2*t + np.pi * 3/4)
signal = signal1 + signal2
# FFT
fft_signal = np.fft.fft(signal)
print(f"fft_siganl is: {fft_signal}")
fft_signal_abs = abs(fft_signal)
argmax_index = np.argmax(fft_signal_abs)
print(f"argmax index: {argmax_index}")
# fft_signal1_freq = np.fft.shift(np.fft.fftfreq(len(fft_signal1), d = 1/fs))
freq = np.fft.fftfreq(len(fft_signal), d = 1/fs) # doesn't need fftshift

threshold = np.max(fft_signal_abs)/10000
for i in np.arange(0, len(fft_signal_abs)):
    if fft_signal_abs[i] < threshold:
        fft_signal[i] = 0


# fft_signal2 = abs(np.fft.fft(signal2))
# fft_signal2_freq = np.fft.shift(np.fft.fftfreq(len(fft_signal2), d = 1/fs))

# plots frequency
plt.figure(0)
plt.plot(t, signal)
plt.title("time domin")
plt.xlabel("t/s")
plt.ylabel("amp")
plt.grid()
plt.show()

plt.figure(1)
plt.plot(freq, fft_signal_abs)
plt.title("FFT")
plt.xlabel("Hz")
plt.ylabel("amp")
# plt.legend(["signal1", "signal2"], loc = "upper right")
plt.grid()
plt.show()


plt.figure(2)
phase = np.angle(fft_signal) / np.pi
# new_signal = np.cos(phase)
print(phase)
plt.plot(freq, phase) 
plt.xlabel("f/hz")
plt.ylabel("phase/pi")
plt.title("phase")
plt.show()
