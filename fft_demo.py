import numpy as np
import matplotlib.pyplot as plt


f1 = 500
f2 = 200

# sampling frequency 
fs = 1200

t = np.arange(0, 5, 1/fs)
print(f"length of t is: {len(t)}")
# signal generate
signal1 = np.cos(2*np.pi*f1*t)
signal2 = np.cos(2*np.pi*f2*t)
signal = signal1 + signal2
# FFT
fft_signal = abs(np.fft.fft(signal)) / len(signal)
# fft_signal1_freq = np.fft.shift(np.fft.fftfreq(len(fft_signal1), d = 1/fs))
freq = np.fft.fftfreq(len(fft_signal), d = 1/fs) # doesn't need fftshift

# fft_signal2 = abs(np.fft.fft(signal2))
# fft_signal2_freq = np.fft.shift(np.fft.fftfreq(len(fft_signal2), d = 1/fs))

# plots
plt.figure(0)
plt.plot(t, signal)
plt.title("time domin")
plt.xlabel("t/s")
plt.ylabel("amp")
plt.grid()
plt.show()

plt.figure(1)
plt.plot(freq, fft_signal)
plt.title("FFT")
plt.xlabel("Hz")
plt.ylabel("amp")
# plt.legend(["signal1", "signal2"], loc = "upper right")
plt.grid()
plt.show()