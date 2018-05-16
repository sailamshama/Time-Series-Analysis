import numpy as np
import matplotlib.pyplot as plt
import scipy.signal


def part2_1():
    data = np.genfromtxt("nwao.vh1")
    plt.plot(data[:,0] / 3600, data[:,1])
    plt.xlabel("time (hours)")
    plt.ylabel("Vertical Displacement")
    plt.savefig("figures/part2_1.png")

def part2_2():
    dt = 10
    data = np.genfromtxt("nwao.vh1")

    velocity = data[:,1]
    velocity_ft = np.fft.fftshift(np.fft.fft(velocity)*dt)
    power_spectrum = np.abs(velocity_ft) ** 2
    freq_axis = np.fft.fftshift(np.fft.fftfreq(len(velocity), dt))*1000.0

    plt.plot(freq_axis, power_spectrum)
    plt.xlabel("frequency")
    plt.ylabel("Power spectrum")
    plt.savefig("figures/part2_2.png")
    plt.clf()

    return power_spectrum, freq_axis

def part2_3():
    dt = 10
    data = np.genfromtxt("nwao.vh1")

    #detrend
    velocity = data[:, 1]
    velocity_detrended = scipy.signal.detrend(velocity)
    #apply hanning function
    n = np.arange(len(velocity_detrended))
    hanning = np.subtract(1, np.cos((2.0*np.pi*n/len(velocity_detrended))))

    velocity_detrended_windowed = np.multiply(velocity_detrended, hanning)
    velocity_detrended_windowed_ft = np.fft.fftshift(np.fft.fft(velocity_detrended_windowed) * dt)
    power_spectrum = np.abs(velocity_detrended_windowed_ft) ** 2
    freq_axis = np.fft.fftshift(np.fft.fftfreq(len(velocity), dt)) * 1000

    plt.plot(freq_axis, power_spectrum)
    plt.xlabel("frequency")
    plt.ylabel("Power spectrum")
    plt.savefig("figures/part2_3.png")
    plt.clf()
    return power_spectrum, freq_axis

def part2_4():
    power_spectrum_unwindowed, freq_axis_unwindowed = part2_2()
    power_spectrum_windowed, freq_axis_windowed = part2_3()

    #xlim(0.1,2.2) doesn't work
    start = np.where(freq_axis_windowed >= 0.1)[0][0]
    stop = np.where(freq_axis_windowed >= 2.2)[0][0]
    plt.plot(freq_axis_unwindowed[start:stop], power_spectrum_unwindowed[start:stop], label="unwindowed")
    plt.plot(freq_axis_windowed[start:stop], power_spectrum_windowed[start:stop], label="windowed")
    plt.legend()
    plt.savefig("figures/part2_4.png")
    plt.show()

def part2_5():
    power_spectrum_unwindowed, freq_axis_unwindowed = part2_2()
    power_spectrum_windowed, freq_axis_windowed = part2_3()
    start = np.where(freq_axis_windowed >= 0.1)[0][0]
    stop = np.where(freq_axis_windowed >= 2.2)[0][0]

    np.savetxt("power_spectrum_unwindowed.txt", power_spectrum_unwindowed[start:stop])
    np.savetxt("power_spectrum_windowed.txt", power_spectrum_windowed[start:stop])
    np.savetxt("freq_axis.txt",freq_axis_windowed[start:stop])

    #annotation done in plot_part2_5.m matlab script" 

if __name__ == "__main__":
    part2_5()