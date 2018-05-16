import numpy as np
import matplotlib.pyplot as plt

def part1_1():
    dt = 1.0
    MLAC_data = np.genfromtxt("MLAC_data.txt")
    PHL_data = np.genfromtxt("PHL_data.txt")

    #assume data from one column continues in next column

    MLAC_data = MLAC_data.flatten()
    PHL_data = PHL_data.flatten()

    #TODO: understand why pad with zeros?
    MLAC_data = np.append(MLAC_data, np.zeros(len(MLAC_data)))
    PHL_data = np.append(PHL_data, np.zeros(len(PHL_data)))

    #TODO: understand fftshift
    MLAC_data_fft = np.fft.fftshift(np.fft.fft(MLAC_data))
    PHL_data_fft = np.fft.fftshift(np.fft.fft(PHL_data))

    cross_correlation_ft = PHL_data_fft.conjugate() * MLAC_data_fft
    cross_correlation = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(cross_correlation_ft / dt)))
    time_axis = np.fft.fftshift(np.fft.fftfreq(len(cross_correlation_ft), dt))

    #TODO: understand why multiply time with len(time_axis)
    plt.plot(time_axis * len(time_axis), cross_correlation)
    plt.xlim((-250,250))
    plt.xlabel("time, (s)")
    plt.ylabel("cross_correlation")
    plt.savefig("figures/part1_1a.png")
    plt.clf()

    cross_correlation_np = np.correlate(PHL_data, MLAC_data, 'same')
    plt.plot(cross_correlation_np)
    plt.xlabel("time, (s)")
    plt.ylabel("cross_correlation")
    plt.savefig("figures/part1_1b.png")

    return cross_correlation

def part1_2():
    dt = 1.0
    MLAC_data = np.genfromtxt("MLAC_data.txt")
    PHL_data = np.genfromtxt("PHL_data.txt")

    #assume data from one column continues in next column

    MLAC_data = MLAC_data.flatten()
    PHL_data = PHL_data.flatten()

    MLAC_data = np.append(MLAC_data, np.zeros(len(MLAC_data)))
    PHL_data = np.append(PHL_data, np.zeros(len(PHL_data)))

    #"bit conversion"
    MLAC_data = np.sign(MLAC_data)
    PHL_data = np.sign(PHL_data)

    MLAC_data_fft = np.fft.fftshift(np.fft.fft(MLAC_data))
    PHL_data_fft = np.fft.fftshift(np.fft.fft(PHL_data))

    cross_correlation_ft = PHL_data_fft.conjugate() * MLAC_data_fft
    cross_correlation = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(cross_correlation_ft / dt)))
    time_axis = np.fft.fftshift(np.fft.fftfreq(len(cross_correlation_ft), dt))

    plt.plot(time_axis * len(time_axis), cross_correlation)
    plt.xlim((-250,250))
    plt.xlabel("time, (s)")
    plt.ylabel("cross_correlation")
    plt.savefig("figures/part1_2a.png")
    plt.clf()


    #conservation of phase information?
    cross_correlation_part1_1 = part1_1()
    plt.plot(time_axis*len(time_axis), np.angle(cross_correlation_part1_1), label="without sign conversion")
    plt.plot(time_axis*len(time_axis), np.angle(cross_correlation), label= "with bit conversion")
    plt.xlim((-250, 250))
    plt.xlabel("time, (s)")
    plt.ylabel("cross_correlation")
    plt.legend()
    plt.savefig("figures/part1_2b.png")
    plt.clf

if __name__ == "__main__":
    part1_2()