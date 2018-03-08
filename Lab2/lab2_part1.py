import numpy as np
import matplotlib.pyplot as plt
import cmath

def gaussian(t, center, tH):
    gaussian = ( 1 / ( np.sqrt(np.pi) * tH) ) * np.exp( - ( (t-center) /tH) **2 )
    return gaussian

def gaussian_FT_analytical(f, tH):
    gaussian_ft = np.exp ( - (f * tH / 2.0 ) ** 2  )
    return gaussian_ft

def part1_1():
    center = 0
    t = np.arange(-100, 100, 1)

    gaussian_th10 = gaussian(t, center, 10)
    gaussian_th20 = gaussian(t, center, 20)

    fig = plt.figure(figsize=(7,7))
    plt.grid(1)
    plt.plot(t, gaussian_th10, label = "tH = 10")
    plt.plot(t, gaussian_th20, label = "tH = 20")

    plt.title("Gaussian distribution functions with half durations 10 and 20")
    plt.xlabel("t")
    plt.ylabel("gaussian functions, g(t)")
    plt.ylim((-0.01,0.06))
    plt.legend()
    fig.savefig("figures/part1_1.png")

def part1_3():
    T = 100
    dt = 1
    t = np.arange(-T, T, dt)
    center = 0

    #th = 10
    gaussian_fn10 = gaussian(t, center, 10)
    gaussian_fn10_fft = np.fft.fft(gaussian_fn10)
    gaussian_fn10_fft_shifted = np.fft.fftshift(gaussian_fn10_fft)
    f_axis_shifted10 = np.fft.fftshift(np.fft.fftfreq(len(gaussian_fn10), dt))
    gaussian_fn10_analytical = gaussian_FT_analytical(2 * cmath.pi * f_axis_shifted10, 10)
    gaussian_fn10_fft_shifted_amplitude = np.absolute(gaussian_fn10_fft_shifted)
    gaussian_fn10_analytical_amplitude = np.absolute(gaussian_fn10_analytical)

    #th = 20
    gaussian_fn20 = gaussian(t, center, 20)
    gaussian_fn20_fft = np.fft.fft(gaussian_fn20)
    gaussian_fn20_fft_shifted = np.fft.fftshift(gaussian_fn20_fft)
    f_axis_shifted20 = np.fft.fftshift(np.fft.fftfreq(len(gaussian_fn20), dt))
    gaussian_fn20_analytical = gaussian_FT_analytical(2 * cmath.pi * f_axis_shifted20, 20)
    gaussian_fn20_fft_shifted_amplitude = np.absolute(gaussian_fn20_fft_shifted)
    gaussian_fn20_analytical_amplitude = np.absolute(gaussian_fn20_analytical)

    fig = plt.figure(figsize=(7,7))

    plt.plot(f_axis_shifted10, gaussian_fn10_analytical_amplitude, label ="th=10, analytical")
    plt.plot(f_axis_shifted10, gaussian_fn10_fft_shifted_amplitude, marker = "x", linestyle = "None", label = "th=10, calculated")
    plt.plot(f_axis_shifted20, gaussian_fn20_analytical_amplitude, label ="th=20, analytical")
    plt.plot(f_axis_shifted20, gaussian_fn20_fft_shifted_amplitude, marker = "x", linestyle = "None", label = "th=20, calculated")

    plt.xlabel("frequency, w")
    plt.ylabel("G(w)")
    plt.legend()
    plt.title("Comparison of analytical and FFT G(w)")

    fig.savefig("figures/part1_3.png")

def part1_4():
    #TODO: split this task and obtain only bottom right graphs
    dt = 0.01
    t = np.arange(-250, 250, dt)
    center = 0
    general_fn = np.sin(t)

    for th in [10, 20]:
        gaussian_fn = gaussian(t, center, th)
        general_fn_fft = np.fft.fft(general_fn)
        general_fn_fft_axis = np.fft.fftfreq(len(general_fn), dt)

        filtered = np.convolve(general_fn, gaussian_fn)
        filtered_axis = np.linspace(-5, 5, filtered.shape[0])
        filtered_fft = np.fft.fft(filtered)
        filtered_fft_axis = np.fft.fftfreq(len(filtered), dt)

        fig = plt.figure(figsize=(15, 15))

        a = fig.add_subplot(2, 2, 1)
        plt.plot(t, general_fn)
        plt.xlabel("time")
        plt.ylabel("f(t)")
        plt.title("f(t) = sin(t)")

        b = fig.add_subplot(2, 2, 2)
        plt.plot(general_fn_fft_axis, np.absolute(general_fn_fft))
        plt.title("FFT Before filtering")
        plt.xlabel("frequency")
        plt.ylabel("amplitude")
        plt.xlim((-0.5, 0.5))

        c = fig.add_subplot(2, 2, 3)
        plt.plot(filtered_axis, filtered)
        plt.xlabel("t")
        plt.ylabel("f(t) * g(t)")
        plt.title("Convolution of f(t) with gaussian")

        d = fig.add_subplot(2, 2, 4)
        plt.plot(filtered_fft_axis, np.absolute(filtered_fft))
        plt.title("FFT after filtering")
        plt.xlabel("frequency")
        plt.ylabel("amplitude")
        plt.xlim((-0.3, 0.3))

        fig.suptitle("Filtering with gaussian with th = %s" % th)
        plt.savefig("figures/part1_4_tH" + str(th) + ".png")



################## MAIN ##################

if __name__ == "__main__":
    part1_4()

