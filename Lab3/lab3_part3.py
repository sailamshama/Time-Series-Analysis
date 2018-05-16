import numpy as np
import matplotlib.pyplot as plt
import cmath

from co2data import co2Values, time, co2Data
from lab3_part2 import  ratFilter, part2_1

def part3_1():

    p = np.polyfit(time, co2Data, 1)
    linear_fit  = time*p[0]+p[1]
    detrended = co2Data - linear_fit

    fig = plt.figure()
    plt.plot(time, co2Data, label = "original data")
    plt.plot(time, linear_fit, label = "linear fit")
    plt.plot(time, co2Data - linear_fit, label = "detrended data")
    plt.xlabel("time, t (yr)")
    plt.ylabel("Co2 levels (ppm)")
    plt.legend()
    fig.savefig("figures/part3_1.png")

    return linear_fit, detrended

def part3_2():
    a, b, c, B, C = part2_1()
    N = np.array([a, b, c])
    D = np.array([1, B, C])

    linear_fit, detrended_data = part3_1()

    y = ratFilter(N, D, detrended_data)
    y_trended = y + linear_fit

    fig = plt.figure()
    plt.plot(time, co2Data, label="original data")
    plt.plot(time, y_trended, label = "filtered data")
    plt.xlabel("time, t (yr)")
    plt.ylabel("Co2 levels (ppm)")
    plt.legend()
    fig.savefig("figures/part3_2.png")

    return y_trended

def part3_3():

    linear_fit, detrended = part3_1()
    f_s = 12.0
    dt = 1/f_s
    detrended_fft = np.fft.fft(detrended) * dt
    detrended_fft_shifted = np.fft.fftshift(detrended_fft)
    f_axis = np.fft.fftshift(np.fft.fftfreq(len(detrended), dt))

    fig = plt.figure()
    plt.plot(2 * cmath.pi * f_axis, np.angle(detrended_fft_shifted), label = "phase")
    plt.xlabel("frequency, f")
    plt.ylabel("phase")
    #plt.legend()
    fig.savefig("figures/part3_3a.png")

    fig = plt.figure()
    plt.plot(2 * cmath.pi * f_axis, np.absolute(detrended_fft_shifted), label="amplitude")
    # plt.legend()
    plt.xlabel("frequency,f")
    plt.ylabel("amplitude")
    fig.savefig("figures/part3_3b.png")

    fig = plt.figure()
    plt.plot(2 * cmath.pi * f_axis, np.absolute(detrended_fft_shifted), label = "amplitude")
    plt.xlabel("frequency,f")
    plt.ylabel("amplitude")
    plt.xlim((0,2.5))
    fig.savefig("figures/part3_3c.png")

    fig = plt.figure()
    plt.plot(2 * cmath.pi * f_axis, np.angle(detrended_fft_shifted), label="phase")
    plt.xlabel("frequency, f")
    plt.ylabel("phase")
    plt.xlim((0, 2.5))
    fig.savefig("figures/part3_3d.png")

    f_filtered_data = detrended_fft_shifted.copy()
    for indx in np.where(abs(f_axis) >= 0.9)[0]:
        f_filtered_data[indx] = 0

    #turn back to time domain
    f_filtered_data_t = np.fft.ifft( np.fft.fftshift(f_filtered_data) / dt)
    f_filtered_data_t_trended = f_filtered_data_t + linear_fit

    fig = plt.figure()
    plt.plot(time, f_filtered_data_t_trended)
    plt.xlabel("time, t")
    plt.ylabel("f-filtered Co2 data")
    fig.savefig("figures/part3_3e.png")

    return f_filtered_data_t_trended

def part3_4():
    notch_filtered_data = part3_2()
    f_filtered_data = part3_3()

    fig = plt.figure(figsize=(15,10))
    plt.plot(time, co2Data, label="original data")
    plt.plot(time, notch_filtered_data, label = "notch filtered data")
    plt.plot(time, f_filtered_data, label = "f_filtered data")
    plt.xlabel("time, t (yr)")
    plt.ylabel("Co2 levels (ppm)")
    plt.legend()
    plt.title("Analysis done by first detrending data")
    fig.savefig("figures/part3_4.png")

def part3_5():

    #### part 2 method ####
    a, b, c, B, C = part2_1()
    N = np.array([a, b, c])
    D = np.array([1, B, C])

    notch_filtered_data = ratFilter(N, D, co2Data)

    #### part 3 method ####
    f_s = 12.0
    dt = 1/f_s
    fft = np.fft.fft(co2Data) * dt
    fft_shifted = np.fft.fftshift(fft)
    f_axis = np.fft.fftshift(np.fft.fftfreq(len(co2Data), dt))

    f_filtered_data = fft_shifted.copy()
    for indx in np.where(abs(f_axis) >= 0.9)[0]:
        f_filtered_data[indx] = 0

    #turn back to time domain
    f_filtered_data_t = np.fft.ifft( np.fft.fftshift(f_filtered_data) / dt)

    fig = plt.figure(figsize=(10,7))
    plt.plot(time, co2Data, label="original data")
    plt.plot(time, notch_filtered_data, label="notch filtered data")
    plt.plot(time, f_filtered_data_t, label="f_filtered data")
    plt.xlabel("time, t (yr)")
    plt.ylabel("Co2 levels (ppm)")
    plt.legend()
    plt.title("Analysis without detrending data")
    fig.savefig("figures/part3_5.png")
####################################################################################################################################
####################################################################################################################################
part3_1()
part3_2()
part3_3()
part3_4()
part3_5()