import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import detrend, correlate

def get_fft(x, dt):
    x_ft = np.fft.fftshift(np.fft.fft(x)*dt)
    freq_axis = np.fft.fftshift(np.fft.fftfreq(len(x), dt))

    return freq_axis, x_ft


def get_ifft(x, dt):
    ifft = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(x / dt)))
    time_axis = np.fft.fftshift(np.fft.fftfreq(len(x), dt))
    return time_axis, ifft

def remove_nan(x):
    xi = np.arange(len(x))

    mask = np.isfinite(x)
    xfiltered = np.interp(xi, xi[mask], x[mask])

    return xfiltered


if __name__ == "__main__":

    ############################ LOAD DATA ##################################

    rain_data = np.genfromtxt("melbrain.txt") #3653 observations
    temp_data_max = np.genfromtxt("melbtemp.txt")[:,0] #10 x 365 = 3650 observations.
    temp_data_min = np.genfromtxt("melbtemp.txt")[:,1]


    ########################### PROCESS DATA #################################

    #fix dimension problem for cross-correlation
    # delete leap year days between in 1984 and 1988 from rain data
    # remove indices 365*3+31+29 = 1155, 365*7+1+31+29 = 2616
    # assume last data point is extra so remove that as well
    rain_data = np.delete(rain_data, [1155,2616, len(rain_data)-1])

    #remove Nan by linear extrapolation so that FT doesn't all give NaN values
    rain_data = remove_nan(rain_data)
    temp_data_max = remove_nan(temp_data_max)
    temp_data_min = remove_nan(temp_data_min)

    np.savetxt("processed_rain_data.txt", rain_data)

    ########################### PLOT RAW DATA #################################

    plt.plot(rain_data, label="rainfall")
    plt.xlabel("time (days since 1st Jan, 1980)")
    plt.ylabel("Rainfall Amount (mm)")
    plt.savefig("figures/raw_data_rain.png")
    plt.clf()

    plt.plot(temp_data_max, label="max temp")
    plt.plot(temp_data_min, label="min temp")
    plt.legend()
    plt.xlabel("time (days since 1st Jan, 1980)")
    plt.ylabel("Temperature (Celsius)")
    plt.savefig("figures/raw_data_temp.png")
    plt.clf()

    ########################### FOURIER ANAYSIS #################################

    #detrend data
    rain_data_detrended = detrend(rain_data)
    temp_data_max_detrended = detrend(temp_data_max)
    temp_data_min_detrended = detrend(temp_data_min)

    #TODO: remove mean
    rain_data_detrended_centered = rain_data_detrended - np.mean(rain_data_detrended)
    temp_data_max_detrended_centered = temp_data_max_detrended - np.mean(temp_data_max_detrended)
    temp_data_min_detrended_centered = temp_data_min_detrended - np.mean(temp_data_min_detrended)

    #TODO: TRY DIFFERENT DT
    dt = 1
    freq_axis_rain, rain_data_ft = get_fft(rain_data_detrended_centered, dt)
    freq_axis_temp_max, temp_data_max_ft = get_fft(temp_data_max_detrended_centered, dt)
    freq_axis_temp_min, temp_data_min_ft = get_fft(temp_data_min_detrended_centered, dt)


    plt.plot(freq_axis_rain, rain_data_ft, label = "rainfall ft")
    plt.xlabel("Frequency (1/ (3650 day) )")
    plt.ylabel("FT(mm)")
    plt.savefig("figures/freq_spectrum_rainfall.png")
    plt.clf()

    plt.plot(freq_axis_temp_max, temp_data_max_ft, label ="temp_max")
    plt.xlabel("Frequency (1/ (3650 day) )")
    plt.ylabel("FT(degree celsius)")
    plt.savefig("figures/freq_spectrum_max_temp.png")
    plt.clf()

    plt.plot(freq_axis_temp_min, temp_data_max_ft, label = "temp_min")
    plt.xlabel("Frequency (1/ (3650 day) )")
    plt.ylabel("FT(degree celsius)")
    plt.savefig("figures/freq_spectrum_min_temp.png")
    plt.clf()

    ########################### POWER SPECTRUM ANAYSIS #################################

    rainfall_power_spectrum = np.abs(rain_data_ft) ** 2
    plt.plot(freq_axis_rain, rainfall_power_spectrum, label="rainfall power spectrum")
    plt.xlabel("Frequency (1/ (3650 day) )")
    plt.ylabel("|FT(mm)|^2")
    plt.savefig("figures/power_spectrum_rainfall.png")
    plt.clf()

    temp_max_power_spectrum = np.abs(temp_data_max_ft) ** 2
    plt.plot(freq_axis_rain, temp_max_power_spectrum, label="temp max power spectrum")
    plt.xlabel("Frequency (1/ (3650 day) )")
    plt.ylabel("|FT(celsius)|^2")
    plt.savefig("figures/power_spectrum_max_temp.png")
    plt.clf()

    temp_min_power_spectrum = np.abs(temp_data_min_ft) ** 2
    plt.plot(freq_axis_rain, temp_min_power_spectrum, label="temp min power spectrum")
    plt.xlabel("Frequency (1/ (3650 day) )")
    plt.ylabel("|FT(celsius)|^2")
    plt.savefig("figures/power_spectrum_min_temp.png")
    plt.clf()

    ##### get frequencies of peaks ######
    print(freq_axis_rain[np.where(rainfall_power_spectrum == max(rainfall_power_spectrum))])

    ########################### AUTO CORRELATION ANAYSIS #################################

    # auto_correlation_ft = temp_data_max_ft.conjugate() * temp_data_max_ft
    # time_axis, auto_correlation = get_ifft(auto_correlation_ft, dt)
    # plt.plot(time_axis, auto_correlation, label="auto_correlation_temp_max")

    auto_correlation = correlate(temp_data_max_detrended_centered, temp_data_max_detrended_centered)
    time_axis = np.arange(-len(auto_correlation)/2, len(auto_correlation)/2, 1)

    plt.plot(time_axis, auto_correlation, label="auto_correlation_temp_max")
    plt.xlabel("time lag (days)")
    plt.ylabel("mm ^ 2")
    plt.savefig("figures/auto_correlation_max_temp.png")
    plt.clf()

    auto_correlation = correlate(temp_data_min_detrended_centered, temp_data_min_detrended_centered)
    time_axis = np.arange(-len(auto_correlation) / 2, len(auto_correlation) / 2, 1)
    plt.plot(time_axis, auto_correlation, label="auto_correlation_temp_min")
    plt.xlabel("time lag (days)")
    plt.ylabel("celsius^2")
    plt.savefig("figures/auto_correlation_min_temp.png")
    plt.clf()

    auto_correlation = correlate(rain_data_detrended_centered, rain_data_detrended_centered)
    time_axis = np.arange(-len(auto_correlation) / 2, len(auto_correlation) / 2, 1)
    plt.plot(time_axis, auto_correlation, label="auto_correlation_rainfall")
    plt.xlabel("time lag (days)")
    plt.ylabel("celsius^2")
    plt.savefig("figures/auto_correlation_rainfall.png")
    plt.clf()


    ########################### CROSS CORRELATION ANAYSIS #################################
    # cross_correlation_ft = temp_data_max_ft.conjugate() * rain_data_ft
    # time_axis, cross_correlation = get_ifft(cross_correlation_ft,dt)
    #plt.plot(time_axis, cross_correlation, label = "cross_correlation_max_temp_rainfall")

    cross_correlation = correlate(temp_data_max_detrended_centered, rain_data_detrended_centered)
    time_axis = np.arange(-len(cross_correlation) / 2, len(cross_correlation) / 2, 1)
    plt.plot(time_axis, cross_correlation, label = "cross_correlation_max_temp_rainfall")
    plt.xlabel("time lag (days)")
    plt.ylabel("mm * celsius")
    plt.savefig("figures/cross_correlation_max_temp_rainfall.png")
    plt.clf()

    cross_correlation = correlate(temp_data_min_detrended_centered, rain_data_detrended_centered)
    time_axis = np.arange(-len(cross_correlation) / 2, len(cross_correlation) / 2, 1)
    plt.plot(time_axis,cross_correlation, label = "cross_correlation_min_temp_rainfall")
    plt.xlabel("time lag (days)")
    plt.ylabel("mm * celsius")
    plt.savefig("figures/cross_correlation_min_temp_rainfall.png")
    plt.clf()

    cross_correlation = correlate(temp_data_max_detrended_centered, temp_data_min_detrended)
    time_axis = np.arange(-len(cross_correlation) / 2, len(cross_correlation) / 2, 1)
    plt.plot(time_axis,cross_correlation, label = "cross_correlation_max_temp_min_temp")
    plt.xlabel("time lag (days)")
    plt.ylabel("celsius^2")
    plt.savefig("figures/cross_correlation_max_temp_min_temp.png")
    plt.clf()

