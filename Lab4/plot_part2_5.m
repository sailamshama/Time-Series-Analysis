load('freq_axis.txt')
load('power_spectrum_windowed.txt')
load('power_spectrum_unwindowed.txt')
plot(freq_axis, power_spectrum_windowed)
hold on
plot(freq_axis, power_spectrum_unwindowed)
legend("windowed", "unwindowed")