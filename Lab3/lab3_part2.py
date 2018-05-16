import numpy as np
import matplotlib.pyplot as plt
from lab3_part1 import part1_2_for_part2

def part2_1():
    f0 = 1
    fs = 12
    M = 1.056
    EPS = 0.053783

    a=M/((1+EPS)**2)
    b=-2*M*np.cos(2*np.pi*f0/fs)/((1+EPS)**2)
    c=M/((1+EPS)**2)
    B=-2*(1+EPS)*np.cos(2*np.pi*f0/fs)/((1+EPS)**2)
    C=1/((1+EPS)**2)

    return a, b, c, B, C

def ratFilter(N,D,x):
    dt = 1.0
    #TODO: multiply g by dt?

    #multiply polynomials x and N (but reverse to be consistent with polynum convention)
    g = np.polymul(x[::-1], N[::-1])
    g = g[::-1] #reverse back

    y = np.arange(0, len(x), dtype=float)

    #build W through recursive filtering
    y[0] = g[0]/D[0]*(1.0/dt)
    y[1] = (g[1]/dt - D[1]*y[0]) / D[0]
    y[2] = (g[2]/dt - D[1]*y[1] - D[2]*y[0]) / D[0]

    #pad g with zeros
    g = np.append(g, np.zeros(abs(len(y) - len(g))))

    for i in range(3, len(y)):
        y[i] = (g[i]/dt - D[1]*y[i-1] - D[2]*y[i-2]) / D[0]

    return y

def part2_3():
    a, b, c, B, C = part2_1()
    N = np.array([a, b,c])
    D = np.array([1, B, C])
    dt = 1/12

    t = np.arange(0,50,dt)
    x  = np.zeros((t.shape))
    x[0] = 1/dt

    y = ratFilter(N,D,x)

    return y

def part2_4():
    dt = 1 / 12
    t = np.arange(0, 50, dt)
    y = part2_3()

    #TODO: label these plots properly
    fig = plt.figure()
    plt.plot(t, y)
    plt.xlabel("time, t")
    plt.ylabel("y")
    fig.savefig("figures/part2_4a.png")

    fig = plt.figure()
    y_fft = np.fft.fft(y) * dt
    y_fft_shifted = np.fft.fftshift(y_fft)
    f_axis = np.fft.fftshift(np.fft.fftfreq(len(y), dt))
    W = np.sqrt(part1_2_for_part2(f_axis))
    plt.plot(f_axis, abs(y_fft_shifted), label = "FFT ")
    plt.plot(f_axis, W, marker = "x", linestyle="None", label = "theoretical")
    plt.xlabel("frequency, f (1/year)")
    plt.ylabel("Frequency response of notch filter")
    plt.legend()
    fig.savefig("figures/part2_4b.png")


######################################################################################################################################
######################################################################################################################################
part2_4()