import numpy as np
import matplotlib.pyplot as plt
import cmath

def part1_2():

    f0 = 1
    fs = 12
    M = 1.056
    EPS = 0.053783

    dt = 1/12
    # TODO: redefine f vector
    f = np.arange(-5,5,dt)
    w = 2 * cmath.pi * f
    z = np.exp(-1j * w * dt)
    q = cmath.exp(-1j * 2 * cmath.pi * f0 / fs)
    p  = (1 + EPS) * q

    W = M * ( (z - q) * (z - q.conjugate()) )  / ( (z - p) * (z - p.conjugate()) )
    amplitude = W * np.conjugate(W)

    fig = plt.figure(figsize=(7,7))
    plt.grid(1)
    plt.plot(f, amplitude)
    plt.xlabel("frequency, f (1/year)")
    plt.ylabel("|W(f)|^2")
    plt.title("Power spectrum of notch filter")
    plt.savefig("figures/part1_2.png")

def part1_2_for_part2(f_axis):

    f0 = 1
    fs = 12
    M = 1.056
    EPS = 0.053783

    dt = 1/12
    # TODO: redefine f vector
    f = f_axis
    w = 2 * cmath.pi * f
    z = np.exp(-1j * w * dt)
    q = cmath.exp(-1j * 2 * cmath.pi * f0 / fs)
    p  = (1 + EPS) * q

    W = M * ( (z - q) * (z - q.conjugate()) )  / ( (z - p) * (z - p.conjugate()) )
    amplitude = W * np.conjugate(W)

    return(amplitude)
######################################################################################################################################
######################################################################################################################################

part1_2()