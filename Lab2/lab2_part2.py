
import numpy as np
import matplotlib.pyplot as plt
import cmath


def plot_figure(x, y_s,figsize, xlabel, ylabel, title, legend_labels,xlim,ylim,filename, markers):
    fig = plt.figure(figsize = figsize)
    plt.grid(1)
    for i in range(len(y_s)):
        plt.plot(x, y_s[i], markers[i], label = legend_labels[i])
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.legend()
    fig.savefig(filename)
    return

def boxcar(t, T):
    start = np.where(t >= 0)[0][0]
    end = np.where(t >= T)[0][0]
    boxcar = np.zeros(t.shape)
    boxcar[start:end] = np.ones(boxcar[start:end].shape)

    return boxcar

def hann(t,T):
    start = np.where(t >= 0 ) [0][0]
    end = np.where(t >= T)[0][0]
    hann = np.zeros(t.shape)
    hann[start:end] = 0.5 * ( np.ones(t[start:end].shape) - np.cos( ( 2 * cmath.pi / T ) * t[start:end] )  )
    return hann

def part2_1():

    #testing boxcar
    T = 10
    dt = 0.01
    t = np.arange(-4, T+10, dt)

    boxcar_fn = boxcar(t, T)
    hann_fn = hann(t, T)
    legend_labels = np.array(["boxcar", "hann"])
    markers = np.array(['-','-'])
    plot_figure(t, np.array([boxcar_fn, hann_fn]), (6,5), "time, t", "windowing function", "Boxcar and Hann function", legend_labels,(-4,T+10),(-0.5,2),"figures/part2_1.png", markers)


#TODO: fix this part, must be smoother?
def part2_2():
    T = 10
    dt = 0.01
    #TODO: see whether to change this range. I know it must be > 10
    t = np.arange(-3, T+1 , dt)

    boxcar_fn = boxcar(t, T)
    hann_fn = hann(t, T)

    boxcar_fft = np.fft.fftshift(np.fft.fft(boxcar_fn)) * dt
    hann_fft = np.fft.fftshift(np.fft.fft(hann_fn)) * dt

    boxcar_axis = np.fft.fftshift(np.fft.fftfreq(len(boxcar_fn) , dt))
    hann_axis = np.fft.fftshift(np.fft.fftfreq(len(hann_fn) , dt))

    fig = plt.figure(figsize = (10,10))
    plt.plot(boxcar_axis, np.absolute(boxcar_fft), label = "boxcar fft")
    plt.plot(hann_axis, np.absolute(hann_fft), label = "hann fft")
    plt.xlabel("frequency, f")
    plt.ylabel("amplitude of fft's")
    plt.title("fft's of boxcar and hann windowing functions")
    plt.legend()
    plt.xlim((-3,3))
    fig.savefig("figures/part2_2.png")

#test functions
def test_boxcar():
    # testing boxcar
    T = 10
    dt = 0.01
    t = np.arange(-4, T + 10, dt)

    boxcar_fn = boxcar(t, T)
    plt.plot(t, boxcar_fn)
    plt.show()


def test_hann():
    T = 10
    dt = 0.01
    t = np.arange(-4, T + 10, dt)

    hann_fn = hann(t, T)
    plt.plot(t, hann_fn)
    plt.show()

if __name__ == "__main__":
    part2_2()

