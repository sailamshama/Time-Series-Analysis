#Just gonna use CSC411_py3 environment#

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy

####### PART 1 ########

def gaussian(t, center, tH):
    gaussian = ( 1 / ( np.sqrt(np.pi) * tH) ) * np.exp( - ( (t-center) /tH) **2 )

    return gaussian

def gaussian_FT_analytical(w, tH):
    gaussian_ft = exp ( - (w * tH / 2.0 ) ** 2  )
    return gaussian_ft

def plot_figure(x, y_s,figsize, xlabel, ylabel, title, legend_labels,xlim,ylim,filename):
    fig = plt.figure(figsize = figsize)
    plt.grid(1)
    for i in range(len(y_s)):
        plt.plot(x, y_s[i], label = legend_labels[i])
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.legend()
    fig.savefig(filename)
    return

def part1_1():

    tH = np.array([10,20])
    center = 0
    t = np.arange(0, 801, 1)

    gaussians = np.zeros((tH.shape[0], t.shape[0]))
    gaussian_labels = np.chararray(tH.shape)
    for i in range(len(tH)):
        gaussians[i, :] = gaussian(t, center, tH[i])
        gaussian_labels[i] = "tH = %s" % tH[i]

    plot_figure(t, gaussians, (10,15), "t", "Gaussian distributions", "Gaussian distribution functions with half durations 10 and 20", gaussian_labels, (0,800), (-0.1,0.1), "figures/part1_guassians.png")

def part1_3():
    tH = np.array([10, 20])
    t = np.arange(0, 801, 1)
    #TODO: figure out how to set w i.e. in what range to plot G(w)
    w = np.arange(0,800,1)

    gaussians = np.zeros((t.shape[0], tH.shape[0]))
    gaussians_fft = np.zeros(gaussians.shape)

    #TODO: figure out how to use the shift property?


    # gaussians_analytical = np.zeros(gaussians.shape)
    # gaussian_labels = np.chararray(tH.shape)
    # gaussian_labels_fft = np.chararray(tH.shape)
    # gaussian_labels_analytical = np.chararray(tH.shape)
    #

    #print(gaussians.shape, gaussians_fft.shape, gaussians_analytical.shape)
    #
    # for i in range(len(tH)):
    #     gaussians[i, :] = gaussian(t, center, tH[i])
    #     gaussians_fft = np.fft.fft(gaussians[i])
    #     gaussians_analytical = gaussian_FT_analytical(gaussians[i])
    #     gaussian_labels[i] = "g(t) with tH = %s" % tH[i]
    #     gaussian_labels_analytical = "analytical G(w) with tH = %s" %tH[i]
    #     gaussians_labels_fft = "FFT of g(t) with tH = %s" % tH[i]
    #


################## PART 2 ##################################

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
    hann[start:end] = 0.5 * ( np.ones(t[start:end].shape) - np.cos( ( 2 * math.pi / T ) * t[start:end] )  )
    return hann

def test_boxcar():
    #testing boxcar
    T = 10
    dt = 0.01
    t = np.arange(-4, T+10, dt)
                                     
    boxcar_fn = boxcar(t, T)
    plt.plot(t, boxcar_fn)
    plt.show()

def test_hann():
    T = 10
    dt = 0.01
    t = np.arange(-4, T+10, dt)

    hann_fn = hann(t, T)
    plt.plot(t, hann_fn)
    plt.show()                      

def part2_1():

    #testing boxcar
    T = 10
    dt = 0.01
    t = np.arange(-4, T+10, dt)

    boxcar_fn = boxcar(t, T)
    hann_fn = hann(t, T)
    legend_labels = np.array(["boxcar", "hann"])
    plot_figure(t, np.array([boxcar_fn, hann_fn]), (6,5), "time, t", "windowing function", "Boxcar and Hann function", legend_labels,(-4,T+10),(-0.5,2),"figures/part2_1.png")




################### MAIN ######################
#part1_1()
#part1_3()
part2_1()