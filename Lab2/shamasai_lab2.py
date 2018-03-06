#Just gonna use CSC411_py3 environment#

import numpy as np
import matplotlib.pyplot as plt
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
    center = 400
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
    w = np.arange(0,800,1)

    gaussians = np.zeros((t.shape[0], tH.shape[0]))
    gaussians_fft = np.zeros(gaussians.shape)
    gaussians_analytical = np.zeros(gaussians.shape)
    gaussian_labels = np.chararray(tH.shape)
    gaussian_labels_fft = np.chararray(tH.shape)
    gaussian_labels_analytical = np.chararray(tH.shape)

    print(gaussians.shape, gaussians_fft.shape, gaussians_analytical.shape)

    # for i in range(len(tH)):
    #     gaussians[i, :] = gaussian(t, center, tH[i])
    #     gaussians_fft = np.fft.fft(gaussians[i])
    #     gaussians_analytical = gaussian_FT_analytical(gaussians[i])
    #     gaussian_labels[i] = "g(t) with tH = %s" % tH[i]
    #     gaussian_labels_analytical = "analytical G(w) with tH = %s" %tH[i]
    #     gaussians_labels_fft = "FFT of g(t) with tH = %s" % tH[i]


    
################### MAIN ######################

#part1_1()
part1_3()