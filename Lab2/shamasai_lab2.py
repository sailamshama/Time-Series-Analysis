#Just gonna use CSC411_py3 environment#

import numpy as np
import matplotlib.pyplot as plt
import cmath

###################### PART 1 #########################

def gaussian(t, center, tH):
    gaussian = ( 1 / ( np.sqrt(np.pi) * tH) ) * np.exp( - ( (t-center) /tH) **2 )
    return gaussian

def gaussian_FT_analytical(f, tH):
    gaussian_ft = np.exp ( - (f * tH / 2.0 ) ** 2  )
    return gaussian_ft

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


def part1_1():

    tH = np.array([10,20])
    center = 0
    t = np.arange(-100,100, 1)
    #TODO: fix this convention. you want data in the columns, not rows :(
    gaussians = np.zeros((tH.shape[0], t.shape[0]))
    gaussian_labels = np.array(["tH = 10", "tH = 20 "])
    markers = np.array(['-', '-'])
    for i in range(len(tH)):
        gaussians[i, :] = gaussian(t, center, tH[i])
    plot_figure(t, gaussians, (7,7), "t", "Gaussian distributions", "Gaussian distribution functions with half durations 10 and 20", gaussian_labels, (-100,100), (-0.01,0.06), "figures/part1_guassians.png", markers)

def part1_3():
    th = np.array([10,20])
    T = 100
    dt = 1
    t = np.arange(-T, T, dt)
    center = 0
    #create gaussian, gaussian_fft, gaussian_FT_analytical, amplitudes for each of them and their axes

    #create (4, len(f_axis) np array of gaussian transforms
    #first computed then analytical

    for i in range(len(th)):
        gaussian_fn = gaussian(t, center, th)
        gaussian_fn_fft = np.fft.fft(gaussian_fn)
        gaussian_fn_fft_shifted = np.fft.fftshift(gaussian_fn_fft)

        #TODO: figure out how this thing works
        f_axis = np.fft.fftfreq(len(gaussian_fn), dt)
        f_axis_shifted = np.fft.fftshift(np.fft.fftfreq(len(gaussian_fn) , dt))

        #TODO: figure out what expnonential to multiply by
        gaussian_fn_analytical = gaussian_FT_analytical(2*cmath.pi*f_axis_shifted, th)

        gaussian_fn_fft_shifted_amplitude = np.absolute(gaussian_fn_fft_shifted)
        gaussian_fn_analytical_amplitude = np.absolute(gaussian_fn_analytical)
        gaussian_fn_amplitude = np.absolute(gaussian_fn)

        if i == 0 :
            gaussian_transforms = np.vstack((gaussian_fn_analytical_amplitude, gaussian_fn_fft_shifted_amplitude))
        else:
            gaussian_transforms = np.vstack((gaussian_transforms, gaussian_fn_analytical_amplitude, gaussian_fn_fft_shifted_amplitude))

    labels = np.array(["analytical, tH = 10", "computed, tH = 10", "analytical, tH = 20", "computed, tH = 20"])
    markers = np.array(["-", "x", "-", "x"])

    #TODO: figure out if frequency axes are same for th = 10 and th = 20
    plot_figure()
def part1_4():
    t = np.arange(-5,5,0.01)
    th = 10
    center = 0
    general_fn = np.sin(t)
    gaussian_fn = gaussian(t, center, th)
    filtered = np.convolve(general_fn, gaussian_fn)
    plt.plot(general_fn, label = "general function - sine")
    plt.plot(gaussian_fn, label= "gaussian with th = 10")
    plt.plot(filtered, label = "convolved sine and gaussian")
    plt.legend()
    print(filtered.shape, gaussian_fn.shape, general_fn.shape)
    plt.show()


#test functions

def test_fft():
    th = 10
    T = 100
    dt = 1
    t = np.arange(-T, T, dt)
    center = 0
    #create gaussian, gaussian_fft, gaussian_FT_analytical, amplitudes for each of them and their axes


    gaussian_fn = gaussian(t, center, th)
    gaussian_fn_fft = np.fft.fft(gaussian_fn)
    gaussian_fn_fft_shifted = np.fft.fftshift(gaussian_fn_fft)

    #TODO: figure out how this thing works
    f_axis = np.fft.fftfreq(len(gaussian_fn), dt)
    f_axis_shifted = np.fft.fftshift(np.fft.fftfreq(len(gaussian_fn) , dt))

    #TODO: figure out what expnonential to multiply by
    gaussian_fn_analytical = gaussian_FT_analytical(2*cmath.pi*f_axis_shifted, th)


    gaussian_fn_fft_shifted_amplitude = np.absolute(gaussian_fn_fft_shifted)
    gaussian_fn_analytical_amplitude = np.absolute(gaussian_fn_analytical)
    gaussian_fn_amplitude = np.absolute(gaussian_fn)

    plt.plot(f_axis_shifted, gaussian_fn_fft_shifted_amplitude, label = "computed G(w), tH = %s" %th)
    plt.plot(f_axis_shifted, gaussian_fn_analytical_amplitude, label = "analytical G(w), tH = %s" %th)
    plt.legend()
    plt.show()


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
    plot_figure(t, np.array([boxcar_fn, hann_fn]), (6,5), "time, t", "windowing function", "Boxcar and Hann function", legend_labels,(-4,T+10),(-0.5,2),"figures/part2_1.png")

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



################## MAIN ##################

if __name__ == "__main__":
    part1_1()











