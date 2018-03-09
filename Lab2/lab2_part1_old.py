#Just gonna use CSC411_py3 environment#

import numpy as np
import matplotlib.pyplot as plt
import cmath

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
    plot_figure(t, gaussians, (7,7), "t", "Gaussian distributions", "Gaussian distribution functions with half durations 10 and 20", gaussian_labels, (-100,100), (-0.01,0.06), "figures/part1_1.png", markers)

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
    dt = 0.01
    t = np.arange(-250,250,dt)
    center = 0
    general_fn = np.sin(t)

    for th in [10,20]:
        gaussian_fn = gaussian(t, center, th)
        general_fn_fft = np.fft.fft(general_fn)
        general_fn_fft_axis =np.fft.fftfreq(len(general_fn), dt)

        filtered = np.convolve(general_fn, gaussian_fn)
        filtered_axis = np.linspace(-5, 5, filtered.shape[0])
        filtered_fft = np.fft.fft(filtered)
        filtered_fft_axis = np.fft.fftfreq(len(filtered), dt)


        fig = plt.figure(figsize=(15,15))

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
        plt.xlim((-0.5,0.5))

        c = fig.add_subplot(2,2,3)
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

        fig.suptitle("Filtering with gaussian with th = %s"%th)
        plt.savefig("figures/part1_4_tH"+str(th)+".png")

    #create the necessary plots
    # fig = plt.figure(figsize=(7,15))
    # for th in [10, 20]:
    #     gaussian_fn = gaussian(t, center, th)
    #     general_fn_fft = np.fft.fft(general_fn)
    #     general_fn_fft_axis = np.fft.fftfreq(len(general_fn), dt)
    #
    #     filtered = np.convolve(general_fn, gaussian_fn)
    #     filtered_axis = np.linspace(-5, 5, filtered.shape[0])
    #     filtered_fft = np.fft.fft(filtered)
    #     filtered_fft_axis = np.fft.fftfreq(len(filtered), dt)
    #
    #     fig.add_subplot(2, 2, 3)
    #     plt.plot(filtered_axis, filtered)
    #     plt.xlabel("t")
    #     plt.ylabel("f(t) * g(t)")
    #     plt.title("Convolution of f(t) with gaussian")
    # fig.savefig("figures/part1_4.png")

#    print(filtered.shape, general_fn.shape, general_fn_fft.shape, general_fn_fft_axis.shape, filtered_fft.shape, filtered_fft_axis.shape)
#    plt.show()


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


################## MAIN ##################

if __name__ == "__main__":
    part1_1()











