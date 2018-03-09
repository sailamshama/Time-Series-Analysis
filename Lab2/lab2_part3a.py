import numpy as np
import matplotlib.pyplot as plt
import cmath

def get_Sk():
    # the structure factor data for Argon
    YanData = np.array([
        0.161000, 0.138000, 0.120000, 0.105000, 0.088000, 0.086000, 0.087000, 0.088000,
        0.090000, 0.094000, 0.105000, 0.125000, 0.152000, 0.277000, 0.468000, 1.248000,
        2.391000, 2.604000, 1.725000, 1.160000, 0.904000, 0.761000, 0.640000, 0.626000,
        0.633000, 0.680000, 0.744000, 0.869000, 1.038000, 1.183000, 1.225000, 1.237000,
        1.231000, 1.184000, 1.099000, 1.011000, 0.932000, 0.870000, 0.801000, 0.771000,
        0.871000, 0.961000, 1.026000, 1.075000, 1.102000, 1.107000, 1.087000, 1.064000,
        1.033000, 0.998000, 0.968000, 0.945000, 0.936000, 0.946000, 0.966000, 0.987000,
        1.004000, 1.018000, 1.027000, 1.024000, 1.022000, 1.015000, 1.006000, 0.998000,
        0.997000, 0.996000, 0.996000, 0.995000, 0.995000, 0.995000, 0.995000, 0.996000,
        1.001000, 1.001000, 1.002000, 1.003000, 1.003000, 1.003000, 1.003000, 1.002000,
        1.001000, 1.000000, 1.001000, 1.000000, 0.999000, 0.998000, 0.997000, 0.998000,
        0.998000, 0.999000, 1.000000, 1.000000, 1.001000, 1.001000, 1.001000, 1.001000,
        1.001000, 1.001000, 1.001000, 1.000000, 1.000000, 0.999000, 0.999000, 0.999000,
        1.000000, 0.999000, 0.999000, 0.999000, 1.000000, 1.000000, 1.000000, 1.002000,
        1.001000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 0.999000,
        0.999000, 0.999000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000
    ])
    dk = 0.12  # inverse Angstroms
    massRho = 1.4273  # grams/cc
    molWeight = 39.948  # grams/mol
    Navogadro = 6.0221367e23  # atoms/mol

    return np.append( np.flip(YanData[1:], 0), YanData )

def get_k():
    dk = 0.12
    k = np.arange(-127 * dk, (127+dk) * dk, dk)
    return k

def get_Pk(S, k):
    return -cmath.pi * 1j * k * (S -1)

def RDFcacl(S, dk, rho):
    N = 255
    dk = 0.12
    molWeight = 39.948  # grams/mol
    Navogadro = 6.0221367e23  # atoms/mol
    rho = rho / molWeight * Navogadro * 1e-24     # you need to make sure the units on rho match
    # in the file, you get them as grams,/cc, but you plot it as number of atoms per A^3
    r_n = np.arange(N) * 2 * np.pi / (N * dk)  # should range from 0 to 52.15454471, array length N=255
    k = get_k()
    gn = 1 + ( 1 / (2 * cmath.pi**2 * rho * r_n) ) * np.fft.ifft( np.fft.fftshift(get_Pk(S, k)) ) * N * dk / (2 * cmath.pi)
    return gn, r_n

def part3_2():
    k = get_k()
    Pk = get_Pk(get_Sk(), k)
    fig = plt.figure()
    plt.plot(k, Pk.imag)
    plt.xlabel("inverse armstrong")
    plt.ylabel("Pk")
    plt.title("Plot of Pk")
    fig.savefig("figures/part3_2.png")

def part3_4():
    dk = 0.12
    rho = 1.4273
    gn, rn = RDFcacl(get_Sk(), dk, rho)

    loc_25 = np.where(rn>25)[0][0]
    fig = plt.figure(figsize=(7,7))
    #TODO: figure out which part to plot imag or real
    plt.plot(rn[0:loc_25], np.absolute(gn[0:loc_25]))
    plt.title("Radial Distribution Function g(r) for argon data")
    plt.xlabel("r (armstrongs)")
    plt.ylabel("g(r)")
    fig.savefig("figures/part3_4.png")

#TODO: figure out 3.5
# def part3_5():
#     dk = 0.12
#     rho = 1.4273
#     gn, rn = RDFcacl(get_Sk(), dk, rho)
#     dr = rn[1] - rn[0]
#     Ravg = ( np.sum(gn) * dr ) / (rn[-1] - rn[0])
#
#     return Ravg

##################################### TEST #####################################

def test_Sk():
    Sk = get_Sk()
    print(Sk.shape, Sk)

def test_k():
    k = get_k()
    print(k[-1], k.shape)

##################################### MAIN #####################################

if __name__ == "__main__":
    part3_2()
    part3_4()