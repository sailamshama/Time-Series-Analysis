import numpy as np
import matplotlib.pyplot as plt
import cmath

from lab2_part3a import RDFcacl, get_Pk, get_Sk, get_k

#TODO: include plt.grid(1) for all plots

def get_k_subsample(factor):
    return get_k()[::factor] # get every 2nd (or 4th or 6th) value of the array


def RDFcacl_p7(S, dk, rho):
    k = get_k_subsample(int(dk/0.12))
    N = int(k.shape[0])
    molWeight = 39.948  # grams/mol
    Navogadro = 6.0221367e23  # atoms/mol
    rho = rho / molWeight * Navogadro * 1e-24
    r_n = np.arange(N) * 2 * np.pi / (N * dk)  # should range from 0 to 52.15454471, array length N=255

    gn = 1 + ( 1 / (2 * cmath.pi**2 * rho * r_n) ) * np.fft.ifft( np.fft.fftshift(get_Pk(S, k)) ) * N * dk / (2 * cmath.pi)
    return gn, r_n


def part3_7():
    dk = 0.12
    rho = 1.4273

    S = get_Sk()
    S_subsampled2 = S[::2] # I wouldn't have done this if we didn't need an argument S in RDFcalc
    S_subsampled4 = S[::4]

    gn, rn = RDFcacl(S, dk, rho)
    gn2, rn2 = RDFcacl_p7(S_subsampled2, dk*2, rho)
    gn4, rn4 = RDFcacl_p7(S_subsampled4, dk*4, rho)

    fig = plt.figure(figsize=(7, 7))
    # TODO: figure out which part to plot imag or real

    plt.plot(rn, np.absolute(gn), label="dk = 0.12")
    plt.plot(rn2, np.absolute(gn2), linestyle=":", label="dk = 0.24")
    plt.plot(rn4, np.absolute(gn4), linestyle=":", label="dk = 0.48")

    plt.xlim((0,25))
    plt.title("Radial Distribution Function g(r) for argon data")
    plt.xlabel("r (armstrongs)")
    plt.ylabel("g(r)")
    plt.legend()
    fig.savefig("figures/part3_7.png")

##################################### TEST #####################################
def test_subsamples():
    factor = 2
    k = get_k_subsample(factor)
    S = get_Sk()[::factor]

    print(k.shape, k[-1])
    print(S.shape)

##################################### MAIN #####################################



if __name__ == "__main__":
    part3_7()