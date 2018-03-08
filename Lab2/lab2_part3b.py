import numpy as np
import matplotlib.pyplot as plt
import cmath

from lab2_part3a import RDFcacl, get_Pk, get_Sk, get_k


def get_k_w_kmax(kmax):
    dk = 0.12
    N = kmax/ dk
    k = np.arange(-N * dk, (N+dk) * dk, dk)
    return k

def extract_S(S, k):
    max_indx = int(k.shape[0] / 2)
    zero_index = int(S.shape[0] / 2)
    S = S[zero_index - max_indx + 1 : zero_index + max_indx + 1]
    return S

def RDFcacl_p6(S, dk, rho, factor):
    k = get_k_w_kmax(15.24 / factor)
    S = extract_S(S, k)

    molWeight = 39.948  # grams/mol
    Navogadro = 6.0221367e23  # atoms/mol
    rho = rho / molWeight * Navogadro * 1e-24

    N = int(k.shape[0])
    r_n = np.arange(N) * 2 * np.pi / (N * dk)  # should range from 0 to 52.15454471, array length N=255

    gn = 1 + ( 1 / (2 * cmath.pi**2 * rho * r_n) ) * np.fft.ifft( np.fft.fftshift(get_Pk(S, k)) ) * N * dk / (2 * cmath.pi)
    return gn, r_n

def part3_6():

    dk = 0.12
    rho = 1.4273
    S = get_Sk()


    gn, rn = RDFcacl(S, dk, rho)
    gn2, rn2 = RDFcacl_p6(S, dk, rho, 2)
    gn3, rn3 = RDFcacl_p6(S, dk, rho, 4)
    loc_25 = np.where(rn > 25)[0][0]
    loc_25_2 = np.where(rn2 > 25)[0][0]
    loc_25_3 = np.where(rn3 > 25)[0][0]

    fig = plt.figure(figsize=(7, 7))
    # TODO: figure out which part to plot imag or real
    plt.plot(rn[0:loc_25], np.absolute(gn[0:loc_25]), label = "kmax = 15.24")
    plt.plot(rn2[0:loc_25_2], np.absolute(gn2[0:loc_25_2]), marker='x', linestyle = "None", label = "kmax = 7.62")
    plt.plot(rn3[0:loc_25_3], np.absolute(gn3[0:loc_25_3]), marker='o', linestyle="None", label="kmax = 3.81")
    plt.title("Radial Distribution Function g(r) for argon data")
    plt.xlabel("r (armstrongs)")
    plt.ylabel("g(r)")
    plt.legend()
    fig.savefig("figures/part3_6.png")


##################################### TEST #####################################

def test_k_w_kmax():
    k = get_k_w_kmax(7.62)
    print(np.where(k<0)[0].shape, np.where(k>=0)[0].shape, k.shape)
    print(k[-1])

def test_RDFcacl_p6(S, dk, rho):
    return

def test_S_extraction():
    # S = get_Sk()
    # k = get_k_w_kmax(15.24/2)
    # #get zero
    # zero_index = int(S.shape[0]/2)
    # k1 = get_k()
    # print(k1[zero_index], S[zero_index])
    # max_indx = int(k.shape[0] /2)
    # print(k.shape, max_indx)
    # new_S = S[zero_index - max_indx + 1 : zero_index + max_indx + 1]
    # print(new_S.shape)
    k = get_k_w_kmax(15.24/2)
    S = get_Sk()
    S = extract_S(S, k)
    print(S.shape)

def test_RDFcalc_p6():
    S = get_Sk()
    dk = 0.12
    rho = 1.4273
    factor = 2
    gn, rn = RDFcacl_p6(S, dk, rho, factor)
    plt.plot(rn, gn)
    plt.show()

def test_factor():
    S = get_Sk()
    dk = 0.12
    rho = 1.4273
    factor = 4
    gn, rn = RDFcacl_p6(S, dk, rho, factor)
    plt.plot(rn, gn)
    plt.show()

##################################### MAIN #####################################

if __name__ == "__main__":
    part3_6()