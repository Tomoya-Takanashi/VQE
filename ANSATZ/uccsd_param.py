import numpy
import sys
sys.path.append('/users/tomoya12/work/M1/VQE')
from openfermionpyscf import prepare_pyscf_molecule
from . import amp_comb
from HAMILTONIAN import hamiltonian

def param_sort(target):
    if target[1][0] > target[1][1]:
        temp = target[1][0]
        target[1][0] = target[1][1]
        target[1][1] = temp
        tag = tag * -1

    occ_f = target[0][0] % 2
    occ_s = target[0][1] % 2
    vir_f = target[1][0] % 2
    vir_s = target[1][1] % 2

    if occ_f != vir_f and occ_s != vir_s:
        temp = target[0][0]
        target[0][0] = target[0][1]
        target[0][1] = temp
    return target

def zero_uccsd_init(e_num,qubit_num):
    single_comb,double_comb = amp_comb.combination(e_num,qubit_num)
    param_len = len(single_comb) + len(double_comb)
    phi = numpy.zeros(param_len)
    return phi

def random_uccsd_init(e_num,qubit_num):
    single_comb,double_comb = amp_comb.combination(e_num,qubit_num)
    param_len = len(single_comb) + len(double_comb)
    phi = numpy.random.random(param_len)*1e-1
    return phi

def mp2_uccsd_init(ele_num,qubit_num,mo_a,mo_b,hcore,aoint2,ds=True):
    ampComb = []
    single_comb,double_comb = amp_comb.combination(ele_num,qubit_num)
    #print(single_comb,double_comb)
    double = []
    for i in range(len(double_comb)):
        double.append(param_sort(double_comb[i]))

    one_body_electron,two_body_electron= hamiltonian.compute_spin_integral(mo_a,mo_b,hcore,aoint2,mp2=1)
    orbital_energy = []
    for i in range(len(one_body_electron)):
        orbital_energy.append(one_body_electron[i][i])

    for i in range(len(double_comb)):
        ampComb.append(double_comb[i][0] + double_comb[i][1])

    #Initialize parameter #
    if ds :
        amp_init = []
        for i in range(len(ampComb)):
            J = two_body_electron[ampComb[i][0]][ampComb[i][3]][ampComb[i][1]][ampComb[i][2]]
            K = two_body_electron[ampComb[i][0]][ampComb[i][2]][ampComb[i][1]][ampComb[i][3]]
            E = orbital_energy[ampComb[i][0]] + orbital_energy[ampComb[i][1]] - orbital_energy[ampComb[i][2]] - orbital_energy[ampComb[i][3]]
            amp_init.append((J - K)/E/8)

        for i in range(len(single_comb)):
            amp_init.append(0)
    else:
        amp_init = []

        for i in range(len(single_comb)):
            amp_init.append(0)

        for i in range(len(ampComb)):
            J = two_body_electron[ampComb[i][0]][ampComb[i][3]][ampComb[i][1]][ampComb[i][2]]
            K = two_body_electron[ampComb[i][0]][ampComb[i][2]][ampComb[i][1]][ampComb[i][3]]
            E = orbital_energy[ampComb[i][0]] + orbital_energy[ampComb[i][1]] - orbital_energy[ampComb[i][2]] - orbital_energy[ampComb[i][3]]
            amp_init.append((J - K)/E/8)

    return amp_init
