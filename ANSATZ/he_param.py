import numpy as np

def zero_he_init(qubit_num,depth):
    return np.zeros(2*qubit_num*(depth+1))*1e-1

def random_he_init(qubit_num,depth):
    return np.random.random(2*qubit_num*(depth+1))*1e-1

def zero_he_type2_init(qubit_num,depth):
    return np.zeros(3*qubit_num*(depth+1))*1e-1

def random_he_type2_init(qubit_num,depth):
    return np.random.random(3*qubit_num*(depth+1))*1e-1
