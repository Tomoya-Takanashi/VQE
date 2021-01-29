from qiskit import Aer,aqua,execute,QuantumCircuit
from math import pi
from . import amp_comb

def param_sort(target):
    tag = 1
    if target[0][0] > target[0][1]:
        temp = target[0][0]
        target[0][0] = target[0][1]
        target[0][1] = temp
        tag = tag * -1

    if target[1][0] > target[1][1]:
        temp = target[1][0]
        target[1][0] = target[1][1]
        target[1][1] = temp
        tag = tag * -1

    return target,tag

def EO_judge(target):
    occ_f = target[0][0] % 2
    occ_s = target[0][1] % 2
    vir_f = target[1][0] % 2
    vir_s = target[1][1] % 2
    tag = 1

    if occ_f == vir_f and occ_s == vir_s:
        tag = 1
    elif occ_f != vir_f and occ_s != vir_s:
        tag = -1

    return tag

def exp_single(r3,r1,uccsd,single_comb,Gate_s,amp_s,amp_d,phi,trotter):

    if(Gate_s[r3][0] == 'H'):
        uccsd.h(single_comb[r1][0])
    elif(Gate_s[r3][0] == 'Y'):
        uccsd.rx(-pi/2,single_comb[r1][1])
         
    if(Gate_s[r3][1] == 'H'):
        uccsd.h(single_comb[r1][0])
    elif(Gate_s[r3][1] == 'Y'):
        uccsd.rx(-pi/2,single_comb[r1][1])
        
    loop_s = single_comb[r1][0]
    loop_num = single_comb[r1][1] - single_comb[r1][0]
    
    for r4 in range(loop_num):
        uccsd.cx(loop_s,loop_s+1)
        loop_s = loop_s + 1

    if(r3 == 0):
        uccsd.rz(phi[amp_s+amp_d]*(-1)/2/trotter,single_comb[r1][1])
    elif(r3 == 1):
        uccsd.rz(phi[amp_s+amp_d]/2/trotter,single_comb[r1][1])

    loop_s = single_comb[r1][1]
            
    for r4 in range(loop_num):
        uccsd.cx(loop_s-1,loop_s)
        loop_s = loop_s - 1

    if(Gate_s[r3][0] == 'H'):
        uccsd.h(single_comb[r1][0])
    elif(Gate_s[r3][0] == 'Y'):
        uccsd.rx(pi/2,single_comb[r1][1])
    
    if(Gate_s[r3][1] == 'H'):
        uccsd.h(single_comb[r1][0])
    elif(Gate_s[r3][1] == 'Y'):
        uccsd.rx(pi/2,single_comb[r1][1])

    return uccsd

def exp_double(r3,r1,uccsd,double_comb,Gate,amp_s,amp_d,phi,tag,trotter):

    if(Gate[r3][0] == 'H'):
        uccsd.h(double_comb[r1][0][0])
    elif(Gate[r3][0] == 'Y'):
        uccsd.rx(-pi/2,double_comb[r1][0][0])
        
    if(Gate[r3][1] == 'H'):
        uccsd.h(double_comb[r1][0][1])
    elif(Gate[r3][1] == 'Y'):
        uccsd.rx(-pi/2,double_comb[r1][0][1])

    if(Gate[r3][2] == 'H'):
        uccsd.h(double_comb[r1][1][0])
    elif(Gate[r3][2] == 'Y'):
        uccsd.rx(-pi/2,double_comb[r1][1][0])

    if(Gate[r3][3] == 'H'):
        uccsd.h(double_comb[r1][1][1])
    elif(Gate[r3][3] == 'Y'):
        uccsd.rx(-pi/2,double_comb[r1][1][1])
    occ_cx_loop = double_comb[r1][0][1] - double_comb[r1][0][0]
    vir_cx_loop = double_comb[r1][1][1] - double_comb[r1][1][0]
        
    occ_loop = double_comb[r1][0][0]
    vir_loop = double_comb[r1][1][0]
    for r5 in range(occ_cx_loop):
        uccsd.cx(occ_loop,occ_loop+1)
        occ_loop = occ_loop + 1

    uccsd.cx(double_comb[r1][0][1],double_comb[r1][1][0])

    for r5 in range(vir_cx_loop):
        uccsd.cx(vir_loop,vir_loop+1)
        vir_loop = vir_loop + 1

    if(r3 == 0 or r3 == 1 or r3 == 2 or r3 == 3):
        uccsd.rz(phi[amp_s+amp_d]*(-1)*tag/8/trotter,double_comb[r1][1][1])
    elif(r3 == 4 or r3 == 5 or r3 == 6 or r3 == 7):
        uccsd.rz(phi[amp_s+amp_d]*tag/8/trotter,double_comb[r1][1][1])

    occ_loop = double_comb[r1][0][1]
    vir_loop = double_comb[r1][1][1]
            
    for r5 in range(vir_cx_loop):
        uccsd.cx(vir_loop-1,vir_loop)
        vir_loop = vir_loop -1

    uccsd.cx(double_comb[r1][0][1],double_comb[r1][1][0])
    
    for r5 in range(occ_cx_loop):
        uccsd.cx(occ_loop-1,occ_loop)
        occ_loop = occ_loop - 1
    if(Gate[r3][0] == 'H'):
        uccsd.h(double_comb[r1][0][0])
    elif(Gate[r3][0] == 'Y'):
        uccsd.rx(pi/2,double_comb[r1][0][0])

    if(Gate[r3][1] == 'H'):
        uccsd.h(double_comb[r1][0][1])
    elif(Gate[r3][1] == 'Y'):
        uccsd.rx(pi/2,double_comb[r1][0][1])

    if(Gate[r3][2] == 'H'):
        uccsd.h(double_comb[r1][1][0])
    elif(Gate[r3][2] == 'Y'):
        uccsd.rx(pi/2,double_comb[r1][1][0])

    if(Gate[r3][3] == 'H'):
        uccsd.h(double_comb[r1][1][1])
    elif(Gate[r3][3] == 'Y'):
        uccsd.rx(pi/2,double_comb[r1][1][1])

    return uccsd

def uccsd_circuit(ele_num,qubit_num,state,phi,trotter=None,ds=True):
    
    if trotter is None:
        trotter = 1

    machine = Aer.get_backend('statevector_simulator')
    
    #Gate_s = [['Y','H'],['H','Y']]
    #Gate = [['H','H','Y','H'],['Y','H','Y','Y'],['H','Y','Y','Y'],['H','H','H','Y'],['Y','H','H','H'],['H','Y','H','H'],['Y','Y','Y','H'],['Y','Y','H','Y']]

    Gate_s = [['H','Y'],['Y','H']]
    Gate = [['Y','Y','H','Y'],['Y','Y','Y','H'],['H','Y','H','H'],['Y','H','H','H'],
            ['H','H','H','Y'],['H','Y','Y','Y'],['Y','H','Y','Y'],['H','H','Y','H']]

    single_comb,double_comb = amp_comb.combination(ele_num,qubit_num)
    #single_comb,double_comb = combination(ele_num,qubit_num)

    uccsd = aqua.circuits.StateVectorCircuit(state)
    uccsd = uccsd.construct_circuit()

    #uccsd = state
    amp_s = 0
    amp_d = 0
    
    if ds:
        for t in range(trotter):
            amp_s = 0
            amp_d = 0
            for r1 in range(len(double_comb)):
                new,tag = param_sort(double_comb[r1])
                double_comb[r1] = new
                for r3 in range(8):
                    uccsd = exp_double(r3,r1,uccsd,double_comb,Gate,amp_s,amp_d,phi,tag,trotter)
                amp_d = amp_d + 1
    
            for r1 in range(len(single_comb)):
                for r3 in range(2):
                    uccsd = exp_single(r3,r1,uccsd,single_comb,Gate,amp_s,amp_d,phi,trotter)
                amp_s = amp_s + 1

    else:
        for t in range(trotter):
            amp_s = 0
            amp_d = 0
            for r1 in range(len(single_comb)):
                for r3 in range(2):
                    uccsd = exp_single(r3,r1,uccsd,single_comb,Gate,amp_s,amp_d,phi,trotter)
                amp_s = amp_s + 1

            for r1 in range(len(double_comb)):
                new,tag = param_sort(double_comb[r1])
                double_comb[r1] = new
                for r3 in range(8):
                    uccsd = exp_double(r3,r1,uccsd,double_comb,Gate,amp_s,amp_d,phi,tag,trotter)
                amp_d = amp_d + 1

    result = execute(uccsd,machine).result()
    statevector = result.get_statevector(uccsd)
    uccsd = 0
    return statevector
    #return uccsd
    
    
