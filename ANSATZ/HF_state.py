from qiskit import Aer,QuantumCircuit,execute
from math import pi,sqrt,degrees,asin,radians,acos

def hf_state(ele_num,qubit_num):

    state = QuantumCircuit(qubit_num,qubit_num)
    machine = Aer.get_backend('statevector_simulator')
    for i in range(ele_num):
        state.x(i)
    
    result = execute(state,machine).result()
    statevector = result.get_statevector(state)
    return statevector
    #return state

def FT(ele_num,qubit_num,occ,vir):

    temp = occ / sqrt(occ**2 + vir**2)

    theta = radians(degrees(acos(temp))*2)

    state = QuantumCircuit(qubit_num,qubit_num)
    machine = Aer.get_backend('statevector_simulator')
    
    for i in range(ele_num):
        state.x(i)
        
    for i in range(qubit_num):
        state.ry(theta,i)

    result = execute(state,machine).result()
    statevector = result.get_statevector(state)

    return statevector
