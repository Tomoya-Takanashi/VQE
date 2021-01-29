import numpy
import copy
from qiskit import Aer,aqua,execute
from math import pi
from qiskit.visualization.utils import _validate_input_state

def statevector_expect(p_state,Opr,Num,qubit_num,h_num,leng):
    exp = []
    measure_circuit = 0
    machine = Aer.get_backend('statevector_simulator')
    measure_circuit = aqua.circuits.StateVectorCircuit(p_state)
    measure_circuit = measure_circuit.construct_circuit()
#    measure_circuit = p_state
#    print(p_state)

    if(h_num != leng-1):
        if(len(Opr) != 1):
            Opr_measure = Opr[0]
            for r in range(len(Opr_measure)):
                if(Opr_measure[r][1] == 'X'):
                    measure_circuit.h(Opr_measure[r][0])
                elif(Opr_measure[r][1] == 'Y'):
                    measure_circuit.rx(-pi/2,Opr_measure[r][0])
                    
            result = execute(measure_circuit,machine).result()
            count = result.get_statevector(measure_circuit,decimals=10)            
            rho = _validate_input_state(count)
            num_rho = int(numpy.log2(len(rho)))
            state = [bin(i)[2:].zfill(num_rho) for i in range(2**num_rho)]

            for r in range(len(count)):
                count[r] = abs(count[r])**2

            for r in range(len(Opr)):
                count_cp = 0
                state_cp = 0
                count_cp = copy.copy(count)
                state_cp = copy.copy(state)
                pauli_index = []
                for r2 in range(len(Opr[r])):
                    pauli_index.append(qubit_num - Opr[r][r2][0] - 1)
            
                                           
                for r2 in range(len(state_cp)):
                    new_state = 0
                    new_state = state_cp[r2][pauli_index[0]]
                    for r3 in range(1,len(pauli_index),1):
                        new_state = new_state + str(state_cp[r2][pauli_index[r3]])
                    state_cp[r2] = new_state

                for r2 in range(len(state_cp)):
                    one = state_cp[r2].count('1')
                    if(one % 2 != 0):
                        count_cp[r2] = count_cp[r2] * (-1)

                exp_val = sum(count_cp)
                exp.append(exp_val * Num[r])
                
        elif(len(Opr) == 1):

            Opr = Opr[0]
            for r in range(len(Opr)):
                if(Opr[r][1] == 'X'):
                    measure_circuit.h(Opr[r][0])
                elif(Opr[r][1] == 'Y'):
                    measure_circuit.rx(-pi/2,Opr[r][0])

                    #p_statevector = numpy.dot(basis_matrix,p_statevector.T)
            result = execute(measure_circuit,machine).result()
            count = result.get_statevector(measure_circuit,decimals=10)
            rho = _validate_input_state(count)
            num_rho = int(numpy.log2(len(rho)))
            state = [bin(i)[2:].zfill(num_rho) for i in range(2**num_rho)]

            for r in range(len(count)):
                count[r] = abs(count[r]) ** 2

            count_cp = copy.copy(count)
            state_cp = copy.copy(state)
            pauli_index = []
            
            for r2 in range(len(Opr)):
                pauli_index.append(qubit_num - Opr[r2][0] - 1)

            for r2 in range(len(state_cp)):
                new_state = 0
                new_state = state_cp[r2][pauli_index[0]]
                for r3 in range(1,len(pauli_index),1):
                    new_state = new_state + str(state_cp[r2][pauli_index[r3]])
                state_cp[r2] = new_state
            
                
            for r2 in range(len(state_cp)):
                one = state_cp[r2].count('1')
                if(one % 2 != 0):
                    count_cp[r2] = count_cp[r2] * (-1)

            exp_val = sum(count_cp)
            exp.append(exp_val * Num[0])
    else:
        #result = execute(measure_circuit,machine).result()
        #count = result.get_statevector(measure_circuit,decimals=10)
        count = copy.copy(p_state)
        rho = _validate_input_state(count)
        num_rho = int(numpy.log2(len(rho)))
        state = [bin(i)[2:].zfill(num_rho) for i in range(2**num_rho)]
        
                #print(output)
        for r in range(len(count)):
            count[r] = abs(count[r]) ** 2
            
        for r in range(len(Opr)):
            count_cp = 0
            state_cp = 0
            count_cp = copy.copy(count)
            state_cp = copy.copy(state)
            pauli_index = []
            if(len(Opr[r]) == 0):
                exp.append(Num[r])
            else:
                for r2 in range(len(Opr[r])):
                    pauli_index.append(qubit_num - Opr[r][r2][0] - 1)
                        #print(qubit_num)

                for r2 in range(len(state_cp)):
                    new_state = 0
                    new_state = state_cp[r2][pauli_index[0]]
                    for r3 in range(1,len(pauli_index),1):
                        new_state = new_state + str(state_cp[r2][pauli_index[r3]])
                    state_cp[r2] = new_state

                        #print(state_cp)
                for r2 in range(len(state_cp)):
                    one = state_cp[r2].count('1')
                    if(one % 2 != 0):
                        count_cp[r2] = count_cp[r2] * (-1)
                exp_val = sum(count_cp)
                exp.append(exp_val * Num[r])

    energy = sum(exp).real    
    
    return energy
