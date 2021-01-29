import numpy
import copy
from qiskit import Aer,aqua,execute,QuantumCircuit
from math import pi

def qasm_expect(state,Opr,Num,shot_num,qubit_num,h_num,leng):
    exp = []
    machine = Aer.get_backend('qasm_simulator')
    measure_circuit = aqua.components.initial_states.Custom(qubit_num,state_vector=state)
    measure_circuit = measure_circuit.construct_circuit(mode='circuit')
    cbit_add = QuantumCircuit(qubit_num,qubit_num)
    if(h_num != leng-1):
        if(len(Opr) != 1):
            Opr_measure = Opr[0]
            for r in range(len(Opr_measure)):
                if(Opr_measure[r][1] == 'X'):
                    measure_circuit.h(Opr_measure[r][0])
                elif(Opr_measure[r][1] == 'Y'):
                    measure_circuit.rx(-pi/2,Opr_measure[r][0])
                cbit_add.measure(Opr_measure[r][0],Opr_measure[r][0])
            measure_circuit = measure_circuit.extend(cbit_add)
            
            result = execute(measure_circuit,machine,shots = shot_num).result()
            output = result.get_counts(measure_circuit)            
            state = list(output.keys())
            count = list(output.values())

            for r in range(len(count)):
                count[r] = count[r]/shot_num

            for r in range(len(Opr)):
                count_cp = copy.copy(count)
                state_cp = copy.copy(state)
                pauli_index = []
                for r2 in range(len(Opr[r])):
                    pauli_index.append(qubit_num - Opr[r][r2][0] - 1)
            
                                           
                for r2 in range(len(state_cp)):
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
                cbit_add.measure(Opr[r][0],Opr[r][0])
                    
            measure_circuit = measure_circuit.extend(cbit_add)
            result = execute(measure_circuit,machine,shots = shot_num).result()
            output = result.get_counts(measure_circuit)
            state = list(output.keys())
            count = list(output.values())
            for r in range(len(count)):
                count[r] = count[r]/shot_num

            for r in range(len(state)):
                one = state[r].count('1')
                if(one % 2 != 0):
                    count[r] = count[r] * (-1)

            exp_val = sum(count)
            exp.append(exp_val *Num[0])

    else:        
        for r in range(qubit_num):
            cbit_add.measure(r,r)
            
        measure_circuit = measure_circuit.extend(cbit_add)
        result = execute(measure_circuit,machine,shots = shot_num).result()
        output = result.get_counts(measure_circuit)
        state = list(output.keys())
        count = list(output.values())
                #print(output)
        for r in range(len(count)):
            count[r] = count[r]/shot_num
                #print(output)
       
        for r in range(len(Opr)):
            count_cp = copy.copy(count)
            state_cp = copy.copy(state)
            pauli_index = []
            if(r == 0):
                exp.append(Num[r])
            else:
                for r2 in range(len(Opr[r])):
                    pauli_index.append(qubit_num - Opr[r][r2][0] - 1)
                        #print(qubit_num)
                for r2 in range(len(state_cp)):
                    new_state = state_cp[r2][pauli_index[0]]
                    for r3 in range(1,len(pauli_index),1):
                        new_state = new_state + str(state_cp[r2][pauli_index[r3]])
                                #print(new_state)
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
