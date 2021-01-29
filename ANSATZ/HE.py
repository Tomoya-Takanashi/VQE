from qiskit import Aer,aqua,execute,QuantumCircuit
from math import pi

def he_circuit(qubit_num,theta_list,depth=1):
    machine = Aer.get_backend('statevector_simulator')
    #he = aqua.circuits.StateVectorCircuit(state)
    #circuit = he.construct_circuit()
    circuit = QuantumCircuit(qubit_num)
    for d in range(depth):
        for i in range(qubit_num):
            circuit.ry(theta_list[2*i+2*qubit_num*d],i)
            circuit.rz(theta_list[2*i+1+2*qubit_num*d],i)
        for i in range(qubit_num//2):
            circuit.cz(2*i, 2*i+1)
        for i in range(qubit_num//2-1):
            circuit.cz(2*i+1, 2*i+2)
    for i in range(qubit_num):
        circuit.ry(theta_list[2*i+2*qubit_num*depth],i)
        circuit.rz(theta_list[2*i+1+2*qubit_num*depth],i)

    result = execute(circuit,machine).result()
    statevector = result.get_statevector(circuit)
    return statevector

def he_circuit_type2(qubit_num,theta_list,depth=1):
    machine = Aer.get_backend('statevector_simulator')
    #he = aqua.circuits.StateVectorCircuit(state)
    #circuit = he.construct_circuit()
    circuit = QuantumCircuit(qubit_num)
    for i in range(qubit_num):
        circuit.rx(theta_list[3*i],i)
        circuit.rz(theta_list[3*i+1],i)
        circuit.rx(theta_list[3*i+2],i)
    for d in range(depth):
        for j in range(qubit_num//2):
            circuit.cz(2*j, 2*j+1)
        for j in range(qubit_num//2-1):
            circuit.cz(2*j+1, 2*j+2)
        for j in range(qubit_num):
            circuit.rx(theta_list[3*qubit_num*(d+1)+j*3],j)
            circuit.rz(theta_list[3*qubit_num*(d+1)+j*3+1],j)
            circuit.rx(theta_list[3*qubit_num*(d+1)+j*3+2],j)
    """
    for i in range(qubit_num):
        circuit.rx(theta_list[2*i+2*qubit_num*d],i)
        circuit.rz(theta_list[2*i+1+2*qubit_num*depth],i)
        circuit.rx(theta_list[2*i+2*qubit_num*d],i)
    """
    result = execute(circuit,machine).result()
    statevector = result.get_statevector(circuit)
    return statevector






