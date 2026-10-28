import numpy as np

#1 QUBIT
GATE_I = np.eye(2, dtype=complex)
GATE_X = np.array([[0, 1],[1, 0]], dtype=complex)
GATE_SR_NOT = (1/2) * np.array([[1+1j, 1-1j],[1-1j, 1+1j]], dtype=complex) #SQUARE ROOT NOT
GATE_Y = np.array([[0, -1j],[1j, 0]], dtype=complex)
GATE_Z = np.array([[1, 0],[0, -1]], dtype=complex)
GATE_S = np.array([[1, 0],[0, 1j]], dtype=complex) #Rotación 90grados al rededor de eje Z. Es la raiz cuadrada de la Pauli-Z
GATE_HADAMARD = (1/np.sqrt(2)) * np.array([[1, 1],[1, -1]], dtype=complex)

#2 QUBITS
GATE_SWAP = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], dtype=complex)
GATE_SR_SWAP = np.array([[1,0,0,0],[0,(1/2)*(1+1j),(1/2)*(1-1j),0],[0,(1/2)*(1-1j),(1/2)*(1+1j),0],[0,0,0,1]], dtype=complex) #SQUARE ROOT SWAP

#Multicontroladas:
GATE_CNOT = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]], dtype=complex)
GATE_MCY = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,-1j],[0,0,1j,0]], dtype=complex)
GATE_MCZ = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]], dtype=complex)

#3 QUBITS
GATE_FREDKIN = np.array([ 
        [1,0,0,0,0,0,0,0], 
        [0,1,0,0,0,0,0,0], 
        [0,0,1,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,1,0], # 101 -> 110
        [0,0,0,0,0,1,0,0], # 110 -> 101
        [0,0,0,0,0,0,0,1] 
    ], dtype=complex)

GATE_TOFFOLI = np.array([ 
        [1,0,0,0,0,0,0,0], 
        [0,1,0,0,0,0,0,0], 
        [0,0,1,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,1,0,0,0],
        [0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,1], # 110 -> 111
        [0,0,0,0,0,0,1,0]  # 111 -> 110
    ], dtype=complex)

GATE_MCZ3 = np.array([ 
        [1,0,0,0,0,0,0,0], 
        [0,1,0,0,0,0,0,0], 
        [0,0,1,0,0,0,0,0],
        [0,0,0,1,0,0,0,0],
        [0,0,0,0,1,0,0,0],
        [0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,-1j], # 110 -> 111
        [0,0,0,0,0,0,1j,0]  # 111 -> 110
    ], dtype=complex)