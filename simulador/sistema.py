import numpy as np
from puertas import GATE_HADAMARD, GATE_I, GATE_SWAP
from bases import *
from scipy import sparse
from scipy.linalg import logm
from scipy.sparse.linalg import norm
import bloch as b

class SistemaMezcla:
    def __init__(self, numQbits = 1, lista = None): # Lista es una lista de listas de 2 elementos: probabilidad - vector de estados
        self._numQbits = numQbits
        dim = 2**self._numQbits
        self._proyectores = {}
        if lista is None:
            self._matrizDensidad = sparse.csr_matrix(([1.0 + 0j], ([0], [0])), 
                shape=(dim, dim), 
                dtype=complex)
            
        else:
            prob_total = 0.0
            rho = np.zeros((dim, dim), dtype=complex)

            for prob, estado in lista:
                prob_total += prob

                estado = np.asarray(estado, dtype=complex)

                if estado.size != dim:
                    raise ValueError(f"Todos los estados deben tener {dim} amplitudes")

                estado = estado.reshape(dim, 1)

                rho += (prob/100.0)*(estado @ estado.conj().T)

            if not np.isclose(prob_total, 100):
                raise ValueError("La probabilidad total debe ser 100")

            self._matrizDensidad = sparse.csr_matrix(rho)
        
    # -- MODULO DE PUERTAS

    def superposicion(self):
        for i in range(self._numQbits):
            self.aplicar_puerta(GATE_HADAMARD, [i])

    def aplicar_puerta(self, GATE, POS_QBITS = []):
        OP = sparse.eye(1)
        GATE_SP = sparse.csr_matrix(GATE)
        I_SP = sparse.csr_matrix(GATE_I)

        # -- EXCEPTION HANDLING FOR "aplicar_puerta"

        if int(np.log2(GATE.shape[0])) != len(POS_QBITS):
            raise ValueError(f"More qubits given than needed for a {np.log2(GATE.shape[0])} qubits gate")
    
        if [elem for elem in POS_QBITS if elem >= self._numQbits]:
            raise ValueError(f"Qubit out of system range")
        
        if len(POS_QBITS) > self._numQbits:
            raise ValueError(f"The gate is to big for this system")
        
        if len(POS_QBITS) != len(set(POS_QBITS)):
            raise ValueError(f"All the qubits must be different")

        if len(POS_QBITS) == 1:			
            for i in range(self._numQbits):
                if i == POS_QBITS[0]:
                    OP = sparse.kron(OP, GATE_SP, format='csr')
                else:
                    OP = sparse.kron(OP, I_SP, format='csr')
                    
            self._matrizDensidad = OP @ self._matrizDensidad @ OP.conj().T
        
        else:			
            MAPA_QUBITS = list(range(self._numQbits)) # Como están colocados los qubits inicialmente(cada uno en su posición correcta, el qubit 0 en el 0....)
            NUM_QBITS = len(POS_QBITS)
            SWAPS = []
            
            POS_FINALES = list(range(self._numQbits - NUM_QBITS, self._numQbits)) ## Posicion que cada qubit tiene que ocupar
            
            for i in range(NUM_QBITS):
                posQ = POS_QBITS[i]
                posF = POS_FINALES[i]
                
                posI = MAPA_QUBITS.index(posQ)
                
                if posI == posF:
                    continue

                self.swap_qbits(posI, posF)
                SWAPS.append([posI, posF])
                MAPA_QUBITS[posI], MAPA_QUBITS[posF] = MAPA_QUBITS[posF], MAPA_QUBITS[posI];
                             
            i = 0
            while i < self._numQbits:
                if i == self._numQbits - NUM_QBITS:
                    OP = sparse.kron(OP, GATE_SP, format='csr')
                    i += len(POS_QBITS)
                else:
                    OP = sparse.kron(OP, I_SP, format='csr')
                    i += 1
                    
            self._matrizDensidad = OP @ self._matrizDensidad @ OP.conj().T
            
            for pQ, pF in reversed(SWAPS):
                self.swap_qbits(pF, pQ)

    def swap_qbits(self, Qi, Qf):
        
        if (Qi < 0 or Qi >= self._numQbits) or (Qf < 0 or Qf >= self._numQbits):
            raise ValueError("Qubit out of system range")

        posAct = Qi

        if Qf == Qi:
            return
        
        direccion = 1
        if Qf-Qi < 0:
            direccion = -1

        while Qf-posAct != 0:
            self._swap_qbits_ady(posAct, posAct+direccion)
            posAct = posAct + direccion
        
        posAct -= direccion
        while posAct-Qi != 0:
            self._swap_qbits_ady(posAct, posAct-direccion)
            posAct = posAct-direccion

    def _swap_qbits_ady(self, Qi, Qf):
        OP = sparse.eye(1)
        SWP_SP = sparse.csr_matrix(GATE_SWAP)
        I_SP = sparse.csr_matrix(GATE_I)
        obj = min(Qi, Qf)
        i = 0
        while i < self._numQbits:
            if i == obj:
                OP = sparse.kron(OP, SWP_SP, format='csr')
                i += 2
            else:
                OP = sparse.kron(OP, I_SP, format='csr')
                i += 1

        self._matrizDensidad = OP @ self._matrizDensidad @ OP.conj().T
    
    # -- MODULO MATRIZ DE DENSIDAD Y OTROS
    
    def print_matriz_densidad(self):
        print(np.round(self.matriz_densidad(), 5))

    def print_matriz_densidad_dispersa(self):
        print(self._matrizDensidad)

    def matriz_densidad(self):
        return self._matrizDensidad.toarray()
    
    def entropia_von_neuman(self):
        # Para calcular el logaritmo de la matriz usamos el modulo logm de linalg 
        rho = self.matriz_densidad() 
        return np.real(-np.trace(rho @ logm(rho)))
    
    # -- MODULO DE MEDICIONES

    def medir(self, base_medicion=BASE_Z):
        mediciones = []
        for i in range(self._numQbits):
            mediciones.append(self.medir_qbit(i, base_medicion))

        return mediciones

    def medir_qbit(self, pos_qbit=0, base_medicion=BASE_Z):

        for GATE in base_medicion:
            self.aplicar_puerta(GATE, [pos_qbit])

        rho = self._matrizDensidad

        OP_P0 = self._generar_proyector(pos_qbit, 0) # Creamos los operadores de proyeccion
        OP_P1 = self._generar_proyector(pos_qbit, 1)
        
        prob_0 = (OP_P0 @ rho).trace().real #  p(v) = tr(Pi,v |Ψi) -> Cogemos solo la parte real
        prob_1 = (OP_P1 @ rho).trace().real

        valor = 1 if np.random.random() < prob_1 else 0

        prob = prob_1 if valor == 1 else prob_0 # Probabilidad calculada
        Piv = OP_P1 if valor == 1 else OP_P0 # Proyector utilizado

        if prob == 0:
            raise ValueError("Division by zero")
        
        self._matrizDensidad = (Piv @ rho @ Piv)/prob # Aplicamos la operación para calcular la nueva matriz de densidad

        return valor
    
    def _generar_proyector(self, pos_qbit, valor):
        clave = (pos_qbit, valor)
        if clave in self._proyectores:
            return self._proyectores[clave]
        
        P0 = sparse.csr_matrix([[1, 0], [0, 0]])
        P1 = sparse.csr_matrix([[0, 0], [0, 1]])
        P = P1 if valor == 1 else P0
        
        res = sparse.eye(1)
        for i in range(self._numQbits):
            if i == pos_qbit:
                res = sparse.kron(res, P, format='csr')
            else:
                res = sparse.kron(res, sparse.eye(2), format='csr')

        self._proyectores[clave] = res

        return res
    

class SistemaBasico:
    def __init__(self, numQbits = 1, estadoInicial = None):
        self._numQbits = numQbits
        dim = 2**numQbits
        if(estadoInicial is not None):
            self._estados = sparse.csr_matrix(estadoInicial, dtype=complex)

            if self._estados.shape[0] == 1:
                self._estados = self._estados.T

            self._numQbits = int(np.log2(self._estados.shape[0]))

            n = norm(self._estados)
            
            if not np.isclose(n, 1.0):
                raise ValueError(f"Error al crear el sistema: El estado no está normalizado. Norma actual: {n}. La suma de |a|^2 debe ser 1.")
        else:
            self._estados = sparse.csr_matrix(([1.0 + 0j], ([0], [0])), 
                shape=(dim, 1), 
                dtype=complex)

    # -- MODULO DE ATRIBUTOS

    def print_estado(self, num_decimales = 5, prev = ""):
        estado_c = 0
        estado_S = ""
        filas, _ = self._estados.nonzero() 
        amplitudes = self._estados.data

        for i, fila in enumerate(filas):
            print(fila)
            if estado_S == "":
                estado_S = f"{prev}{round(np.real(amplitudes[i]), num_decimales)}|{fila:0{self._numQbits}b}>"
            else:
                estado_S = f"{estado_S} + {round(np.real(amplitudes[i]), num_decimales)}|{fila:0{self._numQbits}b}>"

        
        print(estado_S)

    def get_estados(self):
        return self._estados

    def get_numQbits(self):
        return self._numQbits
    
    # -- MODULO DE MEDICIONES

    def _probabilidades(self): # Función privada por convención
        return np.abs(self._estados.data)**2

    def print_probabilidades(self):
        
        filas, _ = self._estados.nonzero()               
        probs = self._probabilidades()
        
        for i, fila in enumerate(filas):
            print(f"|{fila:0{self._numQbits}b}>: {100 * round(probs[i], 5)}%")

    def medir_qbit(self, pos_qbit=0, base_medicion=BASE_Z):
        for GATE in base_medicion:
            self.aplicar_puerta(GATE, [pos_qbit])

        filas, _ = self._estados.nonzero()
        probs = self._probabilidades()
        
        prob_1 = 0.0
        bit_pos = self._numQbits-1-pos_qbit

        for i, idx_estado in enumerate(filas):
            if (idx_estado >> bit_pos) & 1:
                prob_1 += probs[i]
        
        valor = 1 if np.random.random() < prob_1 else 0
        self.colapsar(pos_qbit, valor)
        return valor
    
    def medir(self, base_medicion=BASE_Z):
        mediciones = []
        for i in range(self._numQbits):
            mediciones.append(self.medir_qbit(i, base_medicion))

        return mediciones

    def colapsar(self, pos_qbit=0, valor=0):
        bit_pos = self._numQbits-1-pos_qbit
        
        filas, _ = self._estados.nonzero()
        data = self._estados.data
        
        mask = ((filas >> bit_pos) & 1) == valor
        
        filas_vivas = filas[mask]
        datos_vivos = data[mask]
        
        self._estados = sparse.csr_matrix(
            (datos_vivos, (filas_vivas, np.zeros(len(filas_vivas)))),
            shape=(2**self._numQbits, 1), dtype=complex)
        
        n = norm(self._estados)
        if n == 0:
            raise ValueError("Colapso imposible: la probabilidad del resultado era 0")
            
        self._estados = self._estados/n
        
    # -- MODULO DE PUERTAS

    def superposicion(self):
        for i in range(self._numQbits):
            self.aplicar_puerta(GATE_HADAMARD, [i])

    def aplicar_puerta(self, GATE, POS_QBITS = []):
        OP = sparse.eye(1)
        GATE_SP = sparse.csr_matrix(GATE)
        I_SP = sparse.csr_matrix(GATE_I)
        
        # -- EXCEPTION HANDLING FOR "aplicar_puerta"

        if int(np.log2(GATE.shape[0])) != len(POS_QBITS):
            raise ValueError(f"La puerta de {int(np.log2(GATE.shape[0]))} qubits no coincide con la cantidad de qubits proporcionada ({len(POS_QBITS)})")
    
        if [elem for elem in POS_QBITS if elem >= self._numQbits]:
            raise ValueError(f"Índice de qubit fuera de rango: asegúrate de que todos los índices estén entre 0 y {self._numQbits-1}")
        
        if len(POS_QBITS) > self._numQbits:
            raise ValueError(f"La puerta es demasiado grande para el sistema")
        
        if len(POS_QBITS) != len(set(POS_QBITS)):
            raise ValueError(f"Todos los qubits seleccionados deben ser distintos")
        
        if len(POS_QBITS) == 1:
            for i in range(self._numQbits):
                if i == POS_QBITS[0]:
                    OP = sparse.kron(OP, GATE_SP, format='csr')
                else:
                    OP = sparse.kron(OP, I_SP, format='csr')

            self._estados = OP @ self._estados
        
        else: # A partir de 2 qubits solo se pueden usar puertas condicionales
            MAPA_QUBITS = list(range(self._numQbits)) # Como están colocados los qubits inicialmente(cada uno en su posición correcta, el qubit 0 en el 0....)
            NUM_QBITS = len(POS_QBITS)
            SWAPS = []
            
            POS_FINALES = list(range(self._numQbits - NUM_QBITS, self._numQbits)) ## Posicion que cada qubit tiene que ocupar
            
            for i in range(NUM_QBITS):
                posQ = POS_QBITS[i]
                posF = POS_FINALES[i]
                
                posI = MAPA_QUBITS.index(posQ)
                
                if posI == posF:
                    continue

                self.swap_qbits(posI, posF)
                SWAPS.append([posI, posF])
                MAPA_QUBITS[posI], MAPA_QUBITS[posF] = MAPA_QUBITS[posF], MAPA_QUBITS[posI];
                
            i = 0
            while i < self._numQbits:
                if i == self._numQbits - len(POS_QBITS):
                    OP = sparse.kron(OP, GATE_SP, format='csr')
                    i += len(POS_QBITS)
                else:
                    OP = sparse.kron(OP, I_SP, format='csr')
                    i += 1
            self._estados = OP @ self._estados

            for pQ, pF in reversed(SWAPS):
                self.swap_qbits(pF, pQ)


    def swap_qbits(self, Qi, Qf):
        if (0 > Qi or Qi >= self._numQbits) or (0 > Qf or Qf >= self._numQbits):
            raise ValueError(f"Qubits fuera de rango")

        posAct = Qi

        if Qf == Qi:
            return []
        
        direccion = 1
        if Qf-Qi < 0:
            direccion = -1

        while Qf-posAct != 0:
            self._swap_qbits_ady(posAct, posAct+direccion)
            posAct = posAct + direccion
        
        posAct -= direccion
        while posAct-Qi != 0:
            self._swap_qbits_ady(posAct, posAct-direccion)
            posAct = posAct-direccion


    def _swap_qbits_ady(self, Qi, Qf): # Función privada por convención	
        
        GATE_SWAP_SP = sparse.csr_matrix(GATE_SWAP)
        I_SP = sparse.csr_matrix(GATE_I)

        OP = sparse.eye(1)
        obj = min(Qi, Qf)
        i = 0
        while i < self._numQbits:
            if i == obj:
                OP = sparse.kron(OP, GATE_SWAP_SP, format='csr')
                i += 2
            else:
                OP = sparse.kron(OP, I_SP, format='csr')
                i += 1

        self._estados = OP @ self._estados

    # -- MODULO ESFERA/COORDS BLOCH

    def angulos_bloch(self):
        if self._numQbits > 1:
            raise ValueError("No se puede obtener coordenadas de Bloch para más de un qubit")
                
        A0 = self._estados[0, 0]
        A1 = self._estados[1, 0]

        z = abs(A0)**2 - abs(A1)**2
        theta = np.arccos(z)
        phi = np.angle(A1) - np.angle(A0)
        phi %= 2*np.pi

        return theta, phi
    
    def coords_bloch(self):
        if self._numQbits > 1:
            raise ValueError("No se puede obtener coordenadas de Bloch para más de un qubit")
        
        A0 = self._estados[0, 0]
        A1 = self._estados[1, 0]

        z = abs(A0)**2 - abs(A1)**2
        x = round(2*np.real(np.conj(A0)*A1), 5)
        y = round(2*np.imag(np.conj(A0)*A1), 5)

        return x, y, z
    
    # -- MODULO OTROS

    def producto_tensorial(self, sistema2):
        nuevoEstado = sparse.kron(self._estados, sistema2.get_estados(), format='csr')

        return SistemaBasico(int(np.log2(nuevoEstado.shape[0])), nuevoEstado)