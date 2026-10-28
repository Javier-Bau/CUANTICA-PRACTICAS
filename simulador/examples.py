from sistema import *
from puertas import *
from bases import *
import bloch as b

def prueba_toffoli():
    """
    Dado un circuito de 4 qubits vamos a probar a usar la toffoli con varias combinaciones
    """
    print("############")
    print("############")
    print("TOFFOLI TEST")
    print("############")
    print("############")

    sistema = SistemaBasico(4)
    print("\nSistema inicial:")
    sistema.print_estado()
    sistema.aplicar_puerta(GATE_TOFFOLI, [0,1,2]) # No va a pasar nada
    print("\nSistema trás aplicar Toffoli sobre qubits 0,1(CTRL) y 2:")
    sistema.print_estado()
    sistema.aplicar_puerta(GATE_X, [0])
    sistema.aplicar_puerta(GATE_X, [1])
    print("\nAplicamos GATE_X qubits 0,1:")
    sistema.print_estado()
    sistema.aplicar_puerta(GATE_TOFFOLI, [0,1,2]) 
    print("\nSistema trás aplicar Toffoli sobre qubits 0,1(CTRL) y 2:")
    sistema.print_estado()



def bernstein_vazirani():
    """
    Con un sistema de 5 qubits(+ uno auxiliar) buscamos mediante un oraculo marcar el estado 10100
    """
    print("##################")
    print("##################")
    print("BERNSTEIN VAZIRANI")
    print("##################")
    print("##################")
    sistema = SistemaBasico(6)
    sistema.aplicar_puerta(GATE_X, [5])
    sistema.superposicion()
    print("\nSistema antes de marcar el estado:")
    sistema.print_estado()
    sistema.aplicar_puerta(GATE_CNOT, [0,5])
    sistema.aplicar_puerta(GATE_CNOT, [2,5])
    print("\nSistema después de marcar el estado pero en superposición aún:")
    sistema.print_estado()
    sistema.superposicion()
    print("\nSistema después de quitar la superposición:")
    sistema.print_estado()
    """
    Vemos como el único estado con amplitud es el estado 10100(el qubit de más a la derecha es auxiliar)
    """


def prueba_toffoli_mezcla():
    """
    Dado un circuito de 4 qubits vamos a probar a usar la toffoli con varias combinaciones
    """
    print("############")
    print("############")
    print("TOFFOLI TEST")
    print("############")
    print("############")

    sistema = SistemaMezcla(3) # Es sistema de estados mezcla pero solo tiene un estado puro para que se vea bien el ejemplo
    print("\nSistema inicial:")
    sistema.print_matriz_densidad()
    sistema.aplicar_puerta(GATE_TOFFOLI, [0,1,2]) # No va a pasar nada
    print("\nSistema trás aplicar Toffoli sobre qubits 0,1(CTRL) y 2:")
    sistema.print_matriz_densidad()
    sistema.aplicar_puerta(GATE_X, [0])
    sistema.aplicar_puerta(GATE_X, [1])
    print("\nAplicamos GATE_X qubits 0,1:")
    sistema.print_matriz_densidad()
    sistema.aplicar_puerta(GATE_TOFFOLI, [0,1,2]) 
    print("\nSistema trás aplicar Toffoli sobre qubits 0,1(CTRL) y 2:")
    sistema.print_matriz_densidad()

def bernstein_vazirani_mezcla():
    """
    Con un sistema de 5 qubits(+ uno auxiliar) buscamos mediante un oraculo marcar el estado 10
    """
    print("##################")
    print("##################")
    print("BERNSTEIN VAZIRANI")
    print("##################")
    print("##################")
    sistema = SistemaMezcla(3) # Es sistema de estados mezcla pero solo tiene un estado puro para que se vea bien el ejemplo
    print("Buscamos marcar el estado 10:\n")
    sistema.aplicar_puerta(GATE_X, [2])
    sistema.superposicion()
    print("\nSistema antes de marcar el estado:")
    sistema.print_matriz_densidad()
    sistema.aplicar_puerta(GATE_CNOT, [0,2])
    print("\nSistema después de marcar el estado pero en superposición aún:")
    sistema.print_matriz_densidad()
    sistema.superposicion()
    print("\nSistema después de quitar la superposición(Matriz de densidad en formato disperso):")
    sistema.print_matriz_densidad()
    """
    Vemos como el único estado con amplitud es el estado 101(el qubit de más a la derecha es auxiliar)
    """

def prueba_von_neumann():
    print("##################")
    print("##################")
    print("PRUEBA VON NEUMANN")
    print("##################")
    print("##################")
    s1 = SistemaMezcla(1, [[50, [1,0]],[50, [1,0]]])
    s2 = SistemaMezcla(1, [[50, [1,0]],[50, [0,1]]])
    s3 = SistemaMezcla(1)
    print("\n\nSistema1: ")
    s1.print_matriz_densidad()
    print(f"Entropia s3: {s1.entropia_von_neuman()}")
    print("\n\nSistema2: ")
    s2.print_matriz_densidad()
    print(f"Entropia s3: {s2.entropia_von_neuman()}")
    print("\n\nSistema3: ")
    s3.print_matriz_densidad()
    print(f"Entropia s3: {s3.entropia_von_neuman()}")
    print("\n\nSistema3 + Hadamards: ")
    s3.superposicion()
    s3.print_matriz_densidad()
    print(f"Entropia s3: {s3.entropia_von_neuman()}")

def prueba_deutschJozsa():
    print("####################")
    print("####################")
    print("PRUEBA DEUTSCH JOZSA")
    print("####################")
    print("##### PRUEBA 1 #####")
    print("####################")

    sistema = SistemaBasico(4)
    sistema.aplicar_puerta(GATE_X, [3])
    sistema.superposicion()
    print("\nSistema antes de oraculo: ")
    sistema.print_estado()
    ## ORACULO:
    sistema.aplicar_puerta(GATE_CNOT, [2, 3])
    sistema.aplicar_puerta(GATE_CNOT, [1, 3])
    sistema.aplicar_puerta(GATE_TOFFOLI, [0, 1, 3])
    ## FIN ORACULO
    sistema.aplicar_puerta(GATE_HADAMARD, [0])
    sistema.aplicar_puerta(GATE_HADAMARD, [1])
    sistema.aplicar_puerta(GATE_HADAMARD, [2])
    print("\nSistema después de oraculo y quitar interferencias:")
    sistema.print_estado()
    res = sistema.medir()
    resultado = "Balanceada" if 1 in res[:-1] else "Constante"
    print(f"\nResultado medición: {res} -> Función {resultado}")

    input("\n\nPulsa [ENTER] para siguiente prueba DJ: ")
    print("\n\n####################")
    print("##### PRUEBA 2 #####")
    print("####################")
    sistema = SistemaBasico(4)
    sistema.aplicar_puerta(GATE_X, [3])
    sistema.superposicion()
    print("\nSistema antes de oraculo: ")
    sistema.print_estado()
    ## ORACULO:
    sistema.aplicar_puerta(GATE_CNOT, [0, 3])
    sistema.aplicar_puerta(GATE_X, [0])
    sistema.aplicar_puerta(GATE_CNOT, [0, 3])
    sistema.aplicar_puerta(GATE_CNOT, [1, 3])
    sistema.aplicar_puerta(GATE_X, [1])
    sistema.aplicar_puerta(GATE_CNOT, [1, 3])
    sistema.aplicar_puerta(GATE_CNOT, [2, 3])
    sistema.aplicar_puerta(GATE_X, [2])
    sistema.aplicar_puerta(GATE_CNOT, [2, 3])
    ## FIN ORACULO
    sistema.aplicar_puerta(GATE_HADAMARD, [0])
    sistema.aplicar_puerta(GATE_HADAMARD, [1])
    sistema.aplicar_puerta(GATE_HADAMARD, [2])
    print("\nSistema después de oraculo y quitar interferencias:")
    sistema.print_estado()
    res = sistema.medir()
    resultado = "Balanceada" if 1 in res[:-1] else "Constante"
    print(f"\nResultado medición: {res} -> Función {resultado}")
    
def grover():
    
    ORACULO = np.array( 
        [[1,0,0,0],
         [0,1,0,0],
         [0,0,-1,0],
         [0,0,0,1]  # Queremos marcar estado 10
        ], dtype=complex)

    sistema = SistemaBasico(2)
    sistema.superposicion()
    print(f"\nSistema inicial: ")
    sistema.print_estado()

    sistema.aplicar_puerta(ORACULO, [0,1])
    print("\nSistema tras Oráculo: ")
    sistema.print_estado()
    sistema.superposicion()
    sistema.aplicar_puerta(GATE_X, [0])
    sistema.aplicar_puerta(GATE_X, [1])
    sistema.aplicar_puerta(GATE_MCZ, [0,1])
    sistema.aplicar_puerta(GATE_X, [0])
    sistema.aplicar_puerta(GATE_X, [1])
    sistema.superposicion()
    print(f"\nSistema tras inversión sobre la media: ")
    sistema.print_estado()
    print(f"Medición: {sistema.medir()}")


def main():
    print("SISTEMA BÁSICO:")
    input("\n\nPulsa [ENTER] para ver \"prueba Deutsch Jozsa\"")
    prueba_deutschJozsa()
    input("\n\nPulsa [ENTER] para ver \"prueba Bernstein Vazirani\"")
    bernstein_vazirani()
    input("\n\nPulsa [ENTER] para ver \"prueba grover\"")
    grover()
    input("\n\nPulsa [ENTER] para ver \"prueba toffoli\"")
    prueba_toffoli()
    input("\n\nPulsa [ENTER] para ver \"SistemaMezcla\"")
    print("SISTEMA MEZCLA:")
    input("\n\nPulsa [ENTER] para ver \"prueba Bernstein Vazirani\"")
    bernstein_vazirani_mezcla()
    input("\n\nPulsa [ENTER] para ver \"prueba toffoli\"")
    prueba_toffoli_mezcla()
    input("\n\nPulsa [ENTER] para ver \"prueba Von Neumann\"")
    prueba_von_neumann()

if __name__ == "__main__":
    main()
