import numpy as np
from itertools import combinations


def rango_z2(matriz):
    if matriz is None or matriz.size == 0:
        return 0

    matriz = matriz.copy() % 2
    filas, columnas = matriz.shape
    rango = 0

    for col in range(columnas):
        for fila in range(rango, filas):
            if matriz[fila, col] == 1:
                if fila != rango:
                    matriz[[fila, rango]] = matriz[[rango, fila]]
                break
        else:
            continue

        for fila in range(rango + 1, filas):
            if matriz[fila, col] == 1:
                matriz[fila] = (matriz[fila] + matriz[rango]) % 2

        rango += 1

    return rango

def componentes_conexas_por_aristas_y_puntos(grafo):
    # Paso 1: Crear la lista de adyacencias
    vertices = {v[0] for v in grafo[0]}  # Extraemos el primer elemento de cada tupla
    adyacencias = {v: set() for v in vertices}  # Inicializamos cada vértice con un conjunto vacío

    for u, v in grafo[1]:  # Para cada arista (u, v)
        adyacencias[u].add(v)
        adyacencias[v].add(u)

    # Paso 2: Realizar una búsqueda BEP para encontrar componentes conexas
    visitados = set()
    def bep(v):
        stack = [v]
        while stack:
            nodo = stack.pop()
            if nodo not in visitados:
                visitados.add(nodo)
                stack.extend(adyacencias[nodo] - visitados)

    # Paso 3: Contar las componentes conexas
    componentes = 0
    for vertice in vertices:
        if vertice not in visitados:
            bep(vertice)
            componentes += 1

    return componentes

def calcular_homologia(complejo, dim):
    # Matriz de borde para dimensión actual (dim)
    matriz_borde = complejo.matriz_borde(dim)
    rango_borde = rango_z2(matriz_borde) if matriz_borde is not None else 0

    # Matriz de borde para dimensión superior (dim + 1)
    matriz_borde_superior = complejo.matriz_borde(dim + 1)
    rango_borde_superior = rango_z2(matriz_borde_superior) if matriz_borde_superior is not None else 0

    # Dimensión del núcleo (ker)
    dim_ker = len(complejo.simplices[dim]) - rango_borde

    # Dimensión del grupo de homología
    return dim_ker - rango_borde_superior


def is_simplice_positivo(simplex, N_i_menos_uno):
    dim = len(simplex) - 1
    N_i = {dim: set(simplices) for dim, simplices in N_i_menos_uno.items()}
    N_i[dim].add(simplex)

    if dim == 0:
        return True
    elif dim == 1:
        return componentes_conexas_por_aristas_y_puntos(N_i) == componentes_conexas_por_aristas_y_puntos(N_i_menos_uno)
    else:
        complejo_N_i = ComplejoSimplicial([])
        complejo_N_i.simplices = N_i
        complejo_N_i_menos_uno = ComplejoSimplicial([])
        complejo_N_i_menos_uno.simplices = N_i_menos_uno
        return calcular_homologia(complejo_N_i, dim) != calcular_homologia(complejo_N_i_menos_uno, dim)

class ComplejoSimplicial:
    def __init__(self, simplices_maximales):
        self.simplices = {}
        if simplices_maximales:
            self.generar_complejo(simplices_maximales)
        else:
            self.simplices[0] = []

    def generar_complejo(self, simplices_maximales):
        for simplex in simplices_maximales:
            simplex = tuple(sorted(simplex))  # Asegurar simplices ordenados
            for r in range(1, len(simplex) + 1):  # Generar caras de todos los tamaños
                for cara in combinations(simplex, r):
                    dim = len(cara) - 1
                    if dim not in self.simplices:
                        self.simplices[dim] = []
                    if cara not in self.simplices[dim]: 
                        self.simplices[dim].append(cara)

        # Convertimos los conjuntos a listas para consistencia
        for dim in self.simplices:
            self.simplices[dim] = list(self.simplices[dim])

    def matriz_borde(self, dimension):
        if dimension < 1 or dimension not in self.simplices or (dimension - 1) not in self.simplices:
            return None

        simplices_d = self.simplices[dimension]
        simplices_d_minus_1 = self.simplices[dimension - 1]

        n_rows = len(simplices_d_minus_1)
        n_cols = len(simplices_d)
        matriz = np.zeros((n_rows, n_cols), dtype=int)

        for j, simplex in enumerate(simplices_d):
            for i, face in enumerate(simplices_d_minus_1):
                if set(face).issubset(simplex):
                    matriz[i, j] = 1
        return matriz


    def calcular_betti(self):
        max_dim = max(self.simplices.keys()) if self.simplices else 0
        betti = []

        for dim in range(max_dim + 1):
            matriz_borde_dim = self.matriz_borde(dim)
            matriz_borde_dim_plus_1 = self.matriz_borde(dim + 1)

            n_dim = len(self.simplices.get(dim, []))
            rango_dim = rango_z2(matriz_borde_dim) if matriz_borde_dim is not None else 0
            rango_dim_plus_1 = rango_z2(matriz_borde_dim_plus_1) if matriz_borde_dim_plus_1 is not None else 0

            betti_dim = n_dim - rango_dim - rango_dim_plus_1
            betti.append(betti_dim)

        return betti

    def calcular_betti_incremental(self):
        """
        Calcula los números de Betti del complejo simplicial.
        """
        max_dim = max(self.simplices.keys(), default=-1)
        if max_dim == -1:
            return [0]
        betti = [0] * (max_dim + 1)
        procesados_por_dim = {dim: set() for dim in self.simplices}

        # Procesar símplices dimensión por dimensión
        for dim in range(max_dim + 1):

            for simplex in self.simplices[dim]:
                if is_simplice_positivo(simplex, procesados_por_dim):
                    betti[dim] += 1
                else:
                    betti[dim-1] -= 1

                procesados_por_dim[dim].add(simplex)

        return betti

    def calcular_betti_alfa_complejos(self, puntos, radio):
        simplices = []
        n = len(puntos)

        simplices.extend([(i,) for i in range(n)])

        for i in range(n):
            for j in range(i + 1, n):
                distancia = np.linalg.norm(np.array(puntos[i]) - np.array(puntos[j]))
                if distancia <= radio:
                    simplices.append((i, j))

        self.simplices = {}
        self.generar_complejo(simplices)
        return self.calcular_betti_incremental()


def to_string_p3(complejo:ComplejoSimplicial, nombre="Complejo Simplicial"):
    print("###########################################")
    print("--> ", nombre)
    dim = max(complejo.simplices.keys())
    for p in range(dim):
        print(f"- Matriz de borde dimension {p}:\n{complejo.matriz_borde(p+1)}")
    print("- Numeros betti (metodo normal): ", complejo.calcular_betti())
    print("- Numeros betti (algoritmo incremental): ", complejo.calcular_betti_incremental())
    print("###########################################\n")

if __name__ == "__main__":

    tetraedro = ComplejoSimplicial([(0,1,2,3)])
    to_string_p3(tetraedro, "Tetraedro")

    borde_tetraedro = ComplejoSimplicial([(0,1,2),(0,1,3),(0,2,3),(1,2,3)])
    to_string_p3(borde_tetraedro, "Borde Tetraedro")

    toro1 = ComplejoSimplicial([
    (1, 2, 4), (2, 4, 5), (2, 3, 5), (3, 5, 6), 
    (1, 3, 6), (1, 4, 6), (4, 5, 7), (5, 7, 8),
    (5, 6, 8), (6, 8, 9), (4, 6, 9), (4, 7, 9),
    (1, 7, 8), (1, 2, 8), (2, 8, 9), (2, 3, 9),
    (3, 7, 9), (1, 3, 7)])
    to_string_p3(toro1, "Toro: primera triangulacion")


    toro2 = ComplejoSimplicial([
    (1, 3, 7), (1,2,3), (2,5,7),(3,5,7),(3, 4, 5), (2, 5, 6), (2, 4, 7),
    (1, 2, 4), (1, 4,5), (1, 5, 6), (1, 6, 7),
    (4, 6, 7), (3, 4, 6),(2,3,6),])
    to_string_p3(toro2, "Toro: segunda triangulacion")


    plano_proyectivo = ComplejoSimplicial([(1, 2, 6), (2, 3, 4), (1, 3, 4), (1, 2, 5), (2, 3, 5), 
                        (1, 3, 6), (2, 4, 6), (1, 4, 5), (3, 5, 6), (4, 5, 6)])
    to_string_p3(plano_proyectivo, "Plano Proyectivo")

    botella_de_klein = ComplejoSimplicial([
    (0, 1, 3), (1, 3, 4), (1, 2, 4), (2, 5, 4),(2,0,5),(0,5,6),(3,4,6),(4,6,7),
    (4,5,7),(5,7,8),(5,6,8),(6,8,3),(6,7,0),(7,1,0),(7,8,1),(1,2,8),
    (8, 3, 2), (3, 2, 0)])
    to_string_p3(botella_de_klein, "Botella de Klein")


    anillo = ComplejoSimplicial([(1, 2, 4), (1, 3, 6), (1, 4, 6), (2, 3, 5), (2, 4, 5), (3, 5, 6)])
    to_string_p3(anillo, "Anillo")


    sombrero_de_burro = ComplejoSimplicial([
    (1, 2, 4), (1, 3, 4), (2, 3, 5), (2, 4, 5),
    (1, 3, 5), (1, 6, 5), (1, 3, 6), (3, 6, 7),
    (3, 2, 7), (1, 2, 7), (1, 7, 8), (1, 2, 8),
    (6, 7 ,8), (4, 5, 6), (4, 6, 8), (3, 4, 8), (2, 3, 8)])
    to_string_p3(sombrero_de_burro, "Sombrero de Burro")


    ejemplo_trans_4 = ComplejoSimplicial([(0, 1),(1, 2, 3, 4),(4, 5),(5, 6),(4, 6),(6, 7, 8),(8, 9)])
    to_string_p3(ejemplo_trans_4, "Ejemplo de la transparencia 4 (Homologia Simplicial II)")


    doble_toro = ComplejoSimplicial([
    (1, 2, 3), (1, 3, 4), (2, 3, 5), (3, 4, 5), (2, 5, 6),
    (4, 5, 7), (5, 6, 8), (5, 7, 8), (6, 8, 9), (7, 8, 10),
    (8, 9, 10), (1, 4, 7), (1, 6, 9), (1, 7, 10), (2, 6, 9),
    (2, 7, 10), (3, 8, 10), (3, 9, 10), (4, 9, 10), (5, 9, 10),
    (1, 2, 6), (1, 3, 9), (3, 6, 7), (6, 8, 10), (4, 5, 9),
    (2, 4, 8), (3, 7, 9), (5, 6, 10), (7, 8, 9),
    (4, 6, 8), (3, 5, 9)
    ])
    to_string_p3(doble_toro, "Doble Toro")

    # Ejemplo 1: Conjunto vacío
    alfa_complejo_1 = ComplejoSimplicial([])
    puntos1 = []
    radio1 = 0
    print("Números Betti alfa_complejo_1: ", alfa_complejo_1.calcular_betti_alfa_complejos(puntos1, radio1))

    # Ejemplo 2: Conjunto de puntos simples en línea (2D)
    alfa_complejo_2 = ComplejoSimplicial([])
    puntos2 = [[0, 0], [1, 0], [2, 0], [3, 0]]
    radio2 = 1.5
    print("Números Betti alfa_complejo_2: ", alfa_complejo_2.calcular_betti_alfa_complejos(puntos2, radio2))

    # Ejemplo 3: Cuadrado en 2D
    alfa_complejo_3 = ComplejoSimplicial([])
    puntos3 = [[0, 0], [1, 0], [0, 1], [1, 1]]
    radio3 = 1.0
    print("Números Betti alfa_complejo_3: ", alfa_complejo_3.calcular_betti_alfa_complejos(puntos3, radio3))

    # Ejemplo 4: Conjunto de puntos en 3D formando un tetraedro
    alfa_complejo_4 = ComplejoSimplicial([])
    puntos4 = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
    radio4 = 1.5
    print("Números Betti alfa_complejo_4: ", alfa_complejo_4.calcular_betti_alfa_complejos(puntos4, radio4))

    # Ejemplo 5: Conjunto denso de puntos en un plano 2D
    alfa_complejo_5 = ComplejoSimplicial([])
    puntos5 = [[0, 0], [0.5, 0], [1, 0], [0, 0.5], [0.5, 0.5], [1, 0.5], [0, 1], [0.5, 1], [1, 1]]
    radio5 = 0.75
    print("Números Betti alfa_complejo_5: ", alfa_complejo_5.calcular_betti_alfa_complejos(puntos5, radio5))

    # Ejemplo 6: Conjunto disperso de puntos en 3D
    alfa_complejo_6 = ComplejoSimplicial([])
    puntos6 = [[0, 0, 0], [2, 2, 0], [4, 0, 0], [2, 1, 3], [1, 1, 1]]
    radio6 = 2.5
    print("Números Betti alfa_complejo_6: ", alfa_complejo_6.calcular_betti_alfa_complejos(puntos6, radio6))
