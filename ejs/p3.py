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
        max_dim = max(self.simplices.keys())
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

if __name__ == "__main__":
    sc1 = ComplejoSimplicial([(0,1,2,3)])
    print("Ejemplo 1 - Números de Betti:", sc1.calcular_betti_incremental())

    sc3 = ComplejoSimplicial([(0,1),(1,2,3,4),(4,5),(5,6),(4,6),(6,7,8),(8,9)])
    print("Ejemplo 3 - Números de Betti:", sc3.calcular_betti_incremental())

    sc5 = ComplejoSimplicial([(0,1,2),(2,3),(3,4)])
    print("Ejemplo 5 - Números de Betti:", sc5.calcular_betti_incremental())

    sc6 = ComplejoSimplicial([(1,2,4),(1,3,6),(1,4,6),(2,3,5),(2,4,5),(3,5,6)])
    print("Ejemplo 6 - Números de Betti:", sc6.calcular_betti_incremental())

    sc8 = ComplejoSimplicial([
    (1,2,4), (2,4,5), (2,3,5), (3,5,6), (1,3,6),
    (1,4,6), (4,5,7), (5,7,8),(5,6,8),(6,8,9),
    (4,6,9), (4,7,9), (1,7,8), (1,2,8), (2,8,9),
    (2,3,9), (3,7,9), (1,3,7)])
    print("Ejemplo 8 - Números de Betti:", sc8.calcular_betti_incremental())

    sc9 = ComplejoSimplicial([(1,2,6), (2,3,4), (1,3,4), (1,2,5), (2,3,5), (1,3,6), (2,4,6), (1,4,5), (3,5,6), (4,5,6)])
    print("Ejemplo 9 - Números de Betti:", sc9.calcular_betti_incremental())

    sc10 = ComplejoSimplicial([(0,), (1,), (2,3), (4,5), (5,6), (4,6), (6,7,8,9)])
    print("Ejemplo 10 - Números de Betti:", sc10.calcular_betti_incremental())

    unoesqtoro = ComplejoSimplicial([
    (1, 2), (1, 3), (1, 4), (1, 6), (1, 7), (1, 8),
    (2, 3), (2, 4), (2, 5), (2, 8), (2, 9), (3, 5),
    (3, 6), (3, 7), (3, 9), (4, 5), (4, 6), (4, 7),
    (4, 9), (5, 6), (5, 7), (5, 8), (6, 8), (6, 9),
    (7, 8), (7, 9), (8, 9)
    ])
    print("Numeros betti 1-esqueleto toro: ", unoesqtoro.calcular_betti_incremental())

    toro = ComplejoSimplicial([
    (1, 3, 7), (1,2,3), (2,5,7),(3,5,7),(3, 4, 5), (2, 5, 6), (2, 4, 7),
    (1, 2, 4), (1, 4,5), (1, 5, 6), (1, 6, 7),
    (4, 6, 7), (3, 4, 6),(2,3,6)
    ])
    print("Betti Toro: ", toro.calcular_betti_incremental())