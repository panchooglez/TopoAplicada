import numpy as np
from itertools import combinations

class ComplejoSimplicial:
    def __init__(self, simplices_maximales):
        self.simplices = {}
        if simplices_maximales:
            self.generar_complejo(simplices_maximales)
        else:
            self.simplices[0] = []

    def generar_complejo(self, simplices_maximales):
        if len(simplices_maximales) != 0:
            for simplex in simplices_maximales:
                simplex = tuple(sorted(simplex))
                for r in range(1, len(simplex) + 1):
                    for cara in combinations(simplex, r):
                        dim = len(cara) - 1
                        if dim not in self.simplices:
                            self.simplices[dim] = []
                        if cara not in self.simplices[dim]:
                            self.simplices[dim].append(cara)

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
                    # Contamos el número de vértices que necesitamos quitar para llegar de simplex a face
                    # Si es par, ponemos un 1; si es impar, ponemos un -1 (pero en Z2 es lo mismo)
                    matriz[i, j] = 1
        return matriz

    def rango_z2(self, matriz):
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

    def calcular_betti(self):
        max_dim = max(self.simplices.keys()) if self.simplices else 0
        betti = []

        for dim in range(max_dim + 1):
            matriz_borde_dim = self.matriz_borde(dim)
            matriz_borde_dim_plus_1 = self.matriz_borde(dim + 1)

            n_dim = len(self.simplices.get(dim, []))
            rango_dim = self.rango_z2(matriz_borde_dim) if matriz_borde_dim is not None else 0
            rango_dim_plus_1 = self.rango_z2(matriz_borde_dim_plus_1) if matriz_borde_dim_plus_1 is not None else 0

            betti_dim = n_dim - rango_dim - rango_dim_plus_1
            betti.append(betti_dim)

        return betti

    def calcular_betti_incremental(self):
        componentes = {}
        parent = {}

        def find(v):
            while parent[v] != v:
                parent[v] = parent[parent[v]]
                v = parent[v]
            return v

        ciclos = 0

        for simplex in sorted(self.simplices.get(0, [])):
            parent[simplex[0]] = simplex[0]

        for edge in sorted(self.simplices.get(1, [])):
            u, v = edge
            if u not in parent:
                parent[u] = u
            if v not in parent:
                parent[v] = v
            pu, pv = find(u), find(v)
            if pu != pv:
                parent[pu] = pv
            else:
                ciclos += 1

        b_0 = len({find(v) for v in parent})
        b_1 = ciclos
        return [b_0, b_1]

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
