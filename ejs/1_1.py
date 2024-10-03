from itertools import chain, combinations

class ComplejoSimplicial:

    """
    Inicializa el complejo simplicial tomando como argumento lo que serían los maximales del poset
    """
    def __init__(self, sim_maximales):
        self.simplices = dict()
        sim_derivados = self.derivar_simplices(sim_maximales)
        for sim in sim_derivados:
            if isinstance(sim, int):
                dim = 0
            else:    
                dim = len(sim)-1
            if dim not in self.simplices:
                self.simplices[dim] = set()
            self.simplices[dim].add(sim)

    def derivar_simplices(self, simplices_maximales):
        """
        Genera el complejo simplicial tomando todos los subconjuntos de los símplices maximales.
        
        Args:
        simplices_maximales (list of tuples): Lista de símplices maximales.
        
        Returns:
        set: El conjunto de todos los símplices en el complejo simplicial.
        """
        complejo = set()
        for simplejo in simplices_maximales:
            complejo.update(self.potencia(simplejo))
        return complejo
    
    def potencia(self, simplejo):
        """
        Genera todos los subconjuntos de un simplejo (incluido el propio simplejo).
        
        Args:
        simplejo (tuple): Un simplejo representado como una tupla de vértices.
        
        Returns:
        set: Conjunto de todos los subconjuntos del simplejo.
        """
        return set(chain.from_iterable(combinations(simplejo, r) for r in range(len(simplejo) + 1)))

    """
    Devuelve la dimension del complejo simplicial
    """
    def dim(self):
        return max(self.simplices.keys())

    """
    Devuelve una lista con las caras de un complejo simplicial para una dimension dada
    """
    def caras(self, dim):
        caras = []
        for cara in self.simplices[dim]:
            caras.append(cara)
        return caras

    """
    Devuelve la estrella generada por un simplice dentro del complejo simplicial
    """
    def estrella(self, simplice):
        pass

    """
    Devuelve el link del simplice que toma como parametro
    """
    def link(self, simplice):
        pass

    """
    Permite printear el complejo simplicial
    """
    def __str__(self):
        pass

    def print_caras(self):
        dim = self.dim()
        for i in range(0,dim+1):
            print("Caras de dim: ", i)
            d_caras = self.caras(i)
            for cara in d_caras:
                print(cara)
    
    def carac_euler(self):
        dim = self.dim()
        res = 0
        for i in range(0, dim+1):
            if i % 2 == 0:
                res += len(self.simplices[dim])
            else:
                res -= len(self.simplices[dim])
        return res
    
    def num_componentes_conexas(self):
        n = 0
        vertices = list(self.simplices[0])
        aristas = self.simplices[1]

        v_conexos = []
        v_inconexos = []
        v0 = vertices[0]

        while True:
            print('hola')
            for v in vertices[:]:
                print(f"v0: {v0}, v: {v}, existe_camino: {self.existe_camino(v0[0], v[0], aristas)}")
                if self.existe_camino(v0[0], v[0], aristas):
                    v_conexos.append(v)
                    vertices.remove(v)
                else:
                    v_inconexos.append(v)
            
            n += 1
            if len(vertices) == 0:
                return n
            v0 = vertices[0]
            v_conexos = []
            v_inconexos = []

    def existe_camino(self, v, u, aristas):
        #print(aristas)
        # Crear el grafo como un diccionario de adyacencia
        grafo = dict()
        
        # Construir el diccionario de adyacencia a partir de las aristas
        for v1, v2 in aristas:
            if v1 not in grafo.keys():
                grafo[v1] = []
            if v2 not in grafo:
                grafo[v2] = []
            grafo[v1].append(v2)
            grafo[v2].append(v1)  # Como es un grafo no dirigido, añadir ambas conexiones

        # Función recursiva para la búsqueda en profundidad (DFS)
        def bep(vertice, visitados):
            if vertice == u:
                return True
            visitados.append(vertice)  # Marcar el vértice como visitado
            # Recorrer los vecinos
            if vertice in grafo.keys():
                for vecino in grafo[vertice]:  # Obtener vecinos del vértice
                    if vecino not in visitados:  # Solo visitar los no visitados
                        if bep(vecino, visitados):
                            return True
            
            return False
        
        return bep(v, list())


com_sim = ComplejoSimplicial([(0,), (1,), (2,3), (4,5), (5,6), (4,6), (6,7,8,9)])
com_sim.print_caras()
print(com_sim.carac_euler())
print(com_sim.num_componentes_conexas())