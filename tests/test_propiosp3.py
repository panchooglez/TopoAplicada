import pytest
import numpy as np
from ejs.p3 import ComplejoSimplicial

def test_matriz_borde_dim_dos():
    simplices_maximales = [(0, 1, 2)]
    complejo = ComplejoSimplicial(simplices_maximales)
    matriz_esperada = np.array([
        [1],
        [1],
        [1]
    ])
    matriz_obtenida = complejo.matriz_borde(2)
    assert np.array_equal(matriz_obtenida, matriz_esperada)

def test_rango_z2():
    matriz = np.array([
        [1, 0, 1],
        [1, 1, 0],
        [0, 1, 1]
    ])
    complejo = ComplejoSimplicial([])
    rango = complejo.rango_z2(matriz)
    assert rango == 2  # El rango esperado es 2 para esta matriz

def test_calculo_betti_simple():
    simplices_maximales = [(0, 1, 2)]
    complejo = ComplejoSimplicial(simplices_maximales)
    betti = complejo.calcular_betti()
    assert betti == [1, 0, 0]  # El triángulo relleno no tiene agujeros

def test_calculo_betti_multiple():
    simplices_maximales = [(0, 1), (1, 2), (2, 0)]
    complejo = ComplejoSimplicial(simplices_maximales)
    betti = complejo.calcular_betti()
    assert betti == [1, 1]  # Un ciclo, debería tener b_1 = 1

def test_calcular_betti_incremental_simple():
    simplices_maximales = [(0, 1), (1, 2), (2, 0)]
    complejo = ComplejoSimplicial(simplices_maximales)
    betti = complejo.calcular_betti_incremental()
    assert betti == [1, 1]  # Un ciclo, b_0 = 1, b_1 = 1

def test_calculo_betti_alfa_complejo():
    puntos = [(0, 0), (1, 0), (0, 1), (1, 1)]
    radio = 1.5
    complejo = ComplejoSimplicial([])
    betti = complejo.calcular_betti_alfa_complejos(puntos, radio)
    assert betti == [1, 1]  # Un cuadrado con 1 ciclo

def test_calculo_betti_alfa_complejo_con_ciclo():
    puntos = [(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
    radio = 1.5
    complejo = ComplejoSimplicial([])
    betti = complejo.calcular_betti_alfa_complejos(puntos, radio)
    assert betti == [1, 0]  # Complejo lleno, sin ciclos

def test_dos_ciclos():
    simplices_maximales = [(0, 1, 2), (3, 4, 5)]
    complejo = ComplejoSimplicial(simplices_maximales)
    betti = complejo.calcular_betti()
    assert betti == [2, 0, 0]  # Dos componentes conexas, sin agujeros

def test_sin_simplices():
    complejo = ComplejoSimplicial([])
    betti = complejo.calcular_betti()
    assert betti == [0]  # Sin simplices, sin homología

def test_un_punto():
    simplices_maximales = [(0,)]
    complejo = ComplejoSimplicial(simplices_maximales)
    betti = complejo.calcular_betti()
    assert betti == [1]  # Un componente conexa

def test_matriz_borde_dim_uno():
    simplices_maximales = [(0, 1)]
    complejo = ComplejoSimplicial(simplices_maximales)
    matriz_esperada = np.array([
        [1],
        [1]
    ])
    matriz_obtenida = complejo.matriz_borde(1)
    assert np.array_equal(matriz_obtenida, matriz_esperada)
