import numpy as np
import pytest
from ejs.p3 import ComplejoSimplicial  # Cambia "your_module" al nombre del archivo principal

def assert_matrix_equal(mat1, mat2):
    """
    Compara dos matrices numpy, considerando valores en Z2.
    """
    assert np.array_equal(mat1 % 2, mat2 % 2), f"{mat1} != {mat2} en Z2"

def test_ejemplo_1():
    sc = ComplejoSimplicial([(0, 1, 2, 3)])
    expected_b1 = np.array([
        [1, 1, 1, 0, 0, 0],
        [1, 0, 0, 1, 1, 0],
        [0, 1, 0, 1, 0, 1],
        [0, 0, 1, 0, 1, 1],
    ])
    expected_b2 = np.array([
        [1, 1, 0, 0],
        [1, 0, 1, 0],
        [0, 1, 1, 0],
        [1, 0, 0, 1],
        [0, 1, 0, 1],
        [0, 0, 1, 1],
    ])
    expected_b3 = np.array([[1], [1], [1], [1]])
    expected_betti = [1, 0, 0, 0]

    assert_matrix_equal(sc.matriz_borde(1), expected_b1)
    assert_matrix_equal(sc.matriz_borde(2), expected_b2)
    assert_matrix_equal(sc.matriz_borde(3), expected_b3)
    assert sc.calcular_betti() == expected_betti

def test_ejemplo_2():
    sc = ComplejoSimplicial([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)])
    expected_betti = [1, 0, 1]
    assert sc.calcular_betti() == expected_betti

def test_ejemplo_3():
    sc = ComplejoSimplicial([
        (0, 1),
        (1, 2, 3, 4),
        (4, 5),
        (5, 6),
        (4, 6),
        (6, 7, 8),
        (8, 9)
    ])

    # Números de Betti esperados
    expected_betti = [1, 1, 0, 0]

    # Validación de números de Betti
    assert sc.calcular_betti() == expected_betti


def test_ejemplo_5():
    sc = ComplejoSimplicial([(0, 1, 2), (2, 3), (3, 4)])
    expected_betti = [1, 0, 0]
    assert sc.calcular_betti() == expected_betti

def test_ejemplo_6():
    sc = ComplejoSimplicial([
        (1, 2, 4), (1, 3, 6), (1, 4, 6),
        (2, 3, 5), (2, 4, 5), (3, 5, 6)
    ])
    expected_betti = [1, 1, 0]
    assert sc.calcular_betti() == expected_betti


def test_ejemplo_8():
    sc = ComplejoSimplicial([
        (1, 2, 4), (2, 4, 5), (2, 3, 5), (3, 5, 6),
        (1, 3, 6), (1, 4, 6), (4, 5, 7), (5, 7, 8),
        (5, 6, 8), (6, 8, 9), (4, 6, 9), (4, 7, 9),
        (1, 7, 8), (1, 2, 8), (2, 8, 9), (2, 3, 9),
        (3, 7, 9), (1, 3, 7)
    ])
    expected_betti = [1, 2, 1]
    assert sc.calcular_betti() == expected_betti

def test_ejemplo_9():
    sc = ComplejoSimplicial([
        (1, 2, 6), (2, 3, 4), (1, 3, 4), (1, 2, 5), (2, 3, 5),
        (1, 3, 6), (2, 4, 6), (1, 4, 5), (3, 5, 6), (4, 5, 6)
    ])
    expected_betti = [1, 1, 1]
    assert sc.calcular_betti() == expected_betti

def test_ejemplo_10():
    sc = ComplejoSimplicial([(0,), (1,), (2, 3), (4, 5), (5, 6), (4, 6), (6, 7, 8, 9)])
    expected_betti = [4, 1, 0, 0]
    assert sc.calcular_betti() == expected_betti
