#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 10/11/2021

"""

import pytest

from odelabs import BoundaryCondition


# ========== ========== ========== ========== ========== ========== Dirichlet
@pytest.fixture
def dir_hom_bc_0():
    """Dirichlet homogenous Boundary Condition at x=0"""
    return BoundaryCondition(0)


@pytest.fixture
def dir_nho_bc_0():
    """Dirichlet non-homogenous Boundary Condition at x=0"""
    return BoundaryCondition(0, c=1)


@pytest.fixture
def dir_hom_bc_1():
    """Dirichlet homogenous Boundary Condition at x=1"""
    return BoundaryCondition(1)


@pytest.fixture
def dir_nho_bc_1():
    """Dirichlet non-homogenous Boundary Condition at x=1"""
    return BoundaryCondition(1, c=1)


# ========== ========== ========== ========== ========== ========== Neumann
@pytest.fixture
def neu_hom_bc_0():
    """Neumann homogenous Boundary Condition at x=0"""
    return BoundaryCondition(0, a=0, b=1)


@pytest.fixture
def neu_nho_bc_0():
    """Neumann non-homogenous Boundary Condition at x=0"""
    return BoundaryCondition(0, a=0, b=1, c=1)


@pytest.fixture
def neu_hom_bc_1():
    """Neumann homogenous Boundary Condition at x=1"""
    return BoundaryCondition(1, a=0, b=1)


@pytest.fixture
def neu_nho_bc_1():
    """Neumann non-homogenous Boundary Condition at x=1"""
    return BoundaryCondition(1, a=0, b=1, c=1)


# ========== ========== ========== ========== ========== ========== mixed
@pytest.fixture
def mix_hom_bc_0():
    """Mixed homogenous Boundary Condition at x=0"""
    return BoundaryCondition(0, a=1, b=1)


@pytest.fixture
def mix_nho_bc_0():
    """Mixed non-homogenous Boundary Condition at x=0"""
    return BoundaryCondition(0, a=1, b=1, c=1)


@pytest.fixture
def mix_hom_bc_1():
    """Mixed homogenous Boundary Condition at x=1"""
    return BoundaryCondition(1, a=1, b=1)


@pytest.fixture
def mix_nho_bc_1():
    """Mixed non-homogenous Boundary Condition at x=1"""
    return BoundaryCondition(1, a=1, b=1, c=1)


# ========== ========== ========== ========== ========== ========== combinations
@pytest.fixture
def hom_bc_array(
        dir_hom_bc_0,
        dir_hom_bc_1,
        neu_hom_bc_0,
        neu_hom_bc_1,
        mix_hom_bc_0,
        mix_hom_bc_1):

    return [
        dir_hom_bc_0,
        dir_hom_bc_1,
        neu_hom_bc_0,
        neu_hom_bc_1,
        mix_hom_bc_0,
        mix_hom_bc_1
    ]


@pytest.fixture
def nho_bc_array(
        dir_nho_bc_0,
        dir_nho_bc_1,
        neu_nho_bc_0,
        neu_nho_bc_1,
        mix_nho_bc_0,
        mix_nho_bc_1):

    return [
        dir_nho_bc_0,
        dir_nho_bc_1,
        neu_nho_bc_0,
        neu_nho_bc_1,
        mix_nho_bc_0,
        mix_nho_bc_1
    ]


@pytest.fixture
def hom_bc_array_0(dir_hom_bc_0, neu_hom_bc_0, mix_hom_bc_0):
    return [dir_hom_bc_0, neu_hom_bc_0, mix_hom_bc_0]


@pytest.fixture
def hom_bc_array_1(dir_hom_bc_1, neu_hom_bc_1, mix_hom_bc_1):
    return [dir_hom_bc_1, neu_hom_bc_1, mix_hom_bc_1]


@pytest.fixture
def nho_bc_array_0(dir_nho_bc_0, neu_nho_bc_0, mix_nho_bc_0):
    return [dir_nho_bc_0, neu_nho_bc_0, mix_nho_bc_0]


@pytest.fixture
def nho_bc_array_1(dir_nho_bc_1, neu_nho_bc_1, mix_nho_bc_1):
    return [dir_nho_bc_1, neu_nho_bc_1, mix_nho_bc_1]


# ========== ========== ========== ========== ========== ==========
def _test_nonhomogeneous_polynomial_fit(lbc: BoundaryCondition, ubc: BoundaryCondition):
    poly = BoundaryCondition.fit_polynomial(lbc, ubc)

    print()
    print(lbc, ubc, poly, sep=', ')
    print(lbc.get_polynomial_error(poly))
    print(ubc.get_polynomial_error(poly))
    print('\n\n')

    assert lbc.is_satisfied_by_polynomial(poly)
    assert ubc.is_satisfied_by_polynomial(poly)


def test_nonhomogeneous_polynomial_fit(hom_bc_array_0, hom_bc_array_1, nho_bc_array_0, nho_bc_array_1):

    for lbc in hom_bc_array_0:
        for ubc in nho_bc_array_1:
            _test_nonhomogeneous_polynomial_fit(lbc, ubc)

    for lbc in nho_bc_array_0:
        for ubc in hom_bc_array_1:
            _test_nonhomogeneous_polynomial_fit(lbc, ubc)

    for lbc in nho_bc_array_0:
        for ubc in nho_bc_array_1:
            _test_nonhomogeneous_polynomial_fit(lbc, ubc)

    # ---------- ---------- ---------- ---------- ---------- ----------
    lbc = BoundaryCondition(-1, a=1, b=1, c=1)
    ubc = BoundaryCondition(1, a=1, b=-1, c=0)

    _test_nonhomogeneous_polynomial_fit(lbc, ubc)

    # ---------- ---------- ---------- ---------- ---------- ----------
    lbc = BoundaryCondition(-1, a=1, b=1, c=0)
    ubc = BoundaryCondition(1, a=1, b=-1, c=1)

    _test_nonhomogeneous_polynomial_fit(lbc, ubc)


# ========== ========== ========== ========== ========== ==========
def _test_homogeneous_polynomial_fit(lbc: BoundaryCondition, ubc: BoundaryCondition, degree):
    print()
    print(lbc, ubc, sep=', ')
    poly = BoundaryCondition.fit_polynomial(lbc, ubc, min_degree=degree)
    print(poly)
    print(lbc.get_polynomial_error(poly))
    print(ubc.get_polynomial_error(poly))
    print('\n\n')

    assert lbc.is_satisfied_by_polynomial(poly)
    assert ubc.is_satisfied_by_polynomial(poly)


def test_homogeneous_polynomial_fit(hom_bc_array_0, hom_bc_array_1):

    for lbc in hom_bc_array_0:
        for ubc in hom_bc_array_1:
            _test_homogeneous_polynomial_fit(lbc, ubc, degree=2)
            _test_homogeneous_polynomial_fit(lbc, ubc, degree=3)
            _test_homogeneous_polynomial_fit(lbc, ubc, degree=4)
            _test_homogeneous_polynomial_fit(lbc, ubc, degree=5)
            _test_homogeneous_polynomial_fit(lbc, ubc, degree=6)


































