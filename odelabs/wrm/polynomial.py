#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 11/11/2021

"""

import numpy

from numpy.polynomial import Polynomial

from odelabs import BoundaryCondition, Domain


class PolynomialCoefficientsOperator:

    # ========== ========== ========== ========== ========== class attributes
    ...

    # ========== ========== ========== ========== ========== special methods
    def __init__(self, p, r, q):
        self.__p = Polynomial(p)
        self.__r = Polynomial(r)
        self.__q = Polynomial(q)

    def __call__(self, poly):
        return self.p * poly.deriv(2) + self.r * poly.deriv(1) + self.q * poly

    # ========== ========== ========== ========== ========== private methods
    ...

    # ========== ========== ========== ========== ========== protected methods
    ...

    # ========== ========== ========== ========== ========== public methods
    ...

    # ---------- ---------- ---------- ---------- ---------- properties
    @property
    def p(self):
        return self.__p

    @property
    def q(self):
        return self.__q

    @property
    def r(self):
        return self.__r


class PolynomialNonHomogeneousBVP:

    # ========== ========== ========== ========== ========== class attributes
    ...

    # ========== ========== ========== ========== ========== special methods
    def __init__(self, domain=(0, 1), lbc=(1, 0, 0), ubc=(1, 0, 0), p=(1,), r=(0,), q=(1,), u=(1,)):

        self.__domain = Domain(*domain)
        self.__operator = PolynomialCoefficientsOperator(p, r, q)
        self.__lbc = BoundaryCondition(self.domain.infimum, *lbc)
        self.__ubc = BoundaryCondition(self.domain.supremum, *ubc)
        self.__input = Polynomial(u)

    # ========== ========== ========== ========== ========== private methods
    ...

    # ========== ========== ========== ========== ========== protected methods
    ...

    # ========== ========== ========== ========== ========== public methods
    def inner(self, p1, p2):
        # return (p1*p2).integ(lbnd=0)(1)
        return self.domain.integrate_polynomial(p1*p2)

    def solve_galerkin(self, n=3):

        if self.lbc.is_homogeneous() and self.ubc.is_homogeneous():
            phi_nh = Polynomial((0,))

        else:
            phi_nh = BoundaryCondition.fit_polynomial(self.lbc, self.ubc)

        # print(phi_nh)

        lbc = self.lbc.get_homogenous_copy()
        ubc = self.ubc.get_homogenous_copy()

        phis = []
        for degree in range(2, n+2):
            phis.append(BoundaryCondition.fit_polynomial(lbc, ubc, min_degree=degree))

        # for phi in phis:
        #     print(phi)

        Ann = [[self.inner(w, self.operator(phi)) for phi in phis] for w in phis]
        Ann = numpy.array(Ann)

        bn = [self.inner(w, self.input - self.operator(phi_nh)) for w in phis]
        bn = numpy.array(bn)

        # print(Ann)
        # print(bn)

        cn = numpy.linalg.solve(Ann, bn)

        s = 0
        for c, phi in zip(cn, phis):
            s += c*phi

        return phi_nh + s

    # ---------- ---------- ---------- ---------- ---------- properties
    @property
    def domain(self):
        return self.__domain

    @property
    def operator(self):
        return self.__operator

    @property
    def lbc(self):
        return self.__lbc

    @property
    def ubc(self):
        return self.__ubc

    @property
    def input(self):
        return self.__input


