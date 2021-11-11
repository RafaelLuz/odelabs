#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 10/11/2021

"""

import numpy

from odelabs.utils import parse_float, parse_positive_integer, parse_non_negative_integer


class BoundaryCondition:
    """
    Abstraction for conditions of the form

    a*y(x) + b*y'(x) = c

    where a, b and c are real numbers

    """
    # ========== ========== ========== ========== ========== class attributes
    ...

    # ========== ========== ========== ========== ========== special methods
    def __init__(self, x, a=1, b=0, c=0, theta=None):

        self.__x = parse_float(x)
        self._parse_theta_and_c(a, b, c, theta)

    def __copy__(self):
        return type(self)(self.x, theta=self.theta, c=self.c)

    def __repr__(self):

        if self.sin_theta == 0:
            return f"Boundary Condition: y({self.x:.4g}) = {self.c:.4g}"

        if self.sin_theta == 1:
            return f"Boundary Condition: y'({self.x:.4g}) = {self.c:.4g}"

        if self.sin_theta > 0:
            return f"Boundary Condition: {self.cos_theta:.4g}*y({self.x:.4g}) + {self.sin_theta:.4g}*y'({self.x:.4g}) = {self.c:.4g}"

        return f"Boundary Condition: {self.cos_theta:.4g}*y({self.x:.4g}) - {-self.sin_theta:.4g}*y'({self.x:.4g}) = {self.c:.4g}"

    # ========== ========== ========== ========== ========== private methods
    @classmethod
    def __fit_polynomial_for_nonhomogeneous_bcs(cls, lbc, ubc):
        """
        Ac = B
        :param lbc:
        :param ubc:
        :param degree:
        :return:
        """
        B = numpy.array([lbc.c, ubc.c]).reshape(-1, 1)

        A = numpy.array([
            [lbc._polynomial_regressor(0), lbc._polynomial_regressor(1)],
            [ubc._polynomial_regressor(0), ubc._polynomial_regressor(1)],
        ])

        if not cls._equals_zero(numpy.linalg.det(A)):
            # hence, a polynomial of degree 1 is enough
            coefficients = numpy.linalg.solve(A, B).reshape(-1,)

        elif not (cls._equals_zero(lbc.theta + ubc.theta) and cls._equals_zero(lbc.tan_theta - (ubc.x - lbc.x) / 2)):
            # hence, a polynomial of degree 2 is enough
            col = numpy.array([
                lbc._polynomial_regressor(2),
                ubc._polynomial_regressor(2)
            ]).reshape(-1, 1)

            A = numpy.hstack([A, col])

            coefficients = numpy.linalg.lstsq(A, B)[0].reshape(-1, )

        else:
            # hence, a polynomial of degree 3 is necessary

            cols = numpy.array([
                [lbc._polynomial_regressor(2), lbc._polynomial_regressor(3)],
                [ubc._polynomial_regressor(2), ubc._polynomial_regressor(3)]
            ])

            A = numpy.hstack([A, cols])

            coefficients = numpy.linalg.lstsq(A, B)[0].reshape(-1, )

        return numpy.polynomial.Polynomial(coefficients)

    @classmethod
    def __fit_polynomial_for_homogeneous_bcs(cls, lbc, ubc, degree):

        assert degree >= 2

        M = numpy.array([
            [lbc._polynomial_regressor(0), lbc._polynomial_regressor(1), lbc._polynomial_regressor(degree)],
            [ubc._polynomial_regressor(0), ubc._polynomial_regressor(1), ubc._polynomial_regressor(degree)],
        ])

        A = M[:, [True, True, False]]
        dA = numpy.linalg.det(A)

        B = M[:, [True, False, True]]
        dB = numpy.linalg.det(B)

        C = M[:, [False, True, True]]
        dC = numpy.linalg.det(C)

        if not cls._equals_zero(dA):
            # no problem! planes are not the same, but c2 is not zero
            coefficients = [dC/dA, -dB/dA] + (degree-2)*[0] + [1]

        else:

            if cls._equals_zero(dB) and cls._equals_zero(dC):
                # no problem! planes are the same

                if cls._equals_zero(M[0, 0]):
                    coefficients = [0, -M[0, 2] / M[0, 1]] + (degree - 2) * [0] + [1]

                else:
                    coefficients = [-M[0, 2] / M[0, 0], 0] + (degree - 2) * [0] + [1]

            else:
                # problem!! planes are not the same, but c2 is zero

                if degree == 2:  # fit a polynomial of degree < 2
                    coefficients = numpy.linalg.lstsq(A[1, :].reshape(1, -1), [[0]])[0].reshape(-1, )

                else:
                    return cls.__fit_polynomial_for_homogeneous_bcs_alternative(lbc, ubc, degree)

        return numpy.polynomial.Polynomial(coefficients)

    @classmethod
    def __fit_polynomial_for_homogeneous_bcs_alternative(cls, lbc, ubc, degree):

        assert degree > 2

        M = numpy.array([
            [lbc._polynomial_regressor(0), lbc._polynomial_regressor(2), lbc._polynomial_regressor(degree)],
            [ubc._polynomial_regressor(0), ubc._polynomial_regressor(2), ubc._polynomial_regressor(degree)],
        ])

        A = M[:, [True, True, False]]
        dA = numpy.linalg.det(A)

        B = M[:, [True, False, True]]
        dB = numpy.linalg.det(B)

        C = M[:, [False, True, True]]
        dC = numpy.linalg.det(C)

        if not cls._equals_zero(dA):
            # no problem! planes are not the same, but leading coefficient is not zero
            coefficients = [dC / dA, -dB / dA] + (degree - 2) * [0] + [1]

        else:

            if cls._equals_zero(dB) and cls._equals_zero(dC):
                # no problem! planes are the same

                if cls._equals_zero(M[0, 0]):
                    coefficients = [0, 0, -M[0, 2] / M[0, 1]] + (degree - 1) * [0] + [1]

                else:
                    coefficients = [-M[0, 2] / M[0, 0], 0, 0] + (degree - 1) * [0] + [1]

            else:
                # problem!! planes are not the same, but leading coefficient is zero
                raise NotImplementedError()

        return numpy.polynomial.Polynomial(coefficients)

    # ========== ========== ========== ========== ========== protected methods
    def _parse_theta_and_c(self, a, b, c, theta):

        if theta is not None:
            self.__theta = parse_float(theta)
            assert -numpy.pi/2 < self.__theta <= numpy.pi/2, f"Expected -pi/2 < theta <= pi/2. Given {theta}"
            self.__c = parse_float(c)

        else:
            a, b, c = parse_float(a), parse_float(b), parse_float(c)

            if a == 0:
                if b == 0:
                    error = "'a' and 'b' can not be both zero"
                    raise ValueError(error)

                self.__theta = numpy.pi / 2
                self.__c = c / b

            else:
                self.__theta = numpy.arctan(b / a)
                self.__c = (numpy.cos(self.__theta) / a) * c

    def _polynomial_regressor(self, degree):

        degree = parse_non_negative_integer(degree)

        if degree == 0:
            return self.cos_theta

        return self.x**(degree-1) * (self.x*self.cos_theta + degree * self.sin_theta)

    @classmethod
    def _equals_zero(cls, value):
        return numpy.abs(value) < 1e-15

    # ========== ========== ========== ========== ========== public methods
    def is_homogeneous(self):
        return self.c == 0

    def is_dirichlet(self):
        return self.theta == 0

    def is_neumann(self):
        return self.theta == numpy.pi/2

    def is_mixed(self):
        return not self.is_dirichlet() and not self.is_neumann()

    def get_polynomial_error(self, poly):
        return self.cos_theta * poly(self.x) + self.sin_theta * poly.deriv()(self.x) - self.c

    def is_satisfied_by_polynomial(self, poly):
        return self._equals_zero(self.get_polynomial_error(poly))

    def get_homogenous_copy(self):
        return type(self)(self.x, theta=self.theta, c=0)

    def get_nonhomogeneous_copy(self, c):
        return type(self)(self.x, theta=self.theta, c=parse_float(c))

    @classmethod
    def fit_polynomial(cls, lbc, ubc, degree=None):
        """

        :param lbc:
        :param ubc:
        :return:
        """
        assert isinstance(lbc, BoundaryCondition) and isinstance(ubc, BoundaryCondition), \
            "lbc and ubc must be instances of BoundaryCondition"

        assert lbc.x < ubc.x, f"Expected lbc.x < ubc.x but {lbc.x} >= {ubc.x}"

        # ---------- ---------- ---------- ---------- ---------- homogeneous bc
        if lbc.is_homogeneous() and ubc.is_homogeneous():
            return cls.__fit_polynomial_for_homogeneous_bcs(lbc, ubc, degree)

        # ---------- ---------- ---------- ---------- ---------- non-homogeneous
        else:
            return cls.__fit_polynomial_for_nonhomogeneous_bcs(lbc, ubc)

    # ---------- ---------- ---------- ---------- ---------- properties
    @property
    def x(self):
        return self.__x

    @property
    def theta(self):
        return self.__theta

    @property
    def c(self):
        return self.__c

    @property
    def cos_theta(self):

        try:
            return self.__cos_theta

        except AttributeError:
            self.__cos_theta = 0.0 if self.theta == numpy.pi/2 else numpy.cos(self.theta)

            return self.__cos_theta

    @property
    def sin_theta(self):

        try:
            return self.__sin_theta

        except AttributeError:
            self.__sin_theta = numpy.sin(self.theta)

            return self.__sin_theta

    @property
    def tan_theta(self):

        try:
            return self.__tan_theta

        except AttributeError:
            self.__tan_theta = numpy.inf if self.theta == numpy.pi/2 else numpy.tan(self.theta)

            return self.__tan_theta
