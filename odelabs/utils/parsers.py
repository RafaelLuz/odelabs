#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 10/11/2021

"""


def parse_integer(value):
    try:
        return int(value)

    except TypeError or ValueError:
        error = f"Expected integer number (float). Given {value}"
        raise ValueError(error)


def parse_positive_integer(value):
    value = parse_integer(value)

    assert value > 0, f"Expected positive integer. Given {value}"

    return value


def parse_non_negative_integer(value):
    value = parse_integer(value)

    assert value >= 0, f"Expected positive integer. Given {value}"

    return value


def parse_float(value):

    try:
        return float(value)

    except TypeError or ValueError:
        error = f"Expected real number (float). Given {value.__class__.__name__}"
        raise ValueError(error)
