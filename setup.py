#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 10/11/2021

"""

from setuptools import setup

import odelabs


if __name__ == '__main__':
    setup(
        name='odelabs',
        license='MIT',
        packages=[
            'odelabs'
        ],
        url=odelabs.__url__,
        version=odelabs.__version__,
        author=odelabs.__author__,
        author_email=odelabs.__email__,
        description=odelabs.__description__,
        long_description=odelabs.__doc__,
        classifiers=[
            'Development Status :: 1 - Planning',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: 3.11',
        ],
        install_requires=[
            'numpy',
            'matplotlib',
            'scipy',
            'sympy'
        ],
        python_requires='>=3.9'
    )

