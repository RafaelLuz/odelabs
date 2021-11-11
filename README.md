# odelabs - Ordinary Differential Equations Labs


This package intends to provide tools for analysing and solving second
order linear ordinary differential equations.

----
## Table of Contents

* [Installation](#installation)
* [Usage](#usage)
* [License](#license)


----
## Installation

```shell script
pip install odelabs
```


----
## Usage

Consider the second order linear differential operator given by 

![equation](https://latex.codecogs.com/svg.latex?\large&space;\mathsf{L}y\left(x\right)=p\left(x\right)y^{\prime\prime}\left(x\right)+r\left(x\right)y^\prime\left(x\right)+q\left(x\right)y\left(x\right)),

where the functions are real and defined at the compact interval
![equation](https://latex.codecogs.com/svg.latex?\large&space;\left[x^-,x^+\right]), and mixed
boundary conditions of the form

![equation](https://latex.codecogs.com/svg.latex?\large&space;a^-y\left(x^-\right)+b^-y^{\prime}\left(x^-\right)=c^-) 

![equation](https://latex.codecogs.com/svg.latex?\large&space;a^+y\left(x^+\right)+b^+y^{\prime}\left(x^+\right)=c^+) 

obeying ![equation](https://latex.codecogs.com/svg.latex?\large&space;\left(a^{\pm}\right)^2+\left(b^{\pm}\right)^2\neq0).

### Boundary Conditions

In general case, 

![equation](https://latex.codecogs.com/svg.latex?\large&space;ay\left(x\right)+by^{\prime}\left(x\right)=c) 


```pycon
>>> from odelabs import BoundaryCondition as BC

>>> lbc = BC(x=0, a=0, b=1, c=3)  # Nonhomogeneous Neumann BC at x=0
>>> lbc
Boundary Condition: y'(0) = 3

>>> ubc = BC(x=1, a=1, b=0, c=0)  # Homogeneous Dirichlet BC at x=1
>>> ubc
Boundary Condition: y(1) = 0

>>> poly = BC.fit_polynomial(lbc, ubc)
>>> poly
-3.0 + 3.0·x¹

```

[comment]: <> (<img src="https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />)


[comment]: <> (![\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}]&#40;https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}&#41; )

----
## License

[The MIT License](LICENSE)