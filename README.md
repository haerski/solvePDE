# Linear PDE with Constant Coefficients
This repository contains the code accompanying the paper _Linear PDE with constant coefficients_ by Rida Ait El Manssour, Marc Härkönen, and Bernd Sturmfels.

**Note**: These functions will added to the package `NotherianOperators` and distributed with Macaulay2 from version 1.18 (ETA May 2021). 

## Quickstart
The functions in `solvePDE.m2` require Macaulay2 version 1.17 or newer.

Clone the repository, and change to the newly created directory called `solvePDE`
```bash
git clone https://github.com/haerski/solvePDE.git
cd solvePDE
```
The following command will open Macaulay2 and preload the functions
```
M2 solvePDE.m2
```

Alternatively, the functions can be loaded in an existing Macaulay2 session using the command `load("/path/to/solvePDE/solvePDE.m2")`

## A typical session
The main function is called `solvePDE`. Accepted inputs are ideals or modules. The output is a list of pairs whose first entry is a prime and second entry is a list of Noetherian multipliers. The total number of Noetherian multiplier will be equal to the arithmetic multiplicity, which can be computed using the command `amult`.

**Note**: The ouput can be fed to the command `netList` to make the output more human-readable.

### Example 1: ideals
```macaulay2
$ M2 solvePDE.m2
Macaulay2, version 1.17.2.1
with packages: ConwayPolynomials, Elimination, IntegralClosure, InverseSystems, LLLBases, MinimalPrimes, PrimaryDecomposition, ReesAlgebra, Saturation, TangentCone

i1 : R = QQ[x,y,z]

o1 = R

o1 : PolynomialRing

i2 : I = ideal(x^2*y,x^2*z,x*y^2,x*y*z^2)

             2    2      2       2
o2 = ideal (x y, x z, x*y , x*y*z )

o2 : Ideal of R

i3 : amult I

o3 = 5

o3 : QQ

i4 : solvePDE I

o4 = {{ideal x, {| 1 |}}, {ideal (y, x), {| dx |}}, {ideal (z, y), {| 1 |}}, {ideal (z, y, x), {| dxdy |, | dxdydz |}}}

o4 : List

i5 : netList solvePDE I

     +---------------+----------------------+
o5 = |ideal x        |{| 1 |}               |
     +---------------+----------------------+
     |ideal (y, x)   |{| dx |}              |
     +---------------+----------------------+
     |ideal (z, y)   |{| 1 |}               |
     +---------------+----------------------+
     |ideal (z, y, x)|{| dxdy |, | dxdydz |}|
     +---------------+----------------------+
```

### Example 2: modules
```macaulay2
$ M2 solvePDE.m2
Macaulay2, version 1.17.2.1
with packages: ConwayPolynomials, Elimination, IntegralClosure, InverseSystems, LLLBases, MinimalPrimes, PrimaryDecomposition, ReesAlgebra, Saturation, TangentCone

i1 : R = QQ[x,y,z];

i2 : U = image(matrix {
         {x^2,x*y,x*z},
         {y^2,y*z,z^2}});

i3 : netList solvePDE U

     +---------------+--------+
o3 = |ideal x        |{| 1 |} |
     |               | | 0 |  |
     +---------------+--------+
     |       2       |        |
     |ideal(y  - x*z)|{| -z |}|
     |               | |  x | |
     +---------------+--------+
     |ideal (z, y)   |{|  0 |}|
     |               | | dz | |
     +---------------+--------+
```
