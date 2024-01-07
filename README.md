![LBFGSB](media/logo.png)
============

This is a modern Fortran refactoring of the [L-BFGS-B](http://users.iems.northwestern.edu/~nocedal/lbfgsb.html) limited memory code for solving bound constrained optimization problems [L-BFGS-B version 3.0, march, 2011].

**This is a work in progress.**

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![GitHub release](https://img.shields.io/github/release/jacobwilliams/lbfgsb.svg?style=plastic)](https://github.com/jacobwilliams/lbfgsb/releases/latest)
[![Build Status](https://github.com/jacobwilliams/lbfgsb/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/lbfgsb/actions)
[![codecov](https://codecov.io/gh/jacobwilliams/lbfgsb/branch/master/graph/badge.svg?token=BHtd51oUTE)](https://codecov.io/gh/jacobwilliams/lbfgsb)

### Compiling

The library can be compiled with recent versions the Intel Fortran Compiler and GFortran (and presumably any other Fortran compiler that supports modern standards).

A `fpm.toml` file is provided for compiling `LBFGSB` with the [Fortran Package Manager](https://github.com/fortran-lang/fpm). For example, to build:

```
fpm build --profile release
```

By default, the library is built with double precision (`real64`) real values. Explicitly specifying the real kind can be done using the following processor flags:

Preprocessor flag | Kind  | Number of bytes
----------------- | ----- | ---------------
`REAL32`  | `real(kind=real32)`  | 4
`REAL64`  | `real(kind=real64)`  | 8
`REAL128` | `real(kind=real128)` | 16

For example, to build a single precision version of the library, use:

```
fpm build --profile release --flag "-DREAL32"
```

To run the unit tests:

```
fpm test --profile release
```

### Documentation

  * The latest API documentation can be found [here](https://jacobwilliams.github.io/lbfgsb/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford) (i.e. by running `ford ford.md`).

  * The following four articles were also included in the original package:

     1. `algorithm.pdf`     -  describes the algorithm
     2. `compact.pdf`       -  presents the compact form of matrices
     3. `code.pdf`          -  describes the code
     4. `acm-remark.pdf`    -  describes the modification/correction

   The most useful articles for those who only wish to use the code
   are [3,4]. Users interested in understanding the algorithm should
   read [1] and possibly also [2].

### How to use L-BFGS-B

The simplest way to use the code is to modify one of the
drivers provided in the package.  Most users will only need to make
a few changes to the drivers to run their applications.

The user is required to calculate the function value `f` and its gradient `g`.
In order to allow the user complete control over these computations,
reverse communication is used.  The routine `setulb` must be called
repeatedly under the control of the variable task.

### License

L-BFGS-B is released under the "New BSD License" (aka "Modified BSD License" or
"3-clause license")
Please read attached file License.txt

### References

* R. H. Byrd, P. Lu and J. Nocedal. [A Limited Memory Algorithm for Bound Constrained Optimization](http://www.ece.northwestern.edu/~nocedal/PSfiles/limited.ps.gz), (1995), SIAM Journal on Scientific and Statistical Computing , 16, 5, pp. 1190-1208.
* C. Zhu, R. H. Byrd and J. Nocedal. [L-BFGS-B: Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization](http://www.ece.northwestern.edu/~nocedal/PSfiles/lbfgsb.ps.gz) (1997), ACM Transactions on Mathematical Software, Vol 23, Num. 4, pp. 550 - 560.
* J.L. Morales and J. Nocedal. [L-BFGS-B: Remark on Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization](http://www.ece.northwestern.edu/~morales/PSfiles/acm-remark.pdf) (2011), ACM Transactions on Mathematical Software, Volume 38, Issue 1, No. 7, pp 1 - 4.

