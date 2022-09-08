# The enriched subspace iteration method

## Introduction

The enriched subspace iteration method is to find the smallest $p$ eigenvalues
($\lambda_i; i = 1, \ldots, p$) and the corresponding eigenvectors
($\boldsymbol{\phi}_i; i = 1, \ldots, p$) that satisfy
$$
\mathbf{K} \boldsymbol{\phi}_i = \lambda_i \mathbf{M} \boldsymbol{\phi}_i; 
\quad i = 1, \ldots, p
$$
and
$$
\begin{align*}
\boldsymbol{\phi}_i^T \mathbf{M} \boldsymbol{\phi}_j &= \delta_{ij}, \\
\boldsymbol{\phi}_i^T \mathbf{K} \boldsymbol{\phi}_j &= \lambda_i \delta_{ij}
\end{align*}
$$
where $\mathbf{K} \in \mathbb{R}^{n \times n}$ and $\mathbf{M} \in
\mathbb{R}^{n \times n}$ are real symmetric sparse matrices (usually obtained
after discretization like finite elements), and $\delta_{ij}$ is the Kronecker
delta.

<!-- For a detailed description of the enriched subspace iteration method, see [^1]. -->

<!-- The generalized eigenvalue problem defined above arises in many application areas. -->

## Installation

#### Requirements

To build the enriched subspace iteration library you need to have:
- a Fortran 2008 compliant compiler or newer (GCC Fortran or Intel Fortran compilers)
- CMake version 3.13+

#### Installation on Linux

Clone this repository

```bash
git clone https://github.com/ktkimit/esubspace.git
cd esubspace
```

Configure the build and build it with
```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=/to/your/preferred/install/location
cmake --build build
```

Install the library with
```bash
cmake --install build
```

Then the library will be installed into the location you specified with
`DCMAKE_INSTALL_PREFIX` option in the configuration step.

## How to use the library in your project?

The library exports CMake package files to easily find and use it in other
projects.
The CMake files are located in the library directory in the installation
directory you specified with the installation prefix.

You can find and use the installed enriched subspace iteration library by
simply put the following in your CMake configuration
```bash
find_package(Esspace REQUIRED)
...

target_link_libraries(${PROJECT_NAME} PRIVATE Esspace::Esspace)
```

When configure your project, don't forget to make the installed enriched
subspace iteration library discoverable by adding its CMake package directory
to `CMAKE_PREFIX_PATH`.

For example, please look at the CMake configuration of the example given in
this repository. To build the example type the following at its directory
(./example/frame2d/)
```bash
cmake -B build -DCMAKE_PREFIX_PATH=/to/your/preferred/install/location/lib/cmake
```
