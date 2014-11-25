\name{bigalgebra-package}
\alias{bigalgebra-package}
\alias{bigalgebra}
\docType{package}
\title{
Arithmetic routines for native R matrices and big.matrix objects.
}
\description{
This package provides arithmetic functions for native \R matrices and
\code{\link[bigmemory]{big.matrix}} objects.
}
\details{
This package provides arithmetic functions for native \R matrices and
\code{\link[bigmemory]{big.matrix}} objects. 

The package defines a number of 
global options that begin with \code{bigalgebra}.
They include:
\tabular{ll}{
Option \tab Default value\cr
\code{bigalgebra.temp_pattern} \tab "matrix_"\cr
\code{bigalgebra.tempdir} \tab \code{tempdir}\cr
\code{bigalgebra.mixed_arithmetic_returns_R_matrix} \tab TRUE\cr
\code{bigalgebra.DEBUG} \tab FALSE\cr
}
The \code{bigalgebra.tempdir} option must be a function that returns
a temporary directory path used to big matrix results of BLAS and
LAPACK operations. The deault value is simply the default R \code{tempdir}
function.

The \code{bigalgebra.temp_pattern} is a name prefix for file names of generated
big matrix objects output as a result of BLAS and LAPACK operations.

The \code{bigalgebra.mixed_arithmetic_returns_R_matrix} option determines
whether arithmetic operations involving an R matrix or vector and a big.matrix
matrix or vector return a big matrix (when the option is \code{FALSE}), or
return a normal R matrix (\code{TRUE}).

The package is built, by default, with \R's native BLAS libraries, which use
32-bit signed integer indexing. The default build is limited to vectors of at
most 2**31 - 1 entries and matrices with at most 2**31 - 1 rows and 2**31 - 1
columns (note that standard R matrices are limtied to 2**31 - 1 total entries).

The package includes a reference BLAS implementation that supports 64-bit
integer indexing, relaxing the limitation on vector lengths and matrix
row and column limits. Installation of this package with the 64-bit reference
BLAS implementation may be performed from the command-line install:

\code{REFBLAS=1 R CMD INSTALL bigalgebra}

where "\code{bigalgebra}" is the source package (for example,
\code{bigalgebra_0.8.4.tar.gz}).

The package may also be build with user-supplied external BLAS and LAPACK
libraries, in either 32- or 64-bit varieties. This is an advanced topic
that requires additional Makevars modification, and may include adjustment
of the low-level calling syntax depending on the library used.

Feel free to contact us for help installing and running the package.
}
\author{
Michael J. Kane, Bryan Lewis, John W. Emerson

Maintainers: Mike and Jay <bigmemoryauthors@gmail.com>
}
\references{
\url{http://www.netlib.org/blas/}

\url{http://www.netlib.org/lapack/}
}
\keyword{package}
\seealso{
\code{\link[bigmemory]{bigmemory}}, \code{\link[bigmemory]{big.matrix}} 
}
\examples{

# Testing the development of the user-friendly operators:
# if you have any problems, please email us!  - Jay & Mike 4/29/2010

library("bigmemory")
A <- big.matrix(5,4,init=1)
B <- big.matrix(4,4,init=2)

C <- A %*% B       # Returns a new big.matrix object
D <- A[] %*% B[]   # Compute the same thing in R

print(C - D)       # Compare the results (subtraction of an R matrix from a
                   # big.matrix)

# The next example illustrates mixing R and big.matrix objects. It returns by
# default (see # options("bigalgebra.mixed_arithmetic_returns_R_matrix")
D <- matrix(rnorm(16),4)
E <- A %*% D
}
