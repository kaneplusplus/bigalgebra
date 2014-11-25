\name{daxpy}
\Rdversion{1.1}
\alias{daxpy}
\title{BLAS daxpy functionality}
\description{
This function implements the function Y := A * X + Y where X and Y may be
either native double-precision valued R matrices or numeric vectors, or
double-precision valued \code{\link[bigmemory]{big.matrix}} objects, and A is a
scalar.
}
\usage{
daxpy(A=1, X, Y)
}
\arguments{
  \item{A}{Optional numeric scalar value to scale the matrix \code{X} by,
           with a default value of 1.}
  \item{X}{Requried to be either a native \R \code{\link{matrix}}
           or numeric vector,
           or a \code{\link[bigmemory]{big.matrix}} object}
  \item{Y}{Optional native \R \code{\link{matrix}} or numeric vector,
           or a \code{\link[bigmemory]{big.matrix}} object}
}
\details{
At least one of either \code{X} or \code{Y} must be a \code{big.matrix}.
All values must be of type \code{double} (the only type presently supported
by the bigalgebra package).

This function is rarely necessary to use directly since the bigalgebra
package defines standard arithmetic operations and scalar multiplication.
It is more efficient to use \code{daxpy} directly when both scaling and
matrix addition are required, in which case both operations are performed
in one step.
}
\value{
The output value depends on the classes of input values \code{X} and \code{Y}
and on the value of the global option
\code{bigalgebra.mixed_arithmetic_returns_R_matrix}.

If \code{X} and \code{Y} are both big matrices, or \code{Y} is missing,
\code{options("bigalgebra.mixed_arithmetic_returns_R_matrix")} is \code{FALSE},
then a \code{big.matrix} is returned. The returned \code{big.matrix} is backed
by a temporary file mapping that will be deleted when the returned result is
garbage collected by R (see the examples).

Otherwise, a standard R matrix is returned. The dimensional shape of the output
is taken from \code{X}. If input \code{X} is dimensionless (that is, lacks a
dimension attribute), then the output is a column vector.
}
\references{
\url{http://www.netlib.org/blas/daxpy.f}
}
\author{
Michael J. Kane
}

\seealso{
\code{\link[bigmemory]{bigmemory}}
}
\examples{
require(bigmemory)
A = matrix(1, nrow=3, ncol=2)
B = big.matrix(nrow=3, ncol=2, type='double', init=1)
C = B + B   # C is a new big matrix
D = A + B   # D defaults to a regular R matrix, to change this, set the option:
            # optons(bigalgebra.mixed_arithmetic_returns_R_matrix=FALSE)
E = daxpy(A=1.0, X=B, Y=B)  # Same kind of result as C
print(C[])
print(D)
print(E[])

# The C and E big.matrix file backings will be deleted when garbage collected:
# (We enable debugging to see this explicitly)
options(bigalgebra.DEBUG=TRUE)
rm(C,E)
gc()
}
