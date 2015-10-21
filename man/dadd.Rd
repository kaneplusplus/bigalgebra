\name{dadd}
\Rdversion{1.1}
\alias{dadd}
\title{BLAS dadd functionality}
\description{
This function implements the function Y := ALPHA + Y where Y may be
either native double-precision valued R matrices or numeric vectors, or
double-precision valued \code{\link[bigmemory]{big.matrix}} objects, and ALPHA
is a scalar.
}
\usage{
dadd(Y, ALPHA, SIGN=1, ALPHA_LHS=1)
}
\arguments{
  \item{Y}{Requried to be either a native \R \code{\link{matrix}}
           or numeric vector,
           or a \code{\link[bigmemory]{big.matrix}} object}
  \item{ALPHA}{Required to be a scaler.}
  \item{SIGN}{1 for addition -1 for subtraction.}
  \item{ALPHA_LHS}{Should SIGN be multiplied by Y (default) or alpha.}
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
B = dadd(A, 10)
print(A)
print(B)
}
