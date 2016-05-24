# References
# The Eigenvalues of Random Matrices Experiments with the Classical Ensembles
#   http://web.mit.edu/18.338/www/handouts/handout3.pdf
# Introduction to the Random Matrix Theory: Gaussian Unitary Ensemble and Beyond
#   http://arxiv.org/abs/math-ph/0412017v2
# Tyler's M-Estimator, Random Matrix Theory, and Generalized Elliptical 
# Distributions with Applications to Finance
#   http://papers.ssrn.com/sol3/papers.cfm?abstract_id=1287683
#   http://web.mit.edu/sea06/agenda/talks/Kuijlaars.pdf
# Distributions of the extreme eigenvalues of the complex Jacobi random matrix
# ensemble
#   http://www.math.washington.edu/~dumitriu/kd_submitted.pdf
#   http://www.williams.edu/go/math/sjmiller/public_html/BrownClasses/54/handouts/IntroRMT_Math54.pdf
#   http://web.mit.edu/sea06/agenda/talks/Harding.pdf
#   http://dspace.mit.edu/bitstream/handle/1721.1/39670/180190294.pdf
#   http://projecteuclid.org/DPubS?service=UI&version=1.0&verb=Display&handle=euclid.cmp/1103842703
#   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.10.2630
#   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.5.9005
#
# Resources
#   http://www.math.ucsc.edu/research/rmtg.html
#   The Distribution Functions of Random Matrix Theory, Craig A. Tracy, UC Davis
#
# Search Terms
#   extreme eigenvalues random matrix feature extraction
#
# Trading Ideas
#   http://ideas.repec.org/p/lei/ingber/03ai.html
#   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.7.1041

#' Perform the conjugate transpose of a matrix
#'
#' @section Usage:
#' ct(m) \%::\% matrix : matrix
#'
#' ct(m)
#' 
#' @section Details:
#' This is a convenience function to compute the conjugate transpose. For
#' real-valued matrices, ct(m) == t(m).
#'
#' @name ct
#' @param m A matrix
#' @return THe conjugate transpose of the original matrix
#' @examples
#' x <- rcomp(10)
#' ct(x)
ct(m) %::% matrix : matrix
ct(m) %as% { Conj(t(m)) }


#' Generation of random complex numbers
#'
#' Generate random complex numbers using the specified distribution.
#' By default \code{rnorm} is used.
#'
#' @section Usage:
#' rcomp(n, dist) %::% numeric : Function : complex
#' 
#' rcomp(n, dist=rnorm)
#'
#' @section Details:
#' This function is used primarily to generate random matrices.
#'
#' @name rcomp
#' @param n Length of the output vector
#' @param dist The distribution for the random number genertor
#' @return A vector of random numbers
#' @examples
#' rcomp(10)
#' rcomp(10, runif)
rcomp(n, dist) %::% numeric : Function : complex
rcomp(n, dist=rnorm) %as% {
  complex(real=dist(n), imaginary=dist(n))
}

#' Generation of random matrices
#'
#' Generate various types of random matrices
#'
#' @section Usage:
#' rmatrix(model) %::% WignerModel : matrix
#'
#' rmatrix(model) %::% WishartModel : matrix
#'
#' rmatrix(model) %::% JacobiModel : matrix
#'
#' rmatrix(model)
#'
#' @section Details:
#' Gaussian Orthogonal Ensemble
#' Gaussian Unitary Ensemble
#'
#' @name rmatrix
#' @param model The matrix model to use, which includes the size of
#' the matrix. The model argument must be of type RandomMatrixModel. 
#' Numerous sub-types (e.g. WignerModel, WishartModel) are
#' supported generating the appropriate type of random matrix.
#'
# Also known as a Gaussian Orthogonal Ensemble
rmatrix(model) %::% WignerModel : matrix
rmatrix(model) %when% {
  !model$real
  # result == ct(result)
} %as% {
  n <- model$n
  x <- matrix(rcomp(n^2), nrow=n) / 2^0.5
  (x + ct(x)) / sqrt(2 * n)
}

# Also known as a Gaussian Unitary Ensemble
rmatrix(model) %::% WignerModel : matrix
rmatrix(model) %when% {
  model$real
  # result == t(result)
} %as% {
  n <- model$n
  x <- matrix(rnorm(n^2), nrow=n)
  (x + ct(x)) / sqrt(2 * n)
}

# For definition of self-dual quaternion and matrix representation,
# http://www.aimath.org/conferences/ntrmt/talks/Mezzadri2.pdf
# http://tonic.physics.sunysb.edu/~verbaarschot/lecture/lecture2.ps

rmatrix(model) %::% WishartModel : matrix
rmatrix(model) %when% {
  !model$real
} %as% {
  n <- model$n
  m <- model$m
  dist.fn <- function(x) rnorm(x, sd=model$sd)
  x <- matrix(rcomp(n * m, dist=dist.fn), nrow=n) / 2^0.5
  (x %*% ct(x)) / m
}

rmatrix(model) %::% WishartModel : matrix
rmatrix(model) %when% {
  model$real
} %as% {
  n <- model$n
  m <- model$m
  x <- matrix(rnorm(n * m, sd=model$sd), nrow=n)
  (x %*% t(x)) / m
}


rmatrix(model) %::% JacobiModel : matrix
rmatrix(model) %when% {
  model$real
} %as% {
  n <- model$n
  m1 <- model$m1
  m2 <- model$m2

  x1 <- rmatrix(WishartModel(n,m1, real=model$real))
  x2 <- rmatrix(WishartModel(n,m2, real=model$real))
  solve(x1 + x2) * x1
}


RandomMatrixModel(real=TRUE, ...) %as% list(real=real, ...)

# Random square matrix. Eienvalues form semicircle
WignerModel(n, ...) %as%
{
  RandomMatrixModel(n=n, ...)
}

# n - variables
# m - observations
# model <- WishartModel(100,500, sd=1)
# hist(eigenvalues(rmatrix(model)))
WishartModel(n, m, sd=1, ...) %as%
{
  RandomMatrixModel(n=n, m=m, Q=m/n, sd=sd, ...)
}

JacobiModel(n, m1, m2, ...) %as%
{
  RandomMatrixModel(n=n, m1=m1, m2=m2, ...)
}

Ensemble(count, model) %as%
{
  out <- lapply(seq(count), function(junk) rmatrix(model))
  out@model <- class(model)[1]
  out
}

print.Ensemble <- function(x, ...)
{
  cat("\nClass:", attr(x,'model'))
  cat("\nCount:", length(x))
  cat("\nDimensions:", dim(x[[1]]))
  cat("\n")
  invisible(x)
}



eigenvalues(m) %::% matrix : numeric
eigenvalues(m) %as% {
  o <- eigen(m, only.values=TRUE)
  o$values
}

# Example
# en <- Ensemble(50, WignerModel(200))
# hist(max_eigen(en), freq=FALSE)
max_eigen(ensemble) %::% Ensemble : numeric
max_eigen(ensemble) %as% {
  es <- lapply(ensemble, eigen, only.values=TRUE)
  sapply(es, function(x) x$values[1])
}



# Convenience functions
hermitian(n) %as% rmatrix(WignerModel(n=n, real=FALSE))

symmetric(n) %as% rmatrix(WignerModel(n=n, real=TRUE))

