#' @title Wrapper for the SciMark 2.0 benchmark.
#'
#' @description This function is a simple wrapper around the ANSI C
#' version of the \href{http://math.nist.gov/scimark}{SciMark 2.0 benchmark}
#' which is a benchmark for numerical and scientific computing. Concicely
#' performance measurements for the computational kernels \emph{Fast Fourier
#' Transformation} (FFT), \emph{Gauss-Seidel relaxation}, \emph{Sparse matrix-multiply},
#' \emph{Monte Carlo integration} and \emph{dense LU factorization} are computed.
#'
#' In order to isolate effects of memory hierarchy the problem sizes, e.g.,
#' the size of the matrix fpr the dense LU matrix factorization, are pretty small.
#' However, addressing the performance of the memory subsystem is possible by
#' setting the \code{large} argument to \code{TRUE}.
#'
#' @param large [\code{logical(1)}]\cr
#'  Run large version of benchmark?
#'  Default is \code{FALSE}.
#' @param minimum.time [\code{numeric(1)}]\cr
#'   Minimum time to run each of the benchmarks, in seconds.
#'   Default is 2.
#' @return [\code{numeric}] Named vector of time measurements with the
#' following components:
#' \describe{
#'   \item{Composite}{Mean value of the remaining components.}
#'   \item{FFT}{Performance of the Fast Fourier Transformation (FFT).}
#'   \item{SOR}{Performance of the Jacobi Successive Over-relaxation (SOR).}
#'   \item{MC}{Performance of a Monte Carlo integration.}
#'   \item{SMM}{Performance of a spare matrix multiplication.}
#'   \item{LU}{Performance of a dense LU matrix factorization.}
#' }
#' @export
#' @useDynLib rscimark c_rscimark
rscimark = function(large = FALSE, minimum.time = 2) {
  assertFlag(large)
  assertNumber(minimum.time, lower = 0)
  y = .Call(c_rscimark, large, minimum.time)
  names(y) = c("Composite", "FFT", "SOR", "MC", "SMM", "LU")
  return(y)
}
