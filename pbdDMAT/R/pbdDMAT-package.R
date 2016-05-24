#' Distributed Matrix Methods
#' 
#' A package for dense distributed matrix computations. Includes the use of
#' PBLAS and ScaLAPACK libraries via pbdSLAP, communicating over MPI via the
#' BLACS library and pbdMPI.
#' 
#' \tabular{ll}{ 
#'    Package: \tab pbdDMAT \cr 
#'    Type: \tab Package \cr 
#'    License: \tab GPL \cr 
#'    LazyLoad: \tab yes \cr 
#' } 
#' 
#' This package requires an MPI library (OpenMPI, MPICH2, or LAM/MPI).
#' 
#' @useDynLib pbdDMAT,
#'   R_int_sparse_count_zeros, R_sparse_count_zeros, 
#'   R_convert_dense_to_csr, R_convert_csr_to_dense
#' 
#' @import methods
#' @importFrom stats rnorm runif rweibull rexp
#' @importFrom utils head tail
#' @importFrom pbdMPI comm.cat comm.rank comm.print comm.size
#'    comm.stop comm.warning allreduce barrier comm.match.arg
#'    comm.any comm.all reduce allgather gather comm.max
#' @import pbdBASE
#' 
#' @name pbdDMAT-package
#' @docType package
#' @author Drew Schmidt \email{schmidt AT math.utk.edu}, Wei-Chen Chen, George
#' Ostrouchov, and Pragneshkumar Patel, with contributions from R Core team
#' (some wrappers taken from the base and stats packages).
#' @references Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' @keywords Package
NULL
