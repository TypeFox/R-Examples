#' Demonstrations and Examples for the pbd Project
#' 
#' A set of demos of pbdR packages, together with a useful,
#' unifying vignette.
#' 
#' \tabular{ll}{ Package: \tab pbdDEMO\cr Type: \tab Package\cr License: \tab
#' MPL 2.0\cr LazyLoad: \tab yes\cr } This package requires an MPI library
#' (OpenMPI, MPICH2, or LAM/MPI).
#' 
#' @import methods 
#' @importFrom pbdMPI comm.cat comm.print comm.stop allreduce 
#'   allgather comm.size comm.rank comm.stop spmd.allreduce.integer
#'   comm.any comm.all comm.warning spmd.allgather.integer
#' @import pbdBASE pbdDMAT
#' @importFrom graphics image
#' @importFrom stats uniroot
#' @importFrom utils read.csv
#' 
#' @useDynLib pbdDEMO, pbddemo_linecount
#' 
#' @name pbdDEMO-package
#' @docType package
#' @author Drew Schmidt \email{schmidt AT math.utk.edu}, Wei-Chen Chen, George
#' Ostrouchov, and Pragneshkumar Patel.
#' @references Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' @keywords Package
NULL


