#' ScaLAPACK Wrappers and Utilities
#' 
#' A package contains the basic methods for dealing with distributed data
#' types, as well as the data types themselves.
#' 
#' \tabular{ll}{ Package: \tab pbdBASE\cr Type: \tab Package\cr License: \tab
#' MPL\cr LazyLoad: \tab yes\cr } This package requires an MPI library
#' (OpenMPI, MPICH2, or LAM/MPI).
#' 
#' @import methods pbdSLAP
#' @importFrom pbdMPI allreduce comm.print comm.stop comm.rank comm.warning comm.is.null bcast
#' 
#' @useDynLib pbdBASE,
#'   R_nbd,
#'   R_NUMROC, R_PDLAPRNT, g2l_coords, l2g_coords,
#'   R_optimal_grid, R_blacs_init, R_igsum2d1, R_dgsum2d1, R_igamx2d1, 
#'   R_dgamx2d1, R_igamn2d1, R_dgamn2d1, R_dgesd2d1, R_dgerv2d1,
#'   R_PDGELS,
#'   R_matexp, R_p_matpow_by_squaring, R_p_matexp_pade,
#'   R_PDTRAN, R_PDGEMM, R_PDCROSSPROD, R_PDCHTRI, R_PDCLVAR,
#'   R_MKSUBMAT, R_MKGBLMAT, R_DALLREDUCE, R_PTRI2ZERO, R_PDSWEEP, R_RL2BLAS, 
#'   R_RL2INSERT, R_PDGDGTK, R_PDDIAGMK, R_RCOLCPY, R_RCOLCPY2, R_RROWCPY, 
#'   R_RROWCPY2, R_PDMVSUM, R_DHILBMK, R_PDHILBMK, R_PDMKCPN1, 
#'   R_PDGEQPF, R_PDORGQR, R_PDORMQR, 
#'   R_PDGETRI, R_PDGESV, R_PDGESVD, R_PDSYEV, R_PDPOTRF, R_PDSYEVX, 
#'   R_PDGETRF, R_PDLANGE, R_PDTRCON, R_PDGECON, R_PDGEMR2D,
#'   R_g2lcoord
#' 
#' @name pbdBASE-package
#' @docType package
#' @author Drew Schmidt \email{schmidt AT math.utk.edu}, Wei-Chen Chen, George
#' Ostrouchov, and Pragneshkumar Patel.
#' @references Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' @keywords Package
NULL

