#' The numericMpfr class
#'
#' The S4 class union of numeric and mpfr, primarily used to define slots in ecd class.
#' The use of MPFR does not necessarily increase precision. Its major strength in ecd 
#' is ability to handle very large numbers when studying asymptotic behavior, and 
#' very small numbers caused by small sigma when studying high frequency option data. 
#' Since there are many convergence issues with integrating PDF using native integrateR library,
#' the ecd package adds many algorithms to improve its performance. These additions
#' may decrease precision (knowningly or unknowningly) for the sake of increasing performance. 
#' More research is certainly needed in order to cover a vast range of parameter space!
#' 
#' @name numericMpfr-class
#'
#' @importClassesFrom Rmpfr mpfr mpfrArray mpfrMatrix
#'
#' @exportClass numericMpfr
setClassUnion("numericMpfr", c("numeric", "mpfr", "mpfrArray"))

# end

