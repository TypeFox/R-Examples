#' Tools for manipulating linear systems of (in)equations
#'
#' @section Details:
#' 
#' This package offers a basic and consistent interface to a number of
#' operations on linear systems of (in)equations not available in base R. Except for the projection on
#' the convex polytope, operations are currently supported for dense matrices
#' only.
#'
#'
#' The following operations are implemented.
#' 
#' \itemize{
#'   \item{Split matrices in independent blocks}
#'   \item{Remove spurious rows and columns from a system of (in)equations}
#'   \item{Rewrite equalities in reduced row echelon form}
#'   \item{Eliminate variables through Gaussian or Fourier-Motzkin elimination}
#'   \item{Determine the feasibility of a system of linear (in)equations}
#'   \item{Compute Moore-Penrose Pseudoinverse}
#'   \item{Project a vector onto the convec polytope described by a set of linear (in)equations}
#'   \item{Simplify a system by substituting values}
#' }
#'
#' Most functions assume a system of (in)equations to be stored in a standard form. The \code{\link{normalize}}
#' function can bring any system of equations to this form.
#'
#'
#' @name lintools
#' @useDynLib lintools
#' @importFrom utils combn
#' @docType package
#' 
{}

