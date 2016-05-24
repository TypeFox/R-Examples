#' Package to Facilitate ODE-Based Modeling
#'
#' This package provides methods to
#' \itemize{
#'   \item{} import a conceptual ODE-based model stored in tabular form (i.e.
#'     as text files or spreadsheets).
#'   \item{} generate source code (either R or Fortran) to be passed to an
#'     ODE-solver.
#'   \item{} visualize and export basic information about a model, e.g. for
#'     documentation purposes.
#' }
#'
#' Consult the package vignette for details. The concept of writing an ODE
#'   system in tabular/matrix form is nicely introduced, e. g., in the book of
#'   Reichert, P., Borchardt, D., Henze, M., Rauch, W., Shanahan, P.,
#'   Somlyody, L., and Vanrolleghem, P. A. (2001): River water quality model
#'   No. 1, IWA publishing, ISBN 9781900222822.
#'
#' The current source code repository is \url{https://github.com/dkneis/rodeo}.
#'
#' @name rodeo-package
#' @docType package
#'
#' @section Class and class methods:
#' See \code{\link{rodeo-class}} for the \code{rodeo} reference class
#'   and the corresponding class methods.
#'
#' @section Non-class methods:
#'   Type \code{help(package="rodeo")} or see the links below to access the
#'   documentation of non-class methods contained in the package.
#'
#' \itemize{
#'   \item{\code{\link{solverInterface}}} Generation of Fortran wrapper code for
#'     use with the numerical solvers from packages
#'     \code{\link[deSolve]{deSolve}} and \code{\link[rootSolve]{rootSolve}}.
#'   \item{\code{\link{forcingFunctions}}} Generation of forcing functions
#'     in Fortran.
#'   \item{\code{\link{exportDF}}} Export of data frames as TEX or HTML code.
#' }
#' 
#' @author \email{david.kneis@@tu-dresden.de}
NULL
