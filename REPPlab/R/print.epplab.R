#' Print an epplab Object
#' 
#' Prints an \code{epplab} object.
#' 
#' The print function displays the result with the best value in the objective
#' criterion.
#' 
#' @name print.epplab
#' @aliases print.epplab print-method print,epplab-method
#' @docType methods
#' @param x Object of class \code{epplab}.
#' @param ... Additional parameters
#' @author Daniel Fischer
#' @keywords methods print
#' @examples
#' 
#' library(tourr)
#' data(olive)
#' res <- EPPlab(olive[,3:10],PPalg="PSO",PPindex="KurtosisMin",n.simu=10, maxiter=20)
#' print(res)
#' 
#' @export
`print.epplab` <- function(x,...){
  X <- list()
  X$PPindex <- x$PPindex
  X$PPindexVal <- x$PPindexVal[1]
  X$PPalg <- x$PPalg
  X$PPdir <- x$PPdir[,1]
  X$PPiter <- x$PPiter[1]

  print(X,...)
} 
