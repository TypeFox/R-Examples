#' @title Print method for an object of class \code{orthoProjection}
#' @description Prints the contents of an object of class \code{orthoProjection}
#' @aliases print.orthoProjection
#' @usage \method{print}{orthoProjection}(x, ...)
#' @param x an object of class \code{orthoProjection} (as returned by the \code{orthoProjection} function). 
#' @param ... arguments to be passed to methods (not yet functional).
#' @author Leonardo Ramirez-Lopez
#' @export

print.orthoProjection <- function(x, ...){

  cat("\n", "Method: ",x$method)
  cat("\n", "Number of components retained: ",x$n.components,"\n")
  cat(" Number of observations and number of original variables: ", c(nrow(x$scores), ncol(x$X.loadings)), "\n")
  
  cat( "\n", "Standard deviations, cumulative variance explained, individual variance explained:", "\n")
  if(x$method == "pls")
    print(x$variance$x.var, digits = 3)
  else
    print(x$variance, digits = 3)
}
