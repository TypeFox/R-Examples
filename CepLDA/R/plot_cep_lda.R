#' @name plot.ceplda
#' @title Plot Method for Class 'ceplda'
#' 
#' @description
#' Plots a set of data on one, two or more linear discriminants.
#' @param x an object of class "ceplda".
#' @param ... additional arguments.
#' @method plot ceplda
#' @details This function is a method for the generic function \code{plot()} for class "\code{ceplda}". It can be invoked by calling \code{plot(x)}
#' for an object \code{x} of the appropriate class, or directly by calling \code{plot.ceplda(x)} regardless of the class of the object. Details usage of this function
#' is equivalent to \code{plot.lda(MASS)}.
#' @export
#' @importFrom graphics plot
#' @return NULL
#' @seealso
#'  \code{\link{cep.lda}}, \code{\link{plot.lda}} 

plot.ceplda <- function(x,...){
    if(!inherits(x,"ceplda"))
      stop("Object must be of class 'ceplda'")
    if(length(x$C.lda)!=10)
      stop("Object muse be from lda without cross-validation")
    plot(x$C.lda,...)
}