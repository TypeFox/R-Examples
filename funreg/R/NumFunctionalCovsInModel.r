#' @title Count the functional covariates in a model (for internal use by other package functions)
#' @description A very simple function, mainly for internal use by
#' package code, to count the number of functional 
#' covariates in an object of class \code{funreg}.  
#' @param object An object of class \code{funreg}, representing 
#' a fitted penalized functional regression with one or more
#' functional covariates. 
#' @return The number of functional covariates as an integer.
#'@export
num.functional.covs.in.model <- function(object) {
    stopifnot(class(object)=="funreg");
    p <- NA;
    if (is.matrix(object$betafn.estimate.by.grid)) {
        p <- ncol(object$betafn.estimate.by.grid);
    } else {
        p <- 1;
    }
    return(p);
}