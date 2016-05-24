#' @title coef method for funreg object
#' @description Returns coefficient information on a \code{funreg} object.
#' @param object An object of class \code{funreg}
#' @param digits The number of digits past the decimal place to use when printing numbers
#' @param silent If \code{TRUE}, indicates that the summary should be returned
#' as a list object but not printed to the screen.
#' @param ... Other arguments that may be passed from another method.
#' @return At least for now, this is identical to the \code{summary.funreg}
#' function. 
#'@export 
#'@S3method coef funreg
#'@method coef funreg
coef.funreg <- function(object, 
                           digits=4,
                           silent=FALSE, ...) {
    stopifnot(class(object)=="funreg");
    return(summary.funreg(object,digits=digits,silent=silent));
}