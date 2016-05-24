#' Print CEC
#' @export
#' @rdname print.cec
#' @method print Rcpp_CecModel 
#'
#' @title print
#' 
#' @description Print basic information about clusters found.
#' Presents a structure of the cec results object (clusters found)
#'
#' @docType methods
#'
#' @param x CEC object model.
#' @param ... other arguments not used by this method.
#' 
print.Rcpp_CecModel <- function(x, ...) {
    print(sprintf("CEC clustering; %d clusters with energy = %f",
                  length(x$centers), x$energy))
    print("Centers: ")
    print(x$centers)
    print("Covariances: ")
    print(x$covMatrix)
}

show.Rcpp_CecModel <- function(object){
  print(object)
}

setMethod("show", "Rcpp_CecModel", show.Rcpp_CecModel)
