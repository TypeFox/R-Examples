#' Class "Weights"
#'
#' @description An S4 class representing weights for AHP calculation. Each value in \code{numeric} vector
#' represents one weight.
#'
#' @slot weights Object of class \code{numeric} containing weights.
#'
#' @export
setClass(
  Class="Weights",

  representation(
    weights = "numeric"
  ),
  validity=function(object)
  {
     if(!(all.equal(sum(object@weights),1, tolerance=0.000001))){
       return(paste("Sum of weights is not equal to 1. The sum is ", sum(object@weights), ".", sep=""))
     }

    return(TRUE)
  }
)
