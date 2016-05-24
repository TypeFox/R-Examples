#' Class "FuzzyWeights"
#'
#' @description An S4 class to represent fuzzy weights for fuzzy AHP calculation.
#'
#' @slot fnMin Object of class \code{numeric} containing minimal values of fuzzy weights.
#' @slot fnModal Object of class \code{numeric} containing modal values of fuzzy weights.
#' @slot fnMax Object of class \code{numeric} containing maximal values of fuzzy weights.
#'
#' @export
setClass(
  Class="FuzzyWeights",

  representation(
    fnMin = "numeric",
    fnModal = "numeric",
    fnMax = "numeric"
  ),
  validity=function(object)
  {
    if(length(object@fnMin) != length(object@fnModal) || length(object@fnModal) != length(object@fnMax)){
      return("Lenghts of fuzzy values representing minimal, modal and maximal values must be the same!")
    }


    if(length(which(object@fnMin>object@fnModal)) != 0){
      return("Cannot create fuzzy data set. Minimal and modal values are not aligned correctly.")
    }

    if(length(which(object@fnModal>object@fnMax)) != 0){
      return("Cannot create fuzzy data set. Modal and maximal values are not aligned correctly.")
    }

    if(!(all.equal(sum(object@fnModal),1, tolerance=0.000001))){
      return("Sum of modal values of weights is not equal to 1.")
    }

    return(TRUE)
  }
)
