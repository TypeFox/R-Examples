#' Function that binds two FuzzyData together into one FuzzyData
#'
#' @description
#' This methods construct object \code{\linkS4class{FuzzyData}} based on two \code{\linkS4class{FuzzyData}}.
#' The functions merges the sources into single output. This method should be used in situations when both
#' weights and input data are fuzzy.
#'
#' @param data1 An object of \code{\linkS4class{FuzzyData}}.
#' @param data2 An object of \code{\linkS4class{FuzzyData}}.
#'
#' @return An object of class \code{\linkS4class{FuzzyData}}
#'
#' @export
#' @rdname bindColums-methods
#' @name bindColums
setGeneric("bindColums",
           function(data1, data2) standardGeneric("bindColums"))

#' @rdname bindColums-methods
#' @aliases bindColums,FuzzyData,FuzzyData-method
setMethod(
  f= "bindColums",
  signature(data1 = "FuzzyData", data2 = "FuzzyData"),
  definition=function(data1, data2)
  {
    fnMin = cbind(data1@fnMin, data2@fnMin)
    fnModal = cbind(data1@fnModal, data2@fnModal)
    fnMax = cbind(data1@fnMax, data2@fnMax)

    return(new("FuzzyData", fnMin = fnMin, fnModal = fnModal, fnMax = fnMax))
  }
)
