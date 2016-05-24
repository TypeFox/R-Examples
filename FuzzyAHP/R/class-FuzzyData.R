#' Class "FuzzyData"
#'
#' @description An S4 class to represent fuzzy data.
#'
#' @slot fnMin A numeric vector of minimal values of fuzzy data.
#' @slot fnModal A numeric vector of modal values of fuzzy data.
#' @slot fnMax A numeric vector of maximal values of fuzzy data.
#'
#' @export
setClass(
  Class="FuzzyData",

  representation(
    fnMin = "matrix",
    fnModal = "matrix",
    fnMax = "matrix"
  ),
  validity=function(object)
  {
    if(ncol(object@fnMin)!=ncol(object@fnMax) || ncol(object@fnModal)!=ncol(object@fnMax) || ncol(object@fnMin)!=ncol(object@fnModal)){
      return("The colums do not have the same number of columns.")
    }

    if(nrow(object@fnMin)!=nrow(object@fnMax) || nrow(object@fnModal)!=nrow(object@fnMax) || nrow(object@fnMin)!=nrow(object@fnModal)){
      return("The colums do not have the same number of rows.")
    }

    if(length(which(object@fnMin>object@fnModal)) != 0){
      return("Cannot create fuzzy data set. Minimal and modal values are not aligned correctly.")
    }

    if(length(which(object@fnModal>object@fnMax)) != 0){
      return("Cannot create fuzzy data set. Modal and maximal values are not aligned correctly.")
    }

    return(TRUE)
  }
)

#' Function that creates FuzzyData
#'
#' @description  This methods construct object \code{\linkS4class{FuzzyData}} based on provided \code{matrix}.
#' The matrix needs to be have rows represent individual fuzzy numbers and three colums that
#' represent minimal, modal and maximal value of fuzzy number.
#'
#' @param data A \code{matrix} with 3 colums.
#' @param single.value An optional boolean parameter (default value TRUE) specifying if the data to be
#' turn into fuzzy data is single vector of fuzzy numbers (then it needs to have 3 colums) or if the
#' whole matrix needs to be turn into fuzzy values.
#'
#' @return An object of class \code{\linkS4class{FuzzyData}}
#' @export
#'
#' @seealso \linkS4class{FuzzyData}
#'
#' @rdname fuzzyData-methods
#' @name fuzzyData
setGeneric("fuzzyData",
           function(data, single.value = TRUE) standardGeneric("fuzzyData"))

#' @rdname fuzzyData-methods
#' @aliases fuzzyData,matrix-method
setMethod(
  f="fuzzyData",
  signature(data = "matrix"),
  definition=function(data, single.value)
  {

    if(typeof(single.value)!="logical"){
      stop("Parameter single.value is not of type logical.")
    }

    colnames(data) = NULL
    rownames(data) = NULL

    if (single.value) {
      if(ncol(data)!=3){
        stop("Input dataset has to have three colums!")
      }

      return(new("FuzzyData", fnMin = data[,1], fnModal = data[,2], fnMax = data[,3]))
    }else{
      return(new("FuzzyData", fnMin = data, fnModal = data, fnMax = data))
    }

  }
)
