#' Function to extract specific fuzzy numbers
#'
#' @description
#' This methods helps with extracting fuzzy numbers from \code{\linkS4class{FuzzyData}} and
#' \code{\linkS4class{FuzzyWeights}}.
#'
#' @param object An object of class \code{\linkS4class{FuzzyData}} or \code{\linkS4class{FuzzyWeights}}
#' @param index An object of class \code{integer} that represents one or more indices to extract the data from
#'
#' @return A \code{matrix} where rows
#'
#' @export
#' @rdname getFuzzyNumber-methods
#' @name getFuzzyNumber
setGeneric("getFuzzyNumber",
           function(object,index) standardGeneric("getFuzzyNumber"))

#' @rdname getFuzzyNumber-methods
#' @aliases getFuzzyNumber,FuzzyData,integer-method
setMethod(
  f="getFuzzyNumber",
  signature(object = "FuzzyData", index = "integer"),
  definition=function(object, index)
  {


    if(max(index)>nrow(object@fnMin)){
      stop("Maximal value of index is outside of range of number of fuzzy data!")
    }

    data = cbind(object@fnMin[index,],object@fnModal[index,],object@fnMax[index,])

    colnames(data) = c("minimal","modal","maximal")

    return(data)
  }
)

#' @rdname getFuzzyNumber-methods
#' @aliases getFuzzyNumber,FuzzyWeights,integer-method
setMethod(
  f="getFuzzyNumber",
  signature(object = "FuzzyWeights", index = "integer"),
  definition=function(object, index)
  {

    if(max(index)>nrow(object@fnMin)){
      stop("Maximal value of index is outside of range of number of fuzzy weights!")
    }

    data = cbind(object@fnMin[index,],object@fnModal[index,],object@fnMax[index,])

    colnames(data) = c("minimal","modal","maximal")

    return(data)
  }
)
