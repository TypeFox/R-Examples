#' Function to calculate result of AHP
#'
#' @description  This function calculates output of AHP based on \code{\linkS4class{Weights}}
#' or \code{\linkS4class{FuzzyWeights}} on data represented either by \code{matrix} or
#' \code{\linkS4class{FuzzyData}}.
#'
#' @param weights object of class \linkS4class{Weights} or \linkS4class{FuzzyWeights}. Alternatively objects of
#' classes \linkS4class{PairwiseComparisonMatrix} or \linkS4class{FuzzyPairwiseComparisonMatrix} can be passed to
#' directly calculate weights from these classes.
#'
#' @param data matrix or \linkS4class{FuzzyData} with number of colums equal to number of rows in \code{weights}.
#'
#' @return Either a matrix (if \linkS4class{Weights} and \code{matrix} were used as inputs) or
#' \linkS4class{FuzzyData} (if \linkS4class{FuzzyWeights} were used).
#'
#' @export
#' @rdname calculateAHP-methods
#' @name calculateAHP
setGeneric("calculateAHP",
           function(weights, data) standardGeneric("calculateAHP"))

#' @rdname calculateAHP-methods
#' @aliases calculateAHP,Weights,matrix-method
setMethod(
  f="calculateAHP",
  signature(weights = "Weights", data = "matrix"),
  definition=function(weights, data)
  {
    if(!(ncol(data)==length(weights@weights))){
      stop("Can not multiply the data by the weights, the number of data columns does not match the number of weights.")
    }

    numberRows = nrow(data)

    result = matrix(nrow = numberRows, ncol = 1)

    result[, 1] = rowSums(t(t(data[,]) * weights@weights))

    colnames(result) = c("result")

    return(result)
  }
)

#' @rdname calculateAHP-methods
#' @aliases calculateAHP,FuzzyWeights,matrix-method
setMethod(
  f="calculateAHP",
  signature(weights = "FuzzyWeights", data = "matrix"),
  definition=function(weights, data)
  {
    if(!(ncol(data)==length(weights@fnMin))){
      stop("Can not multiply the data by the weights, the number of data columns does not match the number of weights.")
    }

    numberRows = nrow(data)

    result = matrix(nrow = numberRows, ncol = 3)

    result[, 1] = rowSums(t(t(data[,]) * weights@fnMin))
    result[, 2] = rowSums(t(t(data[,]) * weights@fnModal))
    result[, 3] = rowSums(t(t(data[,]) * weights@fnMax))

    return(new("FuzzyData", fnMin = as.matrix(result[,1]), fnModal = as.matrix(result[,2]), fnMax = as.matrix(result[,3])))
  }
)

#' @rdname calculateAHP-methods
#' @aliases calculateAHP,FuzzyWeights,FuzzyData-method
setMethod(
  f="calculateAHP",
  signature(weights = "FuzzyWeights", data = "FuzzyData"),
  definition=function(weights, data)
  {
    if(!(ncol(data@fnModal)==length(weights@fnModal))){
      stop("Can not multiply the fuzzy data by the weights, the number of data columns does not match the number of weights.")
    }

    numberRows = nrow(data@fnModal)

    result = matrix(nrow = numberRows, ncol = 3)

    result[, 1] = rowSums(t(t(data@fnMin[,]) * weights@fnMin))
    result[, 2] = rowSums(t(t(data@fnModal[,]) * weights@fnModal))
    result[, 3] = rowSums(t(t(data@fnMax[,]) * weights@fnMax))

    return(new("FuzzyData", fnMin = as.matrix(result[,1]), fnModal = as.matrix(result[,2]), fnMax = as.matrix(result[,3])))
  }
)

#' @rdname calculateAHP-methods
#' @aliases calculateAHP,PairwiseComparisonMatrix,matrix-method
setMethod(
  f="calculateAHP",
  signature(weights = "PairwiseComparisonMatrix", data = "matrix"),
  definition=function(weights, data)
  {
    return (calculateAHP(calculateWeights(weights), data))
  }
)

#' @rdname calculateAHP-methods
#' @aliases calculateAHP,FuzzyPairwiseComparisonMatrix,matrix-method
setMethod(
  f="calculateAHP",
  signature(weights = "FuzzyPairwiseComparisonMatrix", data = "matrix"),
  definition=function(weights, data)
  {
    return (calculateAHP(calculateWeights(weights), data))
  }
)

#' @rdname calculateAHP-methods
#' @aliases calculateAHP,FuzzyPairwiseComparisonMatrix,FuzzyData-method
setMethod(
  f="calculateAHP",
  signature(weights = "FuzzyPairwiseComparisonMatrix", data = "FuzzyData"),
  definition=function(weights, data)
  {
    return (calculateAHP(calculateWeights(weights), data))
  }
)
