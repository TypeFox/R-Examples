#' Function to calculate fuzzy weights based on comparison matrix
#'
#' @description This functions calculates \code{\linkS4class{Weights}} or \code{\linkS4class{FuzzyWeights}}
#' based on input pairwise comparison matrix.
#'
#' @param comparisonMatrix object of either \linkS4class{PairwiseComparisonMatrix} or \linkS4class{FuzzyPairwiseComparisonMatrix}
#'
#' @seealso \link{PairwiseComparisonMatrix-class}
#'
#' @export
#' @rdname calculateWeights-methods
#' @name calculateWeights
setGeneric("calculateWeights",
           function(comparisonMatrix) standardGeneric("calculateWeights"))

#' @rdname calculateWeights-methods
#' @aliases calculateWeights,PairwiseComparisonMatrix-method
setMethod(
  f="calculateWeights",
  signature(comparisonMatrix="PairwiseComparisonMatrix"),
  definition=function(comparisonMatrix)
  {
    weightsCount = nrow(comparisonMatrix@values)

    weights = vector(mode = "numeric", length = weightsCount)

    for(i in 1:weightsCount){
      weights[i] = prod(comparisonMatrix@values[i,])^(1/weightsCount)
    }

    sumWeights = sum(weights)

    weights = weights/sumWeights

    rNames = c()
    for (i in 1:length(weights)){
      rNames = append(rNames, paste("w",i, sep = ""))
    }

    names(weights) = rNames

    return (new("Weights", weights = weights))
  }
)

#' @rdname calculateWeights-methods
#' @aliases calculateWeights,FuzzyPairwiseComparisonMatrix-method
setMethod(
  f="calculateWeights",
  signature(comparisonMatrix="FuzzyPairwiseComparisonMatrix"),
  definition=function(comparisonMatrix)
  {
    p = nrow(comparisonMatrix@fnMin)

    wMin = c()
    wModal = c()
    wMax = c()

    for(i in 1:p){
      limits = .weightsLimits(comparisonMatrix,i)
      wMin = append(wMin, limits[1])
      wMax = append(wMax, limits[2])

      sum = 0

      for(k in 1:p){
        sum = sum + prod(comparisonMatrix@fnModal[k,])^(1/p)
      }

      wM = ((prod(comparisonMatrix@fnModal[i,])^(1/p)) / sum)

      wModal = append(wModal, wM)
    }

    return (new("FuzzyWeights", fnMin = wMin, fnModal = wModal, fnMax = wMax))
  }
)
