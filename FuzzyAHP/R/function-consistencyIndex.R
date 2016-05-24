#' Function to determine Consistency Index
#'
#' @description
#' This methods calculates Consistency index for \code{\linkS4class{PairwiseComparisonMatrix}}.
#'
#' @param comparisonMatrix A \code{\linkS4class{PairwiseComparisonMatrix}}
#'
#' @return A numeric value of Consistency index.
#'
#' @export
#' @rdname consistencyIndex-methods
#' @name consistencyIndex
setGeneric("consistencyIndex",
           function(comparisonMatrix) standardGeneric("consistencyIndex"))

#' @rdname consistencyIndex-methods
#' @aliases consistencyIndex,PairwiseComparisonMatrix-method
setMethod(
  f="consistencyIndex",
  signature(comparisonMatrix = "PairwiseComparisonMatrix"),
  definition=function(comparisonMatrix)
  {
    # CI = (Re(eigen(comparisonMatrix@values)$values[1]) - nrow(comparisonMatrix@values)) / (nrow(comparisonMatrix@values) - 1)
    CI = .CI(comparisonMatrix@values)

    return(CI)
  }
)

#' @rdname consistencyIndex-methods
#' @aliases consistencyIndex,FuzzyPairwiseComparisonMatrix-method
setMethod(
  f="consistencyIndex",
  signature(comparisonMatrix = "FuzzyPairwiseComparisonMatrix"),
  definition=function(comparisonMatrix)
  {
    CI = .CI(comparisonMatrix@fnModal)

    return(CI)
  }
)





setGeneric(".CI",
           function(matrix) standardGeneric(".CI"))

setMethod(
  ".CI",
  signature(matrix = "matrix"),
  function (matrix)
  {
    CI = (Re(eigen(matrix)$values[1]) - nrow(matrix)) / (nrow(matrix) - 1)

    return(CI)
  }
)
