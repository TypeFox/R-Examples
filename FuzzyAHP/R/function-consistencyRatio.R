#' Function to determine Consistency Ratio
#'
#' @description
#' This methods calculates Consistency Ratio for \code{\linkS4class{PairwiseComparisonMatrix}}.
#' The consistency ratio can only be provided for \code{\linkS4class{PairwiseComparisonMatrix}} with
#' less than 10 rows. For bigger matrices the value is not known.
#'
#' @param comparisonMatrix A \code{\linkS4class{PairwiseComparisonMatrix}}
#' @param print.report Optional boolean parameter stating if short report should be printed along with determination
#' of Consistency Ratio. Default value is \code{TRUE}.
#'
#' @return A numeric value of Consistency Ratio, for \code{\linkS4class{PairwiseComparisonMatrix}} with more than 10
#' an error is raised.
#'
#' @details Generally pairwise comparison matrixes are considered to be consistent if the value of Consistency Ration
#' is smaller than 0.1. For matrices comparing more then 10 elements then Consistency Ratio is unsuitable, because the
#' values of random index, that is necessary to obtain Consistency Ratio, are only known for matrixes with size smaller
#' than \eqn{10\times10}.
#'
#' @export
#' @rdname consistencyRatio-methods
#' @name consistencyRatio
setGeneric("consistencyRatio",
           function(comparisonMatrix, print.report = TRUE) standardGeneric("consistencyRatio"))

#' @rdname consistencyRatio-methods
#' @aliases consistencyRatio,PairwiseComparisonMatrix-method
setMethod(
  f="consistencyRatio",
  signature(comparisonMatrix = "PairwiseComparisonMatrix"),
  definition=function(comparisonMatrix, print.report)
  {
    # random index for 10 values
    # randomIndex = c(0, 0, 0.52, 0.89, 1.11, 1.25, 1.35, 1.4, 1.45, 1.49)

    # random index for 15 values
    randomIndex = c(0, 0, 0.52, 0.89, 1.11, 1.25, 1.35, 1.4, 1.45, 1.49, 1.52, 1.54, 1.56, 1.58, 1.59)

    CI = consistencyIndex(comparisonMatrix)

    if(nrow(comparisonMatrix@values)<= length(randomIndex)){
      CR = CI / (randomIndex[nrow(comparisonMatrix@values)])

      if(print.report){
        if(CR<0.1){
          cat(paste("Consistency ratio is: ", CR, ". The pairwise comparison matrix is consistent for calculations.", sep=""))
          cat("\n")
        }else{
          cat(paste("Consistency ratio is: ", CR, ". It should be lower than: ", randomIndex[nrow(comparisonMatrix@values)],
                        ". The pairwise comparison matrix is not consistent enough for correct calculations. Please consider redefining the matrix!",
                        sep = ""))
          cat("\n")
        }
      }
      return(CR)

    }else{
      stop(paste("Cannot calculate Consistency Ratio for matrices with more then ", length(randomIndex), " rows! The Consistency index is ",
                    CI, ". The consistency of this matrix cannot be determined by Consistency Ratio.",
                 " Please consider using another method to determined if the matrix is consistent.", sep = ""))
    }
  }
)

#' @rdname consistencyRatio-methods
#' @aliases consistencyRatio,FuzzyPairwiseComparisonMatrix-method
setMethod(
  f="consistencyRatio",
  signature(comparisonMatrix = "FuzzyPairwiseComparisonMatrix"),
  definition=function(comparisonMatrix, print.report)
  {
    randomIndex = c(0, 0, 0.52, 0.89, 1.11, 1.25, 1.35, 1.4, 1.45, 1.49)

    CI = consistencyIndex(comparisonMatrix)

    if(nrow(comparisonMatrix@fnModal)<11){
      CR = CI / (randomIndex[nrow(comparisonMatrix@fnModal)])

      if(print.report){
        if(CR<0.1){
          cat(paste("Consistency ratio of modal values is: ", CR, ". The fuzzy pairwise comparison matrix is consistent for calculations.", sep=""))
          cat("\n")
        }else{
          cat(paste("Consistency ratio of modal values is: ", CR,
                    "The fuzzy pairwise comparison matrix is not consistent enough for correct calculations. Please consider redefining the matrix!",
                    sep = ""))
          cat("\n")
        }
      }
      return(CR)

    }else{
      stop(paste("Cannot calculate Consistency Ratio for matrices with more then 10 rows! The Consistency index of modal values is ",
                 CI, ". The consistency of this fuzzy matrix cannot be determined by Consistency Ratio.",
                 " Please consider using another method to determined if the fuzzy matrix is consistent.", sep = ""))
    }
  }
)
