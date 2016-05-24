#' Function to assess strict consistency of Comparison Matrix
#'
#' @description Check if \deqn{a_{ik} = a_{ij} \times a_{jk}}{a[i,k]==a[i,j]*a[j,k]} applies
#' for all \eqn{i,j,k = 1,2,\dots,n}, where \eqn{n} is size of \eqn{a}.
#'
#' @param PairwiseComparisonMatrix A \linkS4class{FuzzyPairwiseComparisonMatrix} or
#' \linkS4class{PairwiseComparisonMatrix}.
#' @param print.report Optional boolean parameter stating if short report should be printed along with determination
#' of Weak Consistency. Default value is \code{TRUE}.
#'
#' @return Boolean value indicating if Comparison Matrix passed the weak consistency test and a warning message
#' listing the problematic triplets if the matrix is not consisten.
#'
#' @export
#' @rdname strictConsistency-methods
#' @name strictConsistency
setGeneric("strictConsistency",
           function(PairwiseComparisonMatrix, print.report = TRUE) standardGeneric("strictConsistency"))

#' @rdname strictConsistency-methods
#' @aliases strictConsistency,FuzzyPairwiseComparisonMatrix-method
setMethod(
  f="strictConsistency",
  signature(PairwiseComparisonMatrix = "FuzzyPairwiseComparisonMatrix"),
  definition=function(PairwiseComparisonMatrix, print.report)
  {

    violationText = .strictConsistencyMethod(PairwiseComparisonMatrix@fnModal)

    if (violationText != "") {
      if(print.report){
        cat(paste("Fuzzy comparison matrix isn't strictly consistent. These indeces violate the condition: \n", violationText, sep = ""))
        cat("\n")
      }
      return(FALSE)
    }
    else{
      if(print.report){
        cat("The fuzzy comparison matrix is strictly consistent. \n")
        cat("\n")
      }
      return(TRUE)
    }
  }
)

#' @rdname strictConsistency-methods
#' @aliases strictConsistency,PairwiseComparisonMatrix-method
setMethod(
  f="strictConsistency",
  signature(PairwiseComparisonMatrix = "PairwiseComparisonMatrix"),
  definition=function(PairwiseComparisonMatrix, print.report)
  {

    violationText = .strictConsistencyMethod(PairwiseComparisonMatrix@values)

    if (violationText != "") {
      if(print.report){
        cat(paste("Comparison matrix isn't strictly consistent. These indeces violate the condition: \n", violationText, sep = ""))
        cat("\n")
      }
      return(FALSE)
    }
    else{
      if(print.report){
        cat("The comparison matrix is strictly consistent. \n")
        cat("\n")
      }
      return(TRUE)
    }
  }
)





setGeneric(".strictConsistencyMethod",
           function(matrix) standardGeneric(".strictConsistencyMethod"))

setMethod(
  ".strictConsistencyMethod",
  signature(matrix = "matrix"),
  function (matrix)
  {
    violationText = ""
    size = nrow(matrix)

    for (i in 1:size){
      for (j in 1:size){
        for (k in 1:size){

          if(!isTRUE(all.equal(matrix[i,k], (matrix[i,j]*matrix[j,k]), tolerance = 0.0001))){

            text = paste("[",i,",",k,"] != ([",i,",",j,"]*[",j,",",k,"]) -- ",
                         matrix[i,k]," != ", matrix[i,j]*matrix[j,k],"\n", sep = "")
            violationText = paste(violationText, text, sep = "")
          }
        }
      }
    }

    return(violationText)
  }
)
