#' Function to assess Weak Consistency of Comparison Matrix
#'
#' @description Check if for \eqn{a_{ij}>1,a_{jk}>1}{a[i,j]>1,a[j,k]>1} applies that
#' \deqn{a_{ik}>=\max(a_{ij},a_{jk})}{a[i,k]>=max(a[i,j],a[j,k])} for all \eqn{i,j,k = 1,2,\dots,n}, where
#' \eqn{n} is size of \eqn{a}.
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
#' @rdname weakConsistency-methods
#' @name weakConsistency
setGeneric("weakConsistency",
           function(PairwiseComparisonMatrix, print.report = TRUE) standardGeneric("weakConsistency"))

#' @rdname weakConsistency-methods
#' @aliases weakConsistency,FuzzyPairwiseComparisonMatrix-method
setMethod(
  f="weakConsistency",
  signature(PairwiseComparisonMatrix = "FuzzyPairwiseComparisonMatrix"),
  definition=function(PairwiseComparisonMatrix, print.report)
  {

    violationText = .weakConsistencyMethod(PairwiseComparisonMatrix@fnModal)

    if (violationText != "") {
      if(print.report){
        cat(paste("Fuzzy comparison matrix isn't consistent. These indeces violate the condition: \n", violationText, sep = ""))
        cat("\n")
      }
      return(FALSE)
    }
    else{
      if(print.report){
        cat("The fuzzy comparison matrix is weakly consistent. \n")
        cat("\n")
      }
      return(TRUE)
    }
  }
)

#' @rdname weakConsistency-methods
#' @aliases weakConsistency,PairwiseComparisonMatrix-method
setMethod(
  f="weakConsistency",
  signature(PairwiseComparisonMatrix = "PairwiseComparisonMatrix"),
  definition=function(PairwiseComparisonMatrix, print.report)
  {

    violationText = .weakConsistencyMethod(PairwiseComparisonMatrix@values)

    if (violationText != "") {
      if(print.report){
        cat(paste("Comparison matrix isn't consistent. These indeces violate the condition: \n", violationText, sep = ""))
        cat("\n")
      }
      return(FALSE)
    }
    else{
      if(print.report){
        cat("The comparison matrix is weakly consistent. \n")
        cat("\n")
      }
      return(TRUE)
    }
  }
)





setGeneric(".weakConsistencyMethod",
           function(matrix) standardGeneric(".weakConsistencyMethod"))

setMethod(
  ".weakConsistencyMethod",
  signature(matrix = "matrix"),
  function (matrix)
  {
    violationText = ""
    size = nrow(matrix)

    for (i in 1:size){
      for (j in 1:size){
        for (k in 1:size){

          if (matrix[i,j]>1 && matrix[j,k]>1){
            if(matrix[i,k] < max(matrix[i,j], matrix[j,k])){

              text = paste("[",i,",",k,"] < max([",i,",",j,"],[",j,",",k,"]) -- ",
                           matrix[i,k]," < max(",
                           matrix[i,j],",",
                           matrix[j,k],")","\n", sep = "")
              violationText = paste(violationText, text, sep = "")
            }
          }

        }
      }
    }

    return(violationText)
  }
)
