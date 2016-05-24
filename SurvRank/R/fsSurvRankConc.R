#' Concordance ranking function
#'
#' This function ranks the input features with the concordance measure.
#' @param x x matrix or data.frame
#' @param y response y as a survival object, generated with \code{Surv()}
#' @param ... other arguments, not used now
#' @keywords SurvRank
#' @export

fsSurvRankConc = function(x,y,...){
  res=list()
  n = ncol(x)
  score=NULL
  for (i in 1:n){
    fit = survival::coxph(y~x[,i])
    sc = survival::survConcordance(y~predict(fit),x)$concordance
    score = c(score,sc)
  }
  names(score) = colnames(x)
  res$rank = names(sort(score,decreasing=T))
  res$nrcoef = length(res$rank)
  return(res)
}
