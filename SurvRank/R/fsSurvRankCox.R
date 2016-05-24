#' Cox ranking function
#'
#' This function ranks the input features according to their significance in the univariate cox models
#' @param x x matrix or data.frame
#' @param y response y as a survival object, generated with \code{Surv()}
#' @param ... other arguments, not used now
#' @keywords SurvRank
#' @export

fsSurvRankCox = function(x,y,...){
  n = ncol(x)
  score = unlist(lapply(1:n,function(i) summary(survival::coxph(y~.,data=data.frame(x[,i])))$sctest[1]))
  names(score) = colnames(x)
  res=list()
  res$rank = names(sort(score,decreasing=T))
  res$nrcoef = length(res$rank)
  return(res)
}
