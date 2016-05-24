#' Rpart ranking function
#'
#' This function ranks the input features according to their importance in recursive partitioning and regression trees fitted with the function \code{rpart}
#' @param x x matrix or data.frame
#' @param y response y as a survival object, generated with \code{Surv()}
#' @param ... other arguments, not used now
#' @keywords SurvRank
#' @export


fsSurvRankRpart = function(x,y,...){
  res=list()
  mod.rpart = rpart::rpart(y~.,data=data.frame(y=y,x=x),method="exp",...)
  vi.mod = mod.rpart$variable.importance
  nam = names(sort(vi.mod,decreasing = T))
  res$rank = substr(x = nam,3,nchar(x))
  res$nrcoef = length(res$rank)
  return(res)
}
