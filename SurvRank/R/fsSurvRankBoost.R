#' Boost ranking function
#'
#' This function ranks the input features according to their selection probability in additive models via component-wise boosting
#' @param x x matrix or data.frame
#' @param y response y as a survival object, generated with \code{Surv()}
#' @param ... other arguments, not used now
#' @keywords SurvRank
#' @export
fsSurvRankBoost = function(x, y, ...){
  res = list()
  model.b = mboost::gamboost(y~.,data=data.frame(y=y,x=scale(x,center = T,scale = T)),baselearner="bols",control = mboost::boost_control(mstop = 2000),family=mboost::CoxPH())
  names.b = t(sort(summary(model.b)$selprob,decreasing=T))
  rank = c(apply(names.b,1,function(x) substr(names(x),8,nchar(names(x))-1)))
  res$rank = rank
  res$nrcoef = length(rank)
  return(res)
}
