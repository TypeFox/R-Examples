#' Random forest ranking function
#'
#' This function ranks the input features with a random forest according to the variable importance.
#' @param x x matrix or data.frame
#' @param y response y as a survival object, generated with \code{Surv()}
#' @param ... other arguments, not used now
#' @keywords SurvRank
#' @export

fsSurvRankRf <- function(x, y, ...) {
  rank <- list()
  rand.s.for = randomForestSRC::rfsrc(Surv(time,event)~.,data=data.frame(time=y[,1],event=y[,2],x),splitrule="logrankscore")
  vi = randomForestSRC::vimp(rand.s.for)
  rank$rank <- names(sort(vi$importance, decreasing=T))
  rank$nrcoef = length(rank$rank)
  return (rank)
}
