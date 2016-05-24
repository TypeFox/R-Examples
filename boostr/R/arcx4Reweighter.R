
#' @title Reweighter for the arc-x4 algorithm.
#' @description An implementation of the reweighter described in the arc-x4
#' algorithm.
#' 
#' @family arc-x4
#' 
#' @template -Reweighter
#' @param m a vector length equal to \code{nrow(data)} enumerating each time the
#' \eqn{i^{th}}{ith} entry in \code{data} has been misclassified by all the estimators
#' previously built.
#' 
#' @note If you're going to use this reweighter with \code{\link{boost}} you'll want
#' to initialize \code{m} to 0 by including \code{.reweighterArgs=list(m=0)} inside
#' your \code{metadata} list.
#' 
#' @return
#' \item{m}{the updated count of misclassifications.}
#' 
#' @export 
arcx4Reweighter <- function(prediction, response, weights, m, ...) {
  d <- as.numeric(prediction != response)
  
  new_m <- m + d
  weights <- (1 + new_m^4) / sum( 1 + new_m^4 )
  
  list(weights=weights, m=new_m)
}

class(arcx4Reweighter) <- c("reweighter", class(arcx4Reweighter))