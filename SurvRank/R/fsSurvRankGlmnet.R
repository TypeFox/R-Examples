#' L1 (lasso) ranking function
#'
#' This function ranks the input features with the lasso algorithm in glmnet.
#' @param x x matrix or data.frame
#' @param y response y as a survival object, generated with \code{Surv()}
#' @param alp alpha value in \code{glmnet} (elasticnet mixing parameter)
#' @param ... other arguments, not used now
#' @keywords SurvRank
#' @export

fsSurvRankGlmnet <- function(x, y, alp=1,...) {
  glmnet <- glmnet::glmnet(as.matrix(x), y, family="cox", alpha=alp,...)
  fsnet = glmnetRank(glmnet,...)
  rank <- list()
  rank$rank <- fsnet$coef
  rank$nrcoef = fsnet$nrcoef
  return (rank)
}
