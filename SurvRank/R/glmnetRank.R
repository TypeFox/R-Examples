#' Ranks features of a previously fitted \code{glmnet} object
#'
#' This function ranks the input features of a previously fitted \code{glmnet} object. Ranking according to the first occurence in the lambda path
#' or effect sizes at the end of the path.
#' @param glmnet previously fitted \code{glmnet} object
#' @param first Defaults to TRUE. TRUE ranking based on occurence. FALSE based on effect sizes
#' @param names Defaults to TRUE. TRUE returns feature names. FALSE returns coefficients
#' @param ... other arguments, not used now
#' @keywords SurvRank
#' @export

glmnetRank <- function(glmnet, first=T, names=T,...) {
  beta <- glmnet$beta[, order(glmnet$lambda, decreasing=T)]
  if (first == T) {
    # Return variables that pops-up first
    coef <- unique(unlist(apply(beta, 2, function(x) {
      o <- order(abs(x), decreasing=T)
      return (o[abs(x[o]) > 0])
    })))
    nrcoef = length(coef)
    coef <- c(coef, setdiff(1:nrow(glmnet$beta), coef))
  } else {
    coef <- order(abs(beta[, ncol(beta)]), decreasing=T)
  }
  if (names) {
    # Return variable names that are non-zero at the end
    coef <- rownames(glmnet$beta)[coef]
  }
  glmnetR=list()
  glmnetR$coef = coef
  glmnetR$nrcoef = nrcoef
  return (glmnetR)
}
