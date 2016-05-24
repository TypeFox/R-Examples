coef.endogMNP <- function(object, subset = NULL, ...) {
  param <- object$param
#  p <- object$n.alt
#  n.cov <- ncol(param) - p*(p-1)/2
n.cov <- ncol(param) - object$n.dim	
  n <- nrow(param)
  if (is.null(subset))
    return(param[,1:n.cov])
  else if (subset > n)
    stop(paste("invalid input for `subset.' only", nrow(param), "draws are stored.")) 
  else
    return(param[subset, 1:n.cov])                 
}
