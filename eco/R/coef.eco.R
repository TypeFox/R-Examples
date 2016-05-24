coef.eco <- function(object, subset = NULL, ...) {
  mu <- object$mu
  if (is.null(subset))
    subset <- 1:nrow(mu)
  else if (max(subset) > nrow(mu))
    stop(paste("invalid input for `subset.' only", nrow(mu), "draws are stored."))

  return(mu[subset,])
}
