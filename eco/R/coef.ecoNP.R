coef.ecoNP <- function(object, subset = NULL, obs = NULL, ...) {
  mu <- object$mu
  if (is.null(subset))
    subset <- 1:nrow(mu)
  else if (max(subset) > nrow(mu))
    stop(paste("invalid input for `subset.' only", nrow(mu), "draws are stored."))

  if (is.null(obs))
    obs <- 1:dim(object$mu)[3]
  else if (max(subset) > dim(object$mu)[3])
    stop(paste("invalid input for `obs.' only", dim(object$mu)[3], "draws are stored."))
  
  return(mu[subset,,obs])
}
