residuals.aodml <- function(object, ..., type = c("deviance", "pearson", "response")) {
  
  type <- match.arg(type)
  resp <- object$resp
  mu <- fitted(object)
  phi.scale <- object$phi.scale
  zphi <- fitted(object, what = "phi")
	zphi <- switch(phi.scale, identity = zphi, log = exp(zphi), exp = log(zphi), inverse = 1 / zphi)
	
  if(object$family == "bb") {
		m <- resp[, 1]
		n <- rowSums(resp)
		y <- m / n
  	v <- (1 / n) * mu * (1 - mu) * (1 + (n - 1) * zphi)
  }
  
  if(object$family == "nb") {
  	y <- as.vector(resp)
    v <- mu + zphi * (mu^2)
  }
  
  d <- 2 * (object$lmax - object$l)
  
  switch(
    type,
    deviance = sign(y - mu) * sqrt(d), # d should be positive
    pearson = (y - mu) / sqrt(v),
    response = y - mu
    )

}

residuals.aodql <- function(object, ...) residuals(object$fm, ...)