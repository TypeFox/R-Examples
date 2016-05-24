varcov <- function(object, ...)
  UseMethod("varcov")

varcov.eco <- function(object, subset = NULL, ...) {
  if (is.null(subset))
    subset <- 1:nrow(object$Sigma)
  else if (max(subset) > nrow(object$Sigma))
    stop(paste("invalid input for `subset.' only", nrow(object$Sigma),
               "draws are stored.")) 

  p <- ncol(object$mu)
  n <- length(subset)
  Sigma <- array(0, c(p, p, n))
  cov <- object$Sigma[subset,]

  for (i in 1:n) {
    count <- 1
    for (j in 1:p) {
      Sigma[j,j:p,i] <- cov[i,count:(count+p-j)]
      count <- count + p - j + 1
    }
    diag(Sigma[,,i]) <- diag(Sigma[,,i]/2)
    Sigma[,,i] <- Sigma[,,i] + t(Sigma[,,i])
  }
  if (n > 1)
    return(Sigma)
  else
    return(Sigma[,,1])  
}

varcov.ecoNP <- function(object, subset = NULL, obs = NULL, ...) {
  if (is.null(subset))
    subset <- 1:nrow(object$Sigma)
  else if (max(subset) > nrow(object$Sigma))
    stop(paste("invalid input for `subset.' only", nrow(object$Sigma),
               "draws are stored.")) 

  if (is.null(obs))
    obs <- 1:dim(object$Sigma)[3]
  else if (max(subset) > dim(object$Sigma)[3])
    stop(paste("invalid input for `obs.' only", dim(object$Sigma)[3],
               "draws are stored."))
  
  p <- ncol(object$mu)
  n <- length(subset)
  m <- length(obs)
  Sigma <- array(0, c(p, p, n, m))
  cov <- object$Sigma[subset,,obs]

  for (k in 1:m) {
    for (i in 1:n) {
      count <- 1
      for (j in 1:p) {
        Sigma[j,j:p,i,k] <- cov[i,count:(count+p-j),k]
        count <- count + p - j + 1
      }
      diag(Sigma[,,i,k]) <- diag(Sigma[,,i,k]/2)
      Sigma[,,i,k] <- Sigma[,,i,k] + t(Sigma[,,i,k])
    }
  }
  if (n > 1)
    if (m > 1)
      return(Sigma)
    else
      return(Sigma[,,,1])
  else
    if (m > 1)
      return(Sigma[,,1,])
    else
      return(Sigma[,,1,1])
}
