predict.eco <- function(object, newdraw = NULL, subset = NULL,
                        verbose = FALSE, ...){

  if (is.null(newdraw) && is.null(object$mu))
    stop("Posterior draws of mu and Sigma must be supplied")
  else if (!is.null(newdraw)){
    if (is.null(newdraw$mu) && is.null(newdraw$Sigma))
      stop("Posterior draws of both mu and Sigma must be supplied.")
    object <- newdraw
  }

  mu <- coef(object, subset = subset)
  n.draws <- nrow(mu)
  p <- ncol(mu)
  Sigma <- varcov(object, subset = subset)
  
  Wstar <- matrix(NA, nrow=n.draws, ncol=p)
  tmp <- floor(n.draws/10)
  inc <- 1
  for (i in 1:n.draws) {
    Wstar[i,] <- mvrnorm(1, mu = mu[i,], Sigma = Sigma[,,i])
    if (i == inc*tmp & verbose) {
      cat("", inc*10, "percent done.\n")
      inc <- inc + 1
    }
  }
  res <- apply(Wstar, 2, invlogit)
  if (ncol(res) == 2)
    colnames(res) <- c("W1", "W2")
  else # this is called from predict.ecoX
    colnames(res) <- c("W1", "W2", "X")
  class(res) <- c("predict.eco", "matrix")
  return(res)
}
