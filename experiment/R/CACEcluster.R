CACEcluster <- function(Y, D, Z, grp, data = parent.frame(),
                        match = NULL, weights = NULL, ...) {

  cov.internal <- function(x){
    return(cov(x[,1],x[,2]))
  }
  
  call <- match.call()
  Y <- eval(call$Y, envir = data)
  D <- eval(call$D, envir = data)
  Z <- eval(call$Z, envir = data)
  grp <- eval(call$grp, envir = data)
  match <- eval(call$match, envir = data)
  weights <- eval(call$weights, envir = data)
  
  ITTY <- ATEcluster(Y = Y, Z = Z, grp = grp, match = match,
                     weights = weights, ...)
  ITTD <- ATEcluster(Y = D, Z = Z, grp = grp, match = match,
                     weights = weights, ...)

  ## point estimate
  n <- ITTY$n
  ITTY.est <- ITTY$est
  ITTD.est <- ITTD$est
  CACE.est <- ITTY.est/ITTD.est

  ## outputs
  res <- list(est = CACE.est, ITTY = ITTY, ITTD = ITTD)

  ## calculation
  if (is.null(match))
    stop("the estimator is not yet available.")
  else {
    res$n1 <- n1 <- ITTY$n1
    res$n0 <- n0 <- ITTY$n0
    res$N1 <- N1 <- ITTY$N1
    res$N0 <- N0 <- ITTY$N0
    res$m <- m <- ITTY$m
    res$w <- w <- ITTY$w
    diffY <- ITTY$diff
    diffD <- ITTD$diff
    res$cov <- Cov <- m*sum((diffY*w - n*ITTY.est/m) *
                            (diffD*w - n*ITTD.est/m))/((m-1)*(n^2))
  }
  
  res$var <- CACE.var <-
    (ITTY$var*(ITTD.est^2) + ITTD$var*(ITTY.est^2) -
     2*Cov*ITTY.est*ITTD.est)/(ITTD.est^4)
  
  class(res) <- "CACEcluster"
  return(res)
}
