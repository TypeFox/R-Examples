ATEcluster <- function(Y, Z, grp, data = parent.frame(),
                       match = NULL, weights = NULL, fpc = TRUE) {

  call <- match.call()
  Y <- eval(call$Y, envir = data)
  Z <- eval(call$Z, envir = data)
  grp <- eval(call$grp, envir = data)
  match <- eval(call$match, envir = data)
  weights <- eval(call$weights, envir = data)

  n <- length(Y)
  res <- list(call = call, n = n, Y = Y, Z = Z, grp = grp,
              match = match, weights = weights) 
  if (is.null(match))
    stop("This option is not yet available.")
  else {
    res$m <- m <- length(unique(match))
    res$Y1bar <- Y1bar <- tapply(Y[Z==1], match[Z==1], mean)
    res$Y0bar <- Y0bar <- tapply(Y[Z==0], match[Z==0], mean)
    res$diff <- diff <- Y1bar-Y0bar
    res$n1 <- n1 <- tapply(rep(1, sum(Z==1)), match[Z==1], sum)
    res$n0 <- n0 <- tapply(rep(1, sum(Z==0)), match[Z==0], sum)
  }

  if (is.null(weights)) {
    ## variance for PATE1 (sampling of clusters)
    N1 <- w1 <- n1
    N0 <- w0 <- n0
  } else {
    ## variance for PATE2 (double sampling)
    w1 <- N1 <- tapply(weights[Z==1], match[Z==1], mean)
    w0 <- N0 <- tapply(weights[Z==0], match[Z==0], mean)
  }
  w <- w1 + w0
  w <- n*w/sum(w)
  ## estimates
  ATE.est <- weighted.mean(diff, w)
  ATE.var <- m*sum((w*diff-n*ATE.est/m)^2)/((m-1)*(n^2))
  ## donner&klar methods:
  w.dk <- w1*w0/(w1 + w0)
  w.dk <- n*w.dk/sum(w.dk)
  ATEdk.est <- weighted.mean(diff, w.dk)
  ATEdk.var <- sum(w.dk^2)*sum(w.dk*(diff-ATEdk.est)^2)/(n^3)
  ATE.dkvar <- sum(w^2)*sum(w*(diff-ATE.est)^2)/(n^3)
  ## lower bound for CATE variance
  if (!is.null(weights)) {
    Y1var <- tapply(Y[Z==1], match[Z==1], var)/n1
    Y0var <- tapply(Y[Z==0], match[Z==0], var)/n0
    if (fpc) {
      Y1var <- (1-n1/N1)*Y1var
      Y0var <- (1-n0/N0)*Y0var
      if ((sum(n0 > N0)+sum(n1 > N1))>0)
        stop("population size is smaller than sample size")
    }
    res$Y1var <- Y1var
    res$Y0var <- Y0var
    res$var.lb <- sum((w/n)^2*(Y1var+Y0var))
  }
  ## unbiased estimation
  ##Y1sum <- tapply(Y[Z==1], match[Z==1], sum)
  ##Y0sum <- tapply(Y[Z==0], match[Z==0], sum)
  ##ATE.estU <- 2*(sum(Y1sum)-sum(Y0sum))/n
  ##ATE.varU <- 4*m*var(Y1sum-Y0sum)/(n^2)
  
  ## return the resutls
  res$est <- ATE.est
  res$est.dk <- ATEdk.est
  res$var <- ATE.var
  res$dkvar <- ATE.dkvar
  res$var.dk <- ATEdk.var
  res$eff <- 1/(1-2*cov(w*Y1bar, w*Y0bar)/(var(w*Y1bar)+var(w*Y0bar)))
  res$w <- w
  res$w.dk <- w.dk
  if (!is.null(match)) {
    res$w1 <- w1
    res$w0 <- w0
    res$N0 <- N0
    res$N1 <- N1
  }
  class(res) <- "ATEcluster"
  return(res)
}
