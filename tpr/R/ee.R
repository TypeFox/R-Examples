tpr <- function(y, delta, x, xtv=list(), z, ztv=list(), w, tis,
                family = poisson(),
                evstr = list(link = 5, v = 3),
                alpha = NULL,
                theta = NULL,
                tidx = 1:length(tis),
                kernstr = list(kern=1, poly=1, band=range(tis)/50),
                control = list(maxit=25, tol=0.0001, smooth=0, intsmooth=0)
                ) 
#### y, delta are lgtdl
#### x and z are matrix
#### xtv and ztv are list of lgtdl
#### w and tis are vector
{
  nt <- length(tis)
  n <- length(y)
  if (!is.matrix(x)) x <- as.matrix(x)
  pv <- length(xtv)
  xmat <- cbind(x, matrix(0, nrow=n, ncol=pv))
  p <- ncol(xmat)
  if (!is.matrix(z)) z <- as.matrix(z)
  qv <- length(ztv)
  zmat <- cbind(z, matrix(0, nrow=n, ncol=qv))
  q <- ncol(zmat)
  if (is.null(alpha) || is.null(beta)) {
    alpha <- matrix(0, nt, p)
    theta <- rep(0, q)
    dat <- data.frame(yt = I(y), rt = I(delta))
    ## get initial values
    ## on a sparse grid first and then interpolate
    ## tidx <- as.integer(quantile(1:nt, p=c(0, .25, .5, .75, 1)))
    ## tidx <- as.integer(quantile(1:nt, p=seq(0, 1, by=0.1)))
    ## tidx <- 1:nt
    for (i in tidx) {
      tt <- rep(tis[i], n)
      yt <- .Call("myinterp", y, tt, PACKAGE="tpr")
      ari <- .Call("myinterp", delta, tt, PACKAGE="tpr")
      if (pv > 0) {
        for (j in 1:pv) {
          xmat[,p - pv + j] <- .Call("myinterp", xtv[[j]], tt, PACKAGE="tpr")
        }
      }
      if (qv > 0) {
        for (j in 1:qv) {
          zmat[,q - qv + j] <- .Call("myinterp", ztv[[j]], tt, PACKAGE="tpr")
        }
      }
      good <- ari == 1
      ycur <- yt[good]
      xcur <- cbind(xmat, zmat)[good, , drop = FALSE]
      ##if (i == 1) print(xcur)
      foo <- glm.fit(xcur, ycur, family=family)
      param <- coef(foo)
      if (p > 0) alpha[i,] <- param[1:p]
      if (q > 0) theta <- theta + param[(p + 1):(p + q)]
    }
    ## set up initial values
    if (q > 0) theta <- theta / length(tidx)
    if (p > 0)
      for (i in 1:p) alpha[,i] <- approx(tis[tidx], alpha[tidx,i], xout=tis)$y
  }
  ##return (list(alpha))
  ## call c function
  evstr <- lapply(evstr, as.integer)
  kernstr <- kern.str(kernstr)
  control <- tpr.control(control)
  
  err <- rep(0, nt)
  ans <- .Call("pfEst_rap", y, delta, xmat, xtv, zmat, ztv, w, tis, t(alpha), theta, err, evstr, kernstr, control) #, PACKAGE="tpr")
  names(ans) <- c("alpha", "theta", "valpha", "vtheta", "niter", "inflAlpha", "inflTheta", "delAlpha")
  ans$alpha <- t(ans$alpha)
  ans$valpha <- t(matrix(unlist(lapply(ans$valpha, diag)), nrow=p))
  ans$tis <- tis
  
  if (!is.null(dimnames(xmat)[[2]]))
    colnames(ans$alpha) <- colnames(ans$valpha) <- dimnames(xmat)[[2]]
  if (!is.null(dimnames(zmat)[[2]]))
    names(ans$theta) <- names(ans$vtheta) <- dimnames(zmat)[[2]] 
  
  class(ans) <- "tpr"
  ans
}
  
## gofTest <- function(fit1, fit2, nsim=100) {
## #### fit1 is the H0 model, reduced (constant coef) model
## #### fit2 is the full model
##   infl1 <- fit1$inflTheta[1,]
##   n <- length(infl1)
##   p <- ncol(fit2$alpha)
##   infl2 <- sapply(fit2$inflAlpha, function(x) x[p,])
##   obs <- fit2$alpha[,p] - fit1$theta[1]
##   infl <- infl2 - infl1
##   sim <- replicate(nsim, apply(infl * rnorm(n), 2, mean))
##   supobs <- max(abs(obs))
##   supsim = apply(abs(sim), 2, max)
##   list(tis = fit1$tis, obs = obs, sim = sim,
##        supobs = supobs, supsim = supsim,
##        suppval = sum(supobs >= supsim) / nsim)
## }

