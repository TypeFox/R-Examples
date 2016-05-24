copas <- function(x,
                  gamma0.range=NULL,
                  gamma1.range=NULL,
                  ngrid=20,
                  nlevels=10,
                  levels=NULL,
                  slope=NULL,
                  left=NULL,
                  rho.bound=0.9999,
                  sign.rsb=0.1,
                  backtransf=x$backtransf,
                  silent=TRUE,
                  warn=options()$warn){
  
  meta:::chkclass(x, "meta")
  
  if (!is.numeric(rho.bound) && (rho.bound<=0|rho.bound>1))
    stop("no valid value for 'rho.bound'")
  
  if (!is.null(slope)){
    if (!is.numeric(slope))
      stop("Argument 'slope' must be numeric")
    if (length(slope)>1){
      warning(paste("Argument 'slope' must be of length 1;",
                    "first element of vector is used"))
      slope <- slope[1]
    }
  }
  
  
  ## Check significance level for test of residual selection bias
  ##
  meta:::chklevel(sign.rsb)
  
  
  oldopt <- options(warn=warn)
  on.exit(options(oldopt))
  
  
  TE <- x$TE
  seTE <- x$seTE
  sel <- !is.na(TE) & !is.na(seTE)
  if (length(TE) != sum(sel))
    warning(paste(length(TE) - sum(sel),
                  "observation(s) dropped due to missing values"))
  TE <- TE[sel]
  seTE <- seTE[sel]
  ##
  TE.random <- x$TE.random
  seTE.random <- x$seTE.random
  ##
  tau <- x$tau
  ##
  seTE.min <- min(seTE)
  seTE.max <- max(seTE)
  
  
  if (!silent){
    cat("\n\n")
    cat("====================================\n")
    cat("========== COPAS ANALYSIS ==========\n")
    cat("====================================\n\n")
    ##
    ## print:
    ## (a) fixed effect analysis
    ## (b) random effect analysis
    ## (c) test for heterogeneity, using appropriate function
    ##  
    cat("\n")
    cat("1) Summary statistics and test for heterogeneity\n")
    cat("================================================\n")
    print(summary(x))
  }
  
  
  ## calculate and display the range of gamma0,
  ## gamma1 used subsequently,
  ## unless options given by user override these
  ##
  if (is.null(gamma1.range)){
    gamma1.range <- c(0, 1.29/(1/seTE.min - 1/seTE.max))
  }
  if (is.null(gamma0.range))
    gamma0.range <- c(-0.25 - gamma1.range[2]/seTE.max, 2)
  ##
  ##
  gamma0 <- seq(gamma0.range[1], gamma0.range[2], length=ngrid)
  gamma1 <- seq(gamma1.range[1], gamma1.range[2], length=ngrid)
  ##
  gamma0.min <- min(gamma0)
  gamma0.max <- max(gamma0)
  gamma1.min <- min(gamma1)
  gamma1.max <- max(gamma1)
  ##
  ##
  if (!silent){
    ##
    cat("\n")
    cat("2) Ranges of gamma0, gamma1 used for calculating contour plot\n")
    cat("=============================================================\n")
    cat("gamma0 ranges from ",
        round(gamma0.range[1], 2),
        " to ",
        round(gamma0.range[2],2),
        "\n")
    cat("gamma1 ranges from ",
        round(gamma1.range[1], 2),
        " to ",
        round(gamma1.range[2],2),
        "\n")
    ##
    ##
    publprob.seTE.max <- range(pnorm(gamma0+gamma1/seTE.max))
    publprob.seTE.min <- range(pnorm(gamma0+gamma1/seTE.min))
    ##
    cat(paste("Range of probability publishing trial with largest SE:  (",
              round(publprob.seTE.max[1], 3),", ",
              round(publprob.seTE.max[2], 3),")",sep="") ,"\n")
    cat(paste("Range of probability publishing trial with smallest SE: (",
              round(publprob.seTE.min[1], 3),", ",
              round(publprob.seTE.min[2], 3),")",sep="") ,"\n\n")
  }
  
  
  ## calculate the contour plot
  ##
  if (!silent){
    cat("\n")
    cat("3) Starting calculations for contour plot\n")
    cat("=========================================\n")
  }
  ##
  if (is.null(left))
    left <- as.logical(sign(metabias(x, meth="linreg", k.min=3)$estimate[1])==1)
  ##
  if (left)
    rho0 <-  rho.bound/2
  else
    rho0 <- -rho.bound/2
  
  
  ##
  TE.contour <- matrix(NA, nrow=ngrid, ncol=ngrid)
  ##
  for (i in seq(along=gamma0)){
    for (j in seq(along=gamma1)){
      ##
      try(junk0 <- optim(c(TE.random, rho0, tau),
                         fn=copas.loglik.without.beta,
                         gr=copas.gradient.without.beta,
                         lower=c(-Inf, -rho.bound,   0),
                         upper=c( Inf,  rho.bound, Inf),
                         gamma=c(gamma0[i], gamma1[j]),
                         TE=TE, seTE=seTE,
                         method="L-BFGS-B"),
          silent=TRUE)
      ##
      TE.contour[i,j] <- junk0$par[1]
    }
    if (!silent){
      cat(paste(round(100*i*j/(ngrid*ngrid), 0), "%, ",sep="" ))
    }
  }
  ##
  if (!silent){
    cat("Done\n\n")
  }
  
  
  ## calculations to approximate route of orthogonal line
  ##
  if (!silent){
    cat("\n")
    cat("4) Calculating approximate orthogonal line\n")
    cat("==========================================\n")
  }
  ##
  ## the approach of lowess smoothing the smallest distances
  ## gives very bendy curves: try instead to calculate gradient of
  ## each contour, then average (-1/.) and then draw straight line through
  ## top right value.
  ##
  gamma0.rescale <- (gamma0-gamma0.min) / (gamma0.max-gamma0.min)
  gamma1.rescale <- (gamma1-gamma1.min) / (gamma1.max-gamma1.min)
  ##
  if (is.null(levels))
    junk <- contourLines(x=gamma0.rescale,
                         y=gamma1.rescale,
                         z=TE.contour,
                         nlevels=nlevels)
  else
    junk <- contourLines(x=gamma0.rescale,
                         y=gamma1.rescale,
                         z=TE.contour,
                         levels=levels)
  ##
  levels <- as.numeric(unlist(junk)[names(unlist(junk))=="level"])
  nlevels <- length(levels)
  ##
  nobs <- rep(NA, nlevels)
  adj.r.squareds <- rep(NA, nlevels)
  slopes <- rep(NA, nlevels)
  se.slopes <- rep(NA, nlevels)
  intercepts <- rep(NA, nlevels)
  ##
  for(l in 1:(nlevels)){
    lm.op <- lm(junk[[l]]$y ~ junk[[l]]$x)
    nobs[l] <- length(junk[[l]]$x)
    adj.r.squareds[l] <- summary(lm.op)$adj.r.squared
    slopes[l] <- lm.op$coef[2]
    se.slopes[l] <- sqrt(diag(vcov(lm.op))[2])
    intercepts[l] <- lm.op$coef[1]
  }

  ## calculate crossings of contours with approximate orthogonal line
  ## this provides the points to estimate the effect
  ##
  adj.r.squareds[is.nan(adj.r.squareds)] <- -100
  ##
  sel <- adj.r.squareds > 0
  ##
  if (is.null(slope)){
    ##
    if (all(!sel)){
      warning("No contour line with corresponding adjusted r.square larger than zero")
      slope <- NA
    }
    else
      slope <- -1/metagen(slopes[sel],
                          1/sqrt(nobs)[sel])$TE.fixed
  }
  ##
  ##x.slope <- ((1-slope-intercepts )/(slopes-slope))[sel]
  ##
  x.slope <- ((1-slope-intercepts )/(slopes-slope))
  y.slope <- 1+slope*(x.slope-1)
  
  
  if (!silent){
    cat("Done\n\n")
  }
  
  
  ## calculations for plot how mean/se changes as
  ## come down from the maximum
  ##
  if (!silent){
    cat("\n")
    cat("5) Calculating TE and seTE, as come down slope\n")
    cat("==============================================\n")
  }
  ##
  gamma0.slope <- x.slope*(gamma0.max-gamma0.min) + gamma0.min
  gamma1.slope <- y.slope*(gamma1.max-gamma1.min) + gamma1.min
  ##
  ## Select only those values within the range of gamma0 and gamma1
  ##
  sel0 <- (gamma0.slope >= min(gamma0.range) &
           gamma0.slope <= max(gamma0.range))
  sel1 <- (gamma1.slope >= min(gamma1.range) &
           gamma1.slope <= max(gamma1.range))
  ##
  x.slope <- x.slope[sel0&sel1]
  y.slope <- y.slope[sel0&sel1]
  ##
  gamma0.slope <- gamma0.slope[sel0&sel1]
  gamma1.slope <- gamma1.slope[sel0&sel1]
  ##
  ## Reorder x.slope, y.slope, gamma0.slope, gamma1.slope
  ## according to publprob.seTE.max
  ##
  ord <- rev(order(pnorm(gamma0.slope+gamma1.slope/seTE.max)))
  ##
  x.slope      <- x.slope[ord]
  y.slope      <- y.slope[ord]
  gamma0.slope <- gamma0.slope[ord]
  gamma1.slope <- gamma1.slope[ord]
  ##
  ## Add the point (gamma0=10, gamma1=0) representing the usual random
  ## effects model with no allowance for selection
  ## (Copas, Shi, 2001, p.256)
  ##
  gamma0.slope <- c(10, gamma0.slope)
  gamma1.slope <- c( 0, gamma1.slope)
  ##
  ##
  sel2 <- !is.na(x.slope) & !is.na(y.slope)
  sel3 <- !is.na(gamma0.slope) & !is.na(gamma1.slope)
  ##
  x.slope <- x.slope[sel2]
  y.slope <- y.slope[sel2]
  ##
  gamma0.slope <- gamma0.slope[sel3]
  gamma1.slope <- gamma1.slope[sel3]
  ##
  n.gamma0.slope <- length(gamma0.slope)
  ##
  publprob <- pnorm(gamma0.slope+gamma1.slope/seTE.max)
  ##
  ##
  ##
  TE.slope   <- rep(NA, n.gamma0.slope)
  seTE.slope <- rep(NA, n.gamma0.slope)
  rho.slope  <- rep(NA, n.gamma0.slope)
  tau.slope  <- rep(NA, n.gamma0.slope)
  beta.slope <- rep(NA, n.gamma0.slope)
  ##
  loglik1  <- rep(NA, n.gamma0.slope)
  conv1    <- rep(NA, n.gamma0.slope)
  message1 <- rep("", n.gamma0.slope)
  ##
  for (i in seq(along=gamma1.slope)){
    try(junk1 <- optim(c(TE.random, rho0, tau),
                       fn=copas.loglik.without.beta,
                       gr=copas.gradient.without.beta,
                       lower=c(-Inf, -rho.bound,   0),
                       upper=c( Inf,  rho.bound, Inf),
                       gamma=c(gamma0.slope[i], gamma1.slope[i]),
                       TE=TE, seTE=seTE,
                       method="L-BFGS-B"),
        silent=TRUE)
    ##
    TE.slope[i]  <- junk1$par[1]
    rho.slope[i] <- junk1$par[2]
    tau.slope[i] <- junk1$par[3]
    ##
    loglik1[i]  <- junk1$value
    conv1[i]    <- junk1$convergence
    message1[i] <- junk1$message
    ##
    ## in case of singular hessian, do the best we can:
    ##
    try(junk2 <- optim(c(TE.random, rho0, tau),
                       fn=copas.loglik.without.beta,
                       gr=copas.gradient.without.beta,
                       lower=c(-Inf, -rho.bound,   0),
                       upper=c( Inf,  rho.bound, Inf),
                       gamma=c(gamma0.slope[i], gamma1.slope[i]),
                       TE=TE, seTE=seTE,
                       method="L-BFGS-B", hessian=TRUE),
        silent=TRUE)
    ##
    ##print(junk2$hessian)
    ##
    try(seTE.slope[i] <-
        sqrt(solve(junk2$hessian+0.00000001)[1,1]),
        silent=TRUE)
    ##
    ## if this fails, take previous sd
    ##
    if ((i>1 & is.na(seTE.slope[i])) ||
        (i>1 & seTE.slope[i]==0))
      seTE.slope[i] <- seTE.slope[i-1]
    ##
    ## if that fails, try this!
    ##
    if (is.na(seTE.slope[i]) || seTE.slope[i]==0)
      try(seTE.slope[i] <- sqrt(1/junk2$hessian[1,1]),
          silent=TRUE)
  }
  
  
  if (!silent){
    cat("Done\n\n")
  }
  
  
  ## calculations for goodness of fit plot along orthogonal line
  ## (above), calculate log-likelihood in model containing sd
  ##
  if (!silent){  
    cat("\n")
    cat("6) Calculating goodness of fit, as come down orthogonal line\n")
    cat("============================================================\n")
  }
  ##
  if (left)
    rho.lim <- c(0, rho.bound)
  else
    rho.lim <- c(-rho.bound, 0)
  ##
  ## beta constrained:
  ##
  TE.slope.bc   <- rep(NA, n.gamma0.slope)
  rho.slope.bc  <- rep(NA, n.gamma0.slope)
  tau.slope.bc  <- rep(NA, n.gamma0.slope)
  beta.slope.bc <- rep(NA, n.gamma0.slope)
  ##
  loglik2  <- rep(NA, n.gamma0.slope)
  conv2    <- rep(NA, n.gamma0.slope)
  message2 <- rep("", n.gamma0.slope)
  ##
  for (i in seq(along=gamma1.slope)){
    try(junk3 <- optim(c(TE.random, rho0, tau, 0),
                       fn=copas.loglik.with.beta,
                       lower=c(-Inf, rho.lim[1], 0  , -Inf),
                       upper=c( Inf, rho.lim[2], Inf,  Inf),
                       gamma=c(gamma0.slope[i], gamma1.slope[i]),
                       TE=TE, seTE=seTE,
                       method="L-BFGS-B"),
        silent=TRUE)
    ##
    TE.slope.bc[i]   <- try(junk3$par[1])
    rho.slope.bc[i]  <- try(junk3$par[2])
    tau.slope.bc[i]  <- try(junk3$par[3])
    beta.slope.bc[i] <- try(junk3$par[4])
    ##
    loglik2[i]  <- try(junk3$value)
    conv2[i]    <- try(junk3$convergence)
    message2[i] <- try(junk3$message)
  }
  ##
  ## P-value for test of residual selection bias:
  ##
  pval.rsb <- 1-pchisq(2*(loglik1-loglik2), df=1)

  
  if (!silent){
    cat("Done\n\n")
  }
  
  
  ##
  ## compute no of unpublished studies
  ##
  N.unpubl <- rep(NA, length(publprob))
  ##
  for (i in seq(along=N.unpubl)){
    p.si <- pnorm(gamma0.slope[i]+gamma1.slope[i]/seTE)
    N.unpubl[i] <- sum((1-p.si)/p.si)
  }
  
  
  res <- list(TE=TE,
              seTE=seTE,
              TE.random=TE.random,
              seTE.random=seTE.random,
              left=left,
              rho.bound=rho.bound,
              gamma0.range=gamma0.range,
              gamma1.range=gamma1.range,
              slope=slope,
              regr=list(
                levels=levels,
                nobs=nobs,
                adj.r.squareds=adj.r.squareds,
                slopes=slopes,
                se.slopes=se.slopes,
                intercepts=intercepts),
              ngrid=ngrid,
              nlevels=nlevels,
              gamma0=gamma0,
              gamma1=gamma1,
              TE.contour=TE.contour,
              ##
              x.slope=x.slope,
              y.slope=y.slope,
              ##
              TE.slope=TE.slope,
              seTE.slope=seTE.slope,
              rho.slope=rho.slope,
              tau.slope=tau.slope,
              ##
              loglik1=loglik1,
              conv1=conv1,
              message1=message1,
              ##
              ##TE.slope.bc=TE.slope.bc,
              ##rho.slope.bc=rho.slope.bc,
              ##tau.slope.bc=tau.slope.bc,
              ##beta.slope.bc=beta.slope.bc,
              ##
              loglik2=loglik2,
              conv2=conv2,
              message2=message2,
              ##
              publprob=publprob,
              pval.rsb=pval.rsb,
              sign.rsb=sign.rsb,
              N.unpubl=N.unpubl,
              sm=x$sm, call=match.call(), x=x)
  
  res$version <- utils::packageDescription("metasens")$Version
  
  if (!is.null(x$title))
    res$title <- x$title
  if (!is.null(x$complab))
    res$complab <- x$complab
  if (!is.null(x$outclab))
    res$outclab <- x$outclab
  
  res$backtransf <- backtransf
  
  class(res) <- c("copas")
  
  res
}
