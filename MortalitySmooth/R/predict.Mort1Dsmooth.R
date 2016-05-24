predict.Mort1Dsmooth <-
function(object, newdata = NULL,
                                 type=c("link", "response"),
                                 se.fit = FALSE, ...){
  ## Input:
  ## object: a Mort1Dsmooth object
  ## newdata: optionally,
  ##          a vector in which to look for x with
  ##          which to predict. If omitted,
  ##          the fitted linear predictors are used  
  ## type: the type of prediction required. The
  ##       default ("link") is on the scale of the
  ##       linear predictors;
  ##       the alternative "response" is on the scale
  ##       of the response variable.
  ## se.fit: logical switch indicating if
  ##         standard errors are required
  
  ## Output: a list with components 
  ## fit: Predictions
  ## se.fit: Estimated standard errors

  type <- match.arg(type)
  
  if (!se.fit) {
    ## No standard errors
    if(is.null(newdata)) {
      pred <- switch(type,
                     link = object$linear.predictors-
                     object$offset,
                     response = object$fitted.values)
    }else{
      ## original data
      n <- object$n
      x <- object$x
      y <- object$y
      offset <- object$offset
      ## new data
      new.x <- sort(unique(c(x, newdata)))
      new.n <- length(new.x)
      new.w <- rep(0, new.n)
      whi <- new.x%in%x
      new.w[whi] <- 1
      new.y <- rep(0, new.n)
      new.y[whi] <- y
      new.offset <- rep(100, new.n)
      new.offset[whi] <- offset
      ## new B-spline basis
      if(min(new.x)>=min(x) & max(new.x)<=max(x)){
        ## only interpolation
        new.ndx <- object$ndx
        xl <- min(new.x)
        xr <- max(new.x)
        xmax <- xr + 0.01 * (xr - xl)
        xmin <- xl - 0.01 * (xr - xl)
        new.B <- MortSmooth_bbase(x=new.x,
                                  xl=xmin, xr=xmax,
                                  ndx=new.ndx,
                                  deg=object$deg)
      }else{
        ## extrapolation + eventual interpolation
        deg <- object$deg
        xl <- min(x)
        xr <- max(x)
        xmax <- xr + 0.01 * (xr - xl)
        xmin <- xl - 0.01 * (xr - xl)
        dx <- (xmax - xmin)/object$ndx
        xl1 <- min(new.x)
        xr1 <- max(new.x)
        minK0 <- xmin-deg:(deg+100)*dx
        minK <- minK0[which(minK0<=xl1)[deg]]
        maxK0 <- xmax+deg:(deg+100)*dx
        maxK <- maxK0[which(maxK0>=xr1)[deg+1]]
        knots <- seq(minK, maxK, by=dx)
        PP <- outer(new.x, knots, MortSmooth_tpower, deg)
        nn <- dim(PP)[2]
        DD <-diff(diag(nn),diff=deg+1)/(gamma(deg+1)*dx^deg)
        new.B <- (-1)^(deg + 1) * PP %*% t(DD)
      }

      ## matplot(x, object$B, t="l", col=1,
      ##         lty=1, xlim=range(new.x), lwd=2)
      ## matlines(new.x, new.B, t="l", col=2, lty=1, lwd=3)
      
      ## penalty stuff
      nb <- ncol(new.B)
      D. <- diff(diag(nb), diff=object$pord)
      DtD <- t(D.)%*%D.
      ## initialize
      new.y[is.na(new.y)] <- 0
      ## 0) simple poisson-GLM
      fit0 <- glm(round(new.y) ~ new.x + offset(new.offset),
                  family=poisson, weights=new.w)
      ## 1) simple penalized-poisson-GLM with new.B
      eta0 <- log(fit0$fitted) - new.offset
      mu0 <- exp(new.offset + eta0)
      w0 <- new.w*mu0
      z0 <- new.w*((new.y - mu0)/mu0 + eta0)
      BtWB <- t(new.B) %*% (w0 * new.B)
      BtWz <- t(new.B) %*% (w0 * z0)
      a.init <- solve(BtWB + 1e08 * DtD, BtWz)
      ## plot(x, log(death/exposure),
      ##      ylim=c(-5, -3), xlim=range(new.x))
      ## lines(new.x, new.B%*%a.init, col=5)
      ## lines(new.x, cbind(1, new.x)%*%fit0$coef, col=4)
      
      ## estimating
      new.fit <- Mort1Dsmooth_estimate(x=new.x,
                                       y=new.y,
                                       offset=new.offset,
                                       wei=new.w,
                                       psi2=1,
                                       B=new.B,
                                       lambda=object$lambda,
                                       DtD=DtD,
                                       a.init=a.init,
                                       MON=FALSE,
                                       TOL1=1e-6,
                                       MAX.IT=50)
      new.eta0 <- c(new.B %*% new.fit$a)
      ## lines(new.x, new.eta0, col=2, lwd=2)
      new.eta <- new.eta0[new.x%in%newdata]
      ## points(newdata, new.eta, col=4)
      ## points(x, log(object$fitted/e), col=3)
      new.fitted0 <- exp(new.eta0 + new.offset)
      new.fitted0[!whi] <- NA
      new.fitted <- new.fitted0[new.x%in%newdata]
      
      pred <- switch(type,
                     link = new.eta,
                     response = new.fitted)
    }
  }else{
    ## WITH standard errors
    if(is.null(newdata)){
      ori.B <- object$B
      ori.nb <- ncol(ori.B)
      ori.BtWB<- t(ori.B)%*%(as.vector(object$fitted)*ori.B)
      ori.D <- diff(diag(ori.nb), diff=object$pord)
      ori.PtP <- object$lambda * t(ori.D) %*% ori.D
      ori.solBtWB <- solve(ori.BtWB + ori.PtP)
      ori.solBtWB1 <- ori.solBtWB%*%ori.BtWB %*% ori.solBtWB
      ori.se <- sqrt(diag(ori.B%*%ori.solBtWB1%*%t(ori.B)))
      ## check usage of overdispersion
      check.over <- eval(object$call$overdispersion)
      if(is.null(check.over)){
        check.over <- FALSE
      }
      if(check.over){
        se.inflate <- sqrt(object$psi2)
        ori.se <- ori.se * se.inflate
      }
      fit <- switch(type,
                    link = object$linear.predictors - object$offset,
                    response = object$fitted.values)
      se <- switch(type,
                   link = ori.se,
                   response = object$fitted.values * (exp(ori.se) -1))
      pred <- list(fit=fit, se.fit=se)
    }else{
      ## original data
      n <- object$n
      x <- object$x
      y <- object$y
      offset <- object$offset
      ## new data
      new.x <- sort(unique(c(x, newdata)))
      new.n <- length(new.x)
      new.w <- rep(0, new.n)
      whi <- new.x%in%x
      new.w[whi] <- 1
      new.y <- rep(0, new.n)
      new.y[whi] <- y
      new.offset <- rep(100, new.n)
      new.offset[whi] <- offset
      ## new B-spline basis
      if(min(new.x)>=min(x) & max(new.x)<=max(x)){
        ## only interpolation
        new.ndx <- object$ndx
        xl <- min(new.x)
        xr <- max(new.x)
        xmax <- xr + 0.01 * (xr - xl)
        xmin <- xl - 0.01 * (xr - xl)
        new.B <- MortSmooth_bbase(x=new.x,
                                  xl=xmin, xr=xmax,
                                  ndx=new.ndx,
                                  deg=object$deg)
      }else{
        ## extrapolation + eventual interpolation
        deg <- object$deg
        xl <- min(x)
        xr <- max(x)
        xmax <- xr + 0.01 * (xr - xl)
        xmin <- xl - 0.01 * (xr - xl)
        dx <- (xmax - xmin)/object$ndx
        xl1 <- min(new.x)
        xr1 <- max(new.x)
        minK0 <- xmin-deg:(deg+100)*dx
        minK <- minK0[which(minK0<=xl1)[deg]]
        maxK0 <- xmax+deg:(deg+100)*dx
        maxK <- maxK0[which(maxK0>=xr1)[deg+1]]
        knots <- seq(minK, maxK, by=dx)
        PP <- outer(new.x, knots, MortSmooth_tpower, deg)
        nn <- dim(PP)[2]
        DD <-diff(diag(nn),diff=deg+1)/(gamma(deg+1)*dx^deg)
        new.B <- (-1)^(deg + 1) * PP %*% t(DD)
      }

      ## penalty stuff
      nb <- ncol(new.B)
      D. <- diff(diag(nb), diff=object$pord)
      DtD <- t(D.)%*%D.
      ## initialize
      new.y[is.na(new.y)] <- 0
      ## 0) simple poisson-GLM with ages
      fit0 <- glm(round(new.y) ~ new.x + offset(new.offset),
                  family=poisson, weights=new.w)
      ## 1) simple penalized-poisson-GLM with B
      eta0 <- log(fit0$fitted) - new.offset
      mu0 <- exp(new.offset + eta0)
      w0 <- new.w*mu0
      z0 <- new.w*((new.y - mu0)/mu0 + eta0)
      BtWB <- t(new.B) %*% (w0 * new.B)
      BtWz <- t(new.B) %*% (w0 * z0)
      a.init <- solve(BtWB + 1e08 * DtD, BtWz)
      ## plot(x, log(object$fitted/e),
      ##      ylim=c(-6.5, -0.5), xlim=range(new.x))
      ## lines(new.x, new.B%*%a.init, col=5)

      ## estimating
      new.fit <- Mort1Dsmooth_estimate(x=new.x,
                                       y=new.y,
                                       offset=new.offset,
                                       wei=new.w,
                                       psi2=1,
                                       B=new.B,
                                       lambda=object$lambda,
                                       DtD=DtD,
                                       a.init=a.init,
                                       MON=FALSE,
                                       TOL1=1e-6,
                                       MAX.IT=50)
      
      new.eta0 <- c(new.B %*% new.fit$a)
      new.eta <- new.eta0[new.x%in%newdata]
      ## points(newdata, new.eta, col=2)
      
      new.fitted0 <- exp(new.eta0 + new.offset)
      new.fitted0[!whi] <- NA
      new.fitted <- new.fitted0[new.x%in%newdata]

      fit <- switch(type,
                    link = new.eta,
                    response = new.fitted)
      
      ## standard errors
      BtWB <- new.fit$BtWB
      solBtWB <- solve(new.fit$BtWB + new.fit$P)
      solBtWB1 <- solBtWB %*% BtWB %*% solBtWB

      se.eta0 <- sqrt(diag(new.B %*% solBtWB1 %*% t(new.B)))
      se.eta <- se.eta0[new.x%in%newdata]
      ## points(newdata, new.eta+2*se.eta, col=3)
      ## points(new.x, new.eta0+2*se.eta0, col=4)
      
      se.fitted0 <- new.fitted0*(exp(se.eta0)-1)
      se.fitted0[!whi] <- NA
      se.fitted <- se.fitted0[new.x%in%newdata]
      ## check usage of overdispersion
      check.over <- eval(object$call$overdispersion)
      if(is.null(check.over)){
        check.over <- FALSE
      }
      if(check.over){
        se.inflate <- sqrt(object$psi2)
        se.eta <- se.eta * se.inflate
        se.fitted <- se.fitted * se.inflate
      }
      se <- switch(type,
                   link = se.eta,
                   response = se.fitted)
      pred <- list(fit=fit, se.fit=se)

      ## plot(newdata, pred$se)
      ## abline(v=x)
    }
  }
  return(pred)
}
