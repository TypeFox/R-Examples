#*******************************************************************************
#
# Local Approximate Gaussian Process Regression
# Copyright (C) 2013, The University of Chicago
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@chicagobooth.edu)
#
#*******************************************************************************


## auto.rho:
##
## initial rho parameter chosen to balance objective values
## with (sum of) squared C values nearby to the constraint
## boundary

auto.rho <- function(obj, C, init=10)
{
  if(init > length(obj)) init <- length(obj)
  if(length(obj) != nrow(C)) stop("length(obj) != nrow(C)")

  ## completely valid runs
  v <- apply(C, 1, function(x) { all(x <= 0) })
  
  ## smallest sum of squared invalids
  Civ <- C[!v,,drop=FALSE] 
  valid <- Civ <= 0
  Civ[valid] <- 0
  Civs2 <- rowSums(Civ^2)
  cm2 <- min(Civs2)

  ## smallest valid objective restricted to first init (10)
  f <- obj[1:init]; fv <- f[v[1:init]]
  if(length(fv) > 0) fm <- min(fv)
  else fm <- median(f)
  ## or if no valids use the median of first init (10)

  rho <- cm2/(2*fm)
  return(rho)
}


## optim.auglag:
##
## Optimization of known or estimated objective under unknown constraints via
## the Augmented Lagriangian framework with GP surrogate modeling
## of the constraint functions

optim.auglag <- function(fn, B, fhat=FALSE, cknown=NULL, 
  start=10, end=100, Xstart=NULL, sep=FALSE, ab=c(3/2,4), 
  lambda=1, rho=1/2, urate=10, ncandf=function(t) { t }, 
  dg.start=c(0.1,1e-6), dlim=sqrt(ncol(B))*c(1/100,10), 
  Bscale=1, ey.tol=0.05,nomax=FALSE, N=1000, plotprog=FALSE, 
  verb=2, ...)
{
  ## check start
  if(start >= end) stop("must have start < end")

  ## check sep and determine whether to use GP or GPsep commands
  if(sep) { newM <- newGPsep; mleM <- mleGPsep; updateM <- updateGPsep; 
    alM <- alGPsep; deleteM <- deleteGPsep; nd <- nrow(B) }
  else { newM <- newGP; mleM <- mleGP; updateM <- updateGP; 
    alM <- alGP; deleteM <- deleteGP; nd <- 1 }
  formals(newM)$dK <- TRUE;

  ## get initial designwhich.min
  X <- dopt.gp(start, Xcand=lhs(10*start, B))$XX
  ## X <- lhs(start, B)
  X <- rbind(Xstart, X)
  start <- nrow(X)

  ## first run to determine dimensionality of the constraint
  out <- fn(X[1,]*Bscale, ...)
  nc <- length(out$c)

  ## allocate progress objects, and initialize
  prog <- obj <- rep(NA, start)
  C <- matrix(NA, nrow=start, ncol=nc)
  obj[1] <- out$obj; C[1,] <- out$c
  if(all(out$c <= 0)) prog[1] <- out$obj
  else prog[1] <- Inf

  ## remainder of starting run
  for(t in 2:start) { ## now that fn is vectorized we can probably remove for
    out <- fn(X[t,]*Bscale, ...)
    obj[t] <- out$obj; C[t,] <- out$c
    ## update best so far
    if(all(out$c <= 0) && out$obj < prog[t-1]) prog[t] <- out$obj
    else prog[t] <- prog[t-1]
  }

  ## handle initial lambda and rho values
  if(length(lambda) == 1) lambda <- rep(lambda, nc)
  if(is.null(rho)) rho <- auto.rho(obj, C)
  else if(length(rho) != 1 || rho <= 0) stop("rho should be a positive scalar")

  ## best auglag seen so far
  valid <- C <= 0
  if(nomax >= 0) {  ## Original AL
    Cm <- C
    Cm[valid] <- 0
    al <- obj + C %*% lambda + Cm^2 %*% rep(1/(2*rho), nc)
  } else {  ## AL with slack variables
    S <- pmax(- C - rho*matrix(rep(lambda, start), nrow=start, byrow=TRUE), 0)
    Cm <- C+S
    al <- obj + Cm %*% lambda + Cm^2 %*% rep(1/(2*rho), nc)
  }
  ybest <- min(al)
  since <- 0

  ## best valid so far
  m2 <- prog[start]

  ## initializing constraint surrogates
  Cgpi <- rep(NA, nc)
  d <- matrix(NA, nrow=nc, ncol=nd)
  Cnorm <- rep(NA, nc)
  for(j in 1:nc) {
    if(j %in% cknown) { Cnorm[j] <- 1; Cgpi[j] <- -1 }
    else {
      Cnorm[j] <- max(abs(C[,j]))
      Cgpi[j] <- newM(X, C[,j]/Cnorm[j], dg.start[1], dg.start[2])
      d[j,] <- mleM(Cgpi[j], param="d", tmin=dlim[1], tmax=dlim[2], ab=ab, 
                    verb=verb-1)$d
    }
  }
  ds <- matrix(rowMeans(d, na.rm=TRUE), nrow=1)

  ## possibly initialize objective surrogate
  if(fhat) {
    fnorm <- max(abs(obj))
    fgpi <- newM(X, obj/fnorm, dg.start[1], dg.start[2])
    df <- mleM(fgpi, param="d", tmin=dlim[1], tmax=dlim[2], ab=ab, verb=verb-1)$d
    dfs <- matrix(df, nrow=1)
  } else { fgpi <- -1; fnorm <- 1 }

  ## init for loop
  mei <- Inf
  new.params <- FALSE

  ## keeping track
  meis <- mei
  lambdas <- as.numeric(lambda)
  rhos <- rho

  ## iterating over the black box evaluations
  for(t in (start+1):end) {

    ## lambda and rho update
    valid <- C <= 0
    if(nomax >= 0) {  ## Original AL
      Cm <- C
      Cm[valid] <- 0
      al <- obj + C %*% lambda + Cm^2 %*% rep(1/(2*rho), nc)

      if(!is.finite(m2)) { ## update if no valid has been found
        
        lambda.new <- rep(0, length=ncol(C))
        rho.new <- auto.rho(obj, C)
        if(rho.new >= rho) rho.new <- 9*rho/10
      
      } else { ## update if any valid has been found

        ck <- C[which.min(al),] 
        lambda.new <- pmax(0, lambda + (1/rho) * ck)
        if(any(ck > 0)) rho.new <- rho/2
        else rho.new <- rho
      }

    } else {  ## Slack Variable AL
      S <- pmax(- C - rho*matrix(rep(lambda, t-1), nrow=t-1, byrow=TRUE), 0)
      Cm <- C+S
      al <- obj + Cm %*% lambda + Cm^2 %*% rep(1/(2*rho), nc)

      if(!is.finite(m2)) { ## update if no valid has been found
        ## same as in non-slack case above
        
        lambda.new <- rep(0, length=ncol(C))
        rho.new <- auto.rho(obj, C)
        if(rho.new >= rho) rho.new <- 9*rho/10
      
      } else { ## update if any valid has been found
        ## specific to slack variable case

        cmk <- Cm[which.min(al),]
        lambda.new <- lambda + (1/rho) * cmk
        if(mean(cmk^2) > quantile(rowMeans(Cm^2), p=0.05)) rho.new <- rho/2
        else rho.new <- rho
      }
    }

    ## printing progress  
    if(any(lambda.new != lambda) || rho.new != rho) {
      if(verb > 0) {
        cat("updating La:")
        if(rho.new != rho) cat(" rho=", rho.new, sep="")
        if(any(lambda.new != lambda))
          cat(" lambda=(", paste(signif(lambda.new,3), collapse=", "), 
            ")", sep="")
        cat("\n")
      }
      new.params <- TRUE
    } else new.params <- FALSE

    ## confirm update of augmented lagrangian
    lambda <- lambda.new; rho <- rho.new
    if(nomax >= 0) ## original
      ybest <- min(obj + C %*% lambda + Cm^2 %*% rep(1/(2*rho), nc)) ## orig
    else ybest <- min(obj + Cm %*% lambda + Cm^2 %*% rep(1/(2*rho), nc)) ## slack
    valid <- apply(valid, 1, prod)
    
    ## keep track of lambda and rho
    lambdas <- rbind(lambdas, lambda)
    rhos <- cbind(rhos, rho)

    ## rebuild surrogates periodically under new normalized responses
    if(t > (start+1) && (t %% urate == 0)) {
    
      ## constraint surrogates 
      Cnorm <- apply(abs(C), 2, max)
      for(j in 1:nc) {
        if(j %in% cknown) next;
        deleteM(Cgpi[j])
        d[j,d[j,] < dlim[1]] <- 10*dlim[1]
        d[j,d[j,] > dlim[2]] <- dlim[2]/10
        Cgpi[j] <- newM(X, C[,j]/Cnorm[j], d[j,], dg.start[2])
        d[j,] <- mleM(Cgpi[j], param="d", tmin=dlim[1], tmax=dlim[2], 
                      ab=ab, verb=verb-1)$d
      }
      ds <- rbind(ds, rowMeans(d, na.rm=TRUE))

      ## possible objective surrogate
      if(fhat) {
        deleteM(fgpi)
        fnorm <- max(abs(obj))
        fgpi <- newM(X, obj/fnorm, df, dg.start[2])
        df <- mleM(fgpi, param="d", tmin=dlim[1], tmax=dlim[2], ab=ab, 
                   verb=verb-1)$d
        dfs <- rbind(dfs, df)
      } else { df <- NULL }

      new.params <- FALSE
    }
  
    ## random candidate grid
    ncand <- ncandf(t)
    if(!is.finite(m2) || fhat) 
      XX <- lhs(ncand, B)
    else XX <- rbetter(ncand, B, sum(X[which(obj == m2),]))
    ## NOTE: might be a version of rbetter for fhat

    ## calculate composite surrogate, and evaluate EI and/or EY
    eyei <- alM(XX, fgpi, fnorm, Cgpi, Cnorm, lambda, 1/(2*rho), ybest, 
                nomax, N, fn, Bscale)
    eis <- eyei$ei; by <- "ei"
    mei <- max(eis)
    nzei <- sum(eis > 0)
    if(nzei <= ey.tol*ncand) { eis <- -(eyei$ey); by <- "ey"; mei <- Inf }
    meis <- c(meis, mei)

    ## plot progress
    if(!is.logical(plotprog) || plotprog) {
      par(mfrow=c(1,3+fhat))
      plot(prog, type="l", main="progress")
      if(is.logical(plotprog)) {
        if(length(eis) < 30) { span <- 0.5 } else { span <- 0.1 }
        g <- interp.loess(XX[,1], XX[,2], eis, span=span)
      } else g <- plotprog(XX[,1], XX[,2], eis)
      image(g, xlim=range(X[,1]), ylim=range(X[,2]), main="EI")
      if(is.matrix(valid)) valid <- apply(valid, 1, prod)
      points(X[1:start,1:2], col=valid[1:start]+3)
      points(X[-(1:start),1:2, drop=FALSE], col=valid[-(1:start)]+3, pch=19)
      matplot(ds, type="l", lty=1, main="constraint lengthscale")
      if(fhat) matplot(dfs, type="l", lty=1, main="objective lengthscale")
    }

    ## calculate next point
    m <- which.max(eis)
    xstar <- matrix(XX[m,], ncol=ncol(X))
    if(verb > 0) {
      cat("t=", t, " ", sep="")
      cat(by, "=", eis[m]/Bscale, " (", nzei,  "/", ncandf(t), ")", sep="")
      cat("; xbest=[", paste(signif(X[obj == m2,],3), collapse=" "), sep="")
      cat("]; ybest (v=", m2, ", al=", ybest, ", since=", since, ")\n", sep="")
    }

    ## new run
    out <- fn(xstar*Bscale, ...)
    ystar <- out$obj; obj <- c(obj, ystar); C <- rbind(C, out$c)

    ## update GP fits
    X <- rbind(X, xstar)
    for(j in 1:nc) if(Cgpi[j] >= 0) 
      updateM(Cgpi[j], xstar, out$c[j]/Cnorm[j], verb=verb-2)

    ## check if best valid has changed
    since <- since + 1
    valid <- apply(C, 1, function(x) { all(x <= 0) })
    if(all(out$c <= 0) && ystar < prog[length(prog)]) {
      m2 <- ystar; since <- 0 
    } ## otherwise m2 unchanged; should be the same as prog[length(prog)]
    prog <- c(prog, m2)

    ## check if best auglag has changed
    if(nomax >= 0) {
      alstar <- out$obj + lambda %*% drop(out$c)
      alstar <- alstar + rep(1/(2*rho),nc) %*% pmax(0, drop(out$c))^2
    } else {
      cps <- drop(out$c + pmax(-out$c - rho*lambda, 0))
      alstar <- out$obj + lambda %*% cps + rep(1/(2*rho),nc) %*% cps^2
    }
    ## alstar <- out$obj + lambda %*% out$c + rep(1/(2*rho),nc) %*% out$c^2
    if(alstar < ybest) { ybest <- alstar; since <- 0 }
  }

  ## delete GP surrogates
  for(j in 1:nc) if(Cgpi[j] > 0) deleteM(Cgpi[j])
  if(fhat) deleteM(fgpi)

  ## return output objects
  if(!fhat) df <- NULL
  return(list(prog=prog, mei=meis, obj=obj, X=X, C=C, d=d, df=df, 
    lambda=as.matrix(lambdas), rho=as.numeric(rhos)))
}


## optim.eic:
##
## Optimization of known or estimated objective under unknown constraints via
## the EI times the product of validity of the constraint functions, following 
## Schonlau et al 1998

optim.eic <- function(fn, B, fhat=FALSE, cknown=NULL, 
  start=10, end=100, Xstart=NULL, sep=FALSE, ab=c(3/2,4), 
  urate=10, ncandf=function(t) { t }, dg.start=c(0.1,1e-6), 
  dlim=sqrt(ncol(B))*c(1/100,10), Bscale=1, plotprog=FALSE, 
  verb=2, ...)
{
  ## check start
  if(start >= end) stop("must have start < end")

  ## check sep and determine whether to use GP or GPsep commands
  if(sep) { newM <- newGPsep; mleM <- mleGPsep; updateM <- updateGPsep; 
    eicM <- eicGPsep; deleteM <- deleteGPsep; nd <- nrow(B) }
  else { newM <- newGP; mleM <- mleGP; updateM <- updateGP; 
    eicM <- eicGP; deleteM <- deleteGP; nd <- 1 }
  formals(newM)$dK <- TRUE;

  ## get initial design
  X <- dopt.gp(start, Xcand=lhs(10*start, B))$XX
  ## X <- lhs(start, B)
  X <- rbind(Xstart, X)
  start <- nrow(X)

  ## first run to determine dimensionality of the constraint
  out <- fn(X[1,]*Bscale, ...)
  nc <- length(out$c)

  ## allocate progress objects, and initialize
  prog <- obj <- rep(NA, start)
  C <- matrix(NA, nrow=start, ncol=nc)
  obj[1] <- out$obj; C[1,] <- out$c
  if(all(out$c <= 0)) prog[1] <- out$obj
  else prog[1] <- Inf

  ## remainder of starting run
  for(t in 2:start) {
    out <- fn(X[t,]*Bscale, ...)
    obj[t] <- out$obj; C[t,] <- out$c
    ## update best so far
    if(all(out$c <= 0) && out$obj < prog[t-1]) prog[t] <- out$obj
    else prog[t] <- prog[t-1]
  }

  ## best valid so far
  m2 <- prog[start]

  ## initializing constraint surrogates
  Cgpi <- rep(NA, nc)
  d <- matrix(NA, nrow=nc, ncol=nd)
  Cnorm <- rep(NA, nc)
  for(j in 1:nc) {
    if(j %in% cknown) { Cnorm[j] <- 1; Cgpi[j] <- -1 }
    else {
      Cnorm[j] <- max(abs(C[,j]))
      Cgpi[j] <- newM(X, C[,j]/Cnorm[j], dg.start[1], dg.start[2])
      d[j,] <- mleM(Cgpi[j], param="d", tmin=dlim[1], tmax=dlim[2], ab=ab, 
                    verb=verb-1)$d
    }
  }
  ds <- matrix(rowMeans(d, na.rm=TRUE), nrow=1)

  ## possibly initialize objective surrogate
  if(fhat) {
    fnorm <- max(abs(obj))
    fgpi <- newM(X, obj/fnorm, dg.start[1], dg.start[2])
    df <- mleM(fgpi, param="d", tmin=dlim[1], tmax=dlim[2], ab=ab, verb=verb-1)$d
    dfs <- matrix(df, nrow=1)
  } else { fgpi <- -1; fnorm <- 1 }


  ## iterating over the black box evaluations
  since <- 0
  for(t in (start+1):end) {

    ## rebuild surrogates periodically under new normalized responses
    if((t > (start+1)) && (t %% urate == 0)) {
    
      ## constraint surrogates 
      Cnorm <- apply(abs(C), 2, max)
      for(j in 1:nc) {
        if(j %in% cknown) next;
        deleteM(Cgpi[j])
        d[j,d[j,] < dlim[1]] <- 10*dlim[1]
        d[j,d[j,] > dlim[2]] <- dlim[2]/10
        Cgpi[j] <- newM(X, C[,j]/Cnorm[j], d[j,], dg.start[2])
        d[j,] <- mleM(Cgpi[j], param="d", tmin=dlim[1], tmax=dlim[2], 
                      ab=ab, verb=verb-1)$d
      }
      ds <- rbind(ds, rowMeans(d, na.rm=TRUE))

      ## possible objective surrogate
      if(fhat) {
        deleteM(fgpi)
        fnorm <- max(abs(obj))
        fgpi <- newM(X, obj/fnorm, df, dg.start[2])
        df <- mleM(fgpi, param="d", tmin=dlim[1], tmax=dlim[2], ab=ab, 
                   verb=verb-1)$d
        dfs <- rbind(dfs, df)
      } else { df <- NULL }

      new.params <- FALSE
    }
  
    ## random candidate grid
    ncand <- ncandf(t)
    if(!is.finite(m2) || fhat) 
      XX <- lhs(ncand, B)
    else XX <- rbetter(ncand, B, sum(X[which(obj == m2),]))
    ## NOTE: might be a version of rbetter for fhat

    ## calculate composite surrogate, and evaluate EI and/or EY
    eic <- eicM(XX, fgpi, fnorm, Cgpi, Cnorm, m2, fn, Bscale)
    leis <- eic$lei + rowSums(eic[,-1,drop=FALSE]) 

    ## continue plotting 
    if(!is.logical(plotprog) || plotprog) {
      par(mfrow=c(1,3+fhat))
      plot(prog, type="l", main="progress")
      if(is.logical(plotprog)) {
        if(length(leis) < 30) { span <- 0.5 } else { span <- 0.1 }
        g <- interp.loess(XX[,1], XX[,2], leis, span=span)
      } else g <- plotprog(XX[,1], XX[,2], leis)
      image(g, xlim=range(X[,1]), ylim=range(X[,2]), main="EI")
      points(X[,1:2], col=valid+3, pch=19)
      matplot(ds, type="l", lty=1, main="constraint lengthscale")
      if(fhat) matplot(dfs, type="l", lty=1, main="objective lengthscale")
    }

    ## calculate next point
    m <- which.max(leis)
    xstar <- matrix(XX[m,], ncol=ncol(X))
    if(verb > 0) {
      cat("t=", t, " ", sep="")
      cat("eic=", exp(leis[m])/Bscale, sep="")
      cat("; xbest=[", paste(signif(X[obj == m2,],3), collapse=" "), sep="")
      cat("]; ybest=", m2, ", since=", since, "\n", sep="")
    }

    ## new run
    out <- fn(xstar*Bscale, ...)
    ystar <- out$obj; obj <- c(obj, ystar); C <- rbind(C, out$c)

    ## update GP fits
    X <- rbind(X, xstar)
    for(j in 1:nc) if(Cgpi[j] >= 0) 
      updateM(Cgpi[j], xstar, out$c[j]/Cnorm[j], verb=verb-2)

    ## check if best valid has changed
    since <- since + 1
    valid <- apply(C, 1, function(x) { all(x <= 0) })
    if(all(out$c <= 0) && ystar < prog[length(prog)]) {
      m2 <- ystar; since <- 0 
    } ## otherwise m2 unchanged; should be the same as prog[length(prog)]
    prog <- c(prog, m2)
  }

  ## delete GP surrogates
  for(j in 1:nc) if(Cgpi[j] > 0) deleteM(Cgpi[j])
  if(fhat) deleteM(fgpi)

  ## return output objects
  if(!fhat) df <- NULL
  return(list(prog=prog, obj=obj, X=X, C=C, d=d, df=df))
}
