## file csampling/R/csampling.R, v 1.2-2 2014-03-31
##
##  Copyright (C) 2000-2014 Alessandra R. Brazzale 
##
##  This file is part of the "csampling" package for R.  This program  
##  is free software; you can redistribute it and/or modify it under 
##  the terms of the GNU General Public License as published by the 
##  Free Software Foundation; either version 2 of the License, or (at 
##  your option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
##  MA 02111-1307 USA or look up the web page 
##  http://www.gnu.org/copyleft/gpl.html.
##
##  Please send any comments, suggestions or errors found to:
##  Alessandra R. Brazzale, Department of Statistics, University of
##  Padova, Via C. Battisti 241/243, 35121 Padova (PD), Italy.
##  Email: alessandra.brazzale@unipd.it 
##  Web: http://www.stat.unipd.it/~brazzale

make.sample.data <- function(rsmObject)
{
  m <- match.call()
  nas <- is.na(coef(rsmObject))
  wzero <- (rsmObject$weights == 0)
  cf <- coef(rsmObject)[!nas]
  anc <- rsmObject$resid[!wzero]
  X <- (if( is.null(rsmObject$X) ) model.matrix(rsmObject)
        else rsmObject$X)[!wzero, !nas]
  if( is.null(dim(X)) )  X <- as.matrix(X)
  disp <- rsmObject$dispersion
  family <- rsmObject$family
  data <- list( anc = anc, X = X, coef = cf, disp = disp, 
                family = family, fixed = rsmObject$fixed )
  data
}

rsm.sample <- function(data = stop("no data given"), R = 10000, 
                       ran.gen = stop("candidate distribution is missing, with no default"), 
                       trace = TRUE, step = 100, ...)
{
  m <- match.call()
  if(R < 3)
    stop("at least 3 iterations needed")
  if( is.null(data$anc) || is.null(data$X) || is.null(data$coef) || 
      is.null(data$family) )
    stop("data not complete")
  anc <- data$anc
  X <- data$X
  beta <- data$coef
  ncf <- length(beta)
  sigma <- if(!is.null(data$disp)) data$disp else 1
  fixed <- if(!is.null(data$fixed)) data$fixed else 
           if(!is.null(data$disp)) FALSE
  dist <- data$family
  if( (length(anc) != dim(X)[1]) || (length(beta) != dim(X)[2]) )
    stop("dimension mismatch")
  if( !is.null(data$fixed) && !data$fixed && is.null(data$disp) )
    stop("no scale parameter")
  else TRUE
  simul <- matrix(nrow = R, ncol = ncf + !fixed)
  rho <- vector(mode="numeric", length=R)
#  seed <- .Random.seed  
  seed <- set.seed( trunc(runif(1, -10000, 10000)) )  
  simul.tmp <- ran.gen(data=data, R=1, ...)
  beta.s <- simul.tmp[1:ncf]
  sigma.s <- if(!fixed) simul.tmp[ncf+1] else sigma
  dens.s <- simul.tmp[length(simul.tmp)]
  arg <- (X %*% (beta.s - beta) + as.vector(anc) * sigma.s)/sigma
  dens <- switch(dist[[1]],
                 logWeibull     = matrix(exp(arg) * 
                                       dweibull(exp(arg), shape = 1)),
                 logExponential = matrix(exp(arg) * 
                                       dweibull(exp(arg), shape = 1)),
                 logRayleigh    = matrix(exp(arg) * 
                                       dweibull(exp(arg), shape = 1)),
                 extreme        = matrix(exp(-arg) * 
                                      dweibull(exp(-arg), shape = 1)),
                 student        = matrix(dt(arg, df = dist$df)),
                 logistic       = matrix(dlogis(arg)),
                 Huber          = dHuber(arg, k = dist$k),
                 matrix(do.call(paste("d", dist[[1]], sep = ""), 
                                      list(arg))))
  dens <- apply(dens, 2, "prod")
  dens <- dens * if(!fixed) 
                   sigma.s^(length(anc)-ncf-1)/sigma^(length(anc)) 
                 else 1
  w0 <- dens/dens.s
  rho[1] <- 1
  simul[1,] <- simul.tmp[-length(simul.tmp)]
  for(i in 2:R)
  {
    simul.tmp <- ran.gen(data=data, R=1, ...) 
    beta.s <- simul.tmp[1:ncf]
    sigma.s <- if(!fixed) simul.tmp[ncf+1] else sigma
    dens.s <- simul.tmp[length(simul.tmp)]
    arg <- (X %*% (beta.s - beta) + as.vector(anc) * sigma.s)/sigma
    dens <- switch(dist[[1]],
                   logWeibull     = matrix(exp(arg) * 
                                       dweibull(exp(arg), shape = 1)),
                   logExponential = matrix(exp(arg) * 
                                       dweibull(exp(arg), shape = 1)),
                   logRayleigh    = matrix(exp(arg) * 
                                       dweibull(exp(arg), shape = 1)),
                   extreme        = matrix(exp(-arg) * 
                                      dweibull(exp(-arg), shape = 1)),
                   student        = matrix(dt(arg, df = dist$df)),
                   logistic       = matrix(dlogis(arg)),
                   Huber          = dHuber(arg, k = dist$k),
                   matrix(do.call(paste("d", dist[[1]], sep = ""), 
                                  list(arg))))
    dens <- apply(dens, 2, "prod")
    dens <- dens * if(!fixed) 
                 sigma.s^(length(anc)-ncf-1)/sigma^(length(anc)) 
                   else 1
    w1 <- dens/dens.s
    rho.tmp <- min(1, exp(log(w1)-log(w0)))
    if(is.na(rho.tmp))  rho.tmp <- 0
    rho[i] <- rho.tmp
    ru <- runif(1)
    if(ru < rho.tmp)
    {
      simul[i,] <- simul.tmp[-length(simul.tmp)]
      w0 <- w1
    }
    else  simul[i,] <- simul[i-1,]		  
    if(trace && (i%%step == 0))
    print(i)
  }
  simul <- list(sim = simul, rho = rho, seed = seed, data = data, 
                call = m)
  attr(simul, "class") <- "cs"
  simul
}

Laplace <- function(which = stop("no choice made"), 
                    data = stop("data are missing"),
                    val1, idx1, val2, idx2, log.scale = TRUE)
{
  if( is.null(data$anc) || is.null(data$X) || is.null(data$coef) || 
      is.null(data$disp) || is.null(data$family) ) 
    stop("data not complete")
  if( (length(data$anc) != dim(data$X)[1]) || 
      (length(data$coef) != dim(data$X)[2]) )
    stop("dimension mismatch")
  if( (!missing(val1) && (length(val1) < 3)) || 
      (!missing(val2) && (length(val2) < 3)) )
    stop("at last 3 values needed")
  if( (!missing(idx1) && !any((1:length(data$coef))==idx1)) ||
      (!missing(idx2) && !any((1:length(data$coef))==idx2)) )
    stop("index out of range")   
  fixed <- if( !is.null(data$fixed) ) data$fixed else FALSE
  ww <- as.character(substitute(which))
  plot.pos <- switch(ww, 
                     "s"  = Laplace.s(val1, data, 
                                      log.scale = log.scale),
                     "c"  = Laplace.c(val1, idx1, data, 
                                      fixed = fixed),
                     "cs" = Laplace.cs(val1, val2, idx1, data, 
                                       log.scale = log.scale),
                     "cc" = Laplace.cc(val1, val2, idx1, idx2, data, 
                                       fixed),
                     stop("\n method not implemented"))
  ret <- switch(ww, 
                "s"  = if(log.scale) list( ls = log(val1) ) 
                       else list( s = val1 ),
                "c"  = list( b = val1 ),
                "cs" = if(log.scale) list( b = val1, 
                                          ls = log(val2) ) 
                       else list( b = val1, s = val2 ),
                "cc" = list( b1 = val1, b2 = val2 )) 
  ret <- c(ret, list(dens=plot.pos))
  class(ret) <- if( (ww=="s") || (ww=="c") ) "Lapl.spl" 
                else "Lapl.cont"
  ret
}

Laplace.s <- function(val, data, log.scale = TRUE)      
{
  anc <- data$anc ; X <- data$X ; n <- length(anc) 
  beta <- data$coef ; sigma <- data$disp ; p <- length(beta) 
  family <- data$family 
  g0 <- family$g0 ; g2 <- family$g2 ; df <- family$df ; k <- family$k
  ret <- matrix(0, nrow=length(val), ncol=2)
  ret[,1] <- val
  for( i in seq(along=val) )
  {
    min.beta <- optim(par=beta, fn=h.s, log.s=log(val[i]), cf=beta, 
                      disp=sigma, anc=anc, X=X, family=family,
                      method="BFGS")$par
    arg <- ( X%*%(min.beta-beta)+val[i]*anc ) / sigma
    h.2 <- t(X) %*% diag(as.vector(g2(arg, df=df, k=k))) %*% X 
    ret[i,2] <- val[i]^(n-p) * exp( - sum(g0(arg, df=df ,k=k)) ) /
                               sqrt(det(h.2))
  }
  if(!log.scale)
    ret[,2] <- ret[,2]/val
  else
    ret[,1] <- log(ret[,1])      
  nas <- is.na(ret[,2])
  ret <- spline(ret[!nas,])  
  cnorm <- sum( (ret$y[-1] + ret$y[-length(ret$y)]) * 
                abs(diff(ret$x)) )/2
  ret$y <- ret$y/cnorm 
  ret
}

h.s <- function( x, log.s, cf, disp, anc, X, family)  
{
  arg <- ( X %*% (x-cf) + exp(log.s) * anc ) / disp
  g0 <- family$g0 ; df <- family$df ; k <- family$k
  sum( g0(arg, df=df, k=k) )
}

Laplace.c <- function(val, idx, data, fixed = FALSE)
{
  anc <- data$anc ; X <- data$X ; n <- length(anc) 
  beta <- data$coef ; sigma <- data$disp ; p <- length(beta) 
  family <- data$family
  g0 <- family$g0 ; g2 <- family$g2 ; df <- family$df ; k <- family$k
  ret <- matrix(0, nrow=length(val), ncol=2)
  ret[,1] <- val
  for( i in seq(along=val) )
  {
    if( fixed & (dim(X)[2]==1) )
    {
#      arg <- (X%*%(val[i]-beta) + disp*anc) / disp
      arg <- (X%*%(val[i]-beta) + sigma*anc) / sigma
      ret[i,2] <- exp( - sum( g0(arg, df=df, k=k)) )
    }
    else
    {
      p.start <- beta[-idx]
      if(!fixed) p.start <- c(p.start,log(sigma))
      min.par <- optim(par=p.start, fn=h.c, bb=val[i], idx=idx, 
                       cf=beta, disp=sigma, anc=anc, X=X, 
                       family=family, fixed=fixed, method="BFGS")$par
      min.beta <- c( if(idx!=1) min.par[1:(idx-1)], val[i], 
                     if(idx!=p) min.par[idx:(p-1)] )
      min.sigma <- if(!fixed) exp(min.par[p]) else sigma
      arg <- ( X%*% (min.beta-beta) + min.sigma*anc) / sigma
      X.t <- X[,-idx] 
      if(!fixed) X.t <- cbind(X.t, min.sigma*anc)
      h.2 <- t(X.t) %*% diag(as.vector(g2(arg, df=df, k=k))) %*% X.t
      if(!fixed) h.2[p,p] <- h.2[p,p] + (n-p)*sigma^2
      ret[i,2] <- exp( - sum(g0(arg, df=df ,k=k)) ) * 
                     ( if(!fixed) min.sigma^(n-p) else 1 ) /
                     sqrt(det(h.2))
    }
  }
  nas <- is.na(ret[,2])
  ret <- spline(ret[!nas,])
  cnorm <- sum( (ret$y[-1] + ret$y[-length(ret$y)]) * 
                abs(diff(ret$x)) )/2
  ret$y <- ret$y/cnorm 
  ret
}

h.c <- function(x, bb, idx, cf, disp, anc, X, family, fixed)
{
  p <- length(cf) ; n <- length(anc)
  cf.temp <- c( if(idx!=1) x[1:(idx-1)], bb, if(idx!=p) x[idx:(p-1)] )
  arg <- ( X %*% (cf.temp-cf) + ( if(fixed) disp 
                                  else exp(x[p]) )*anc) / disp
  g0 <- family$g0 ; df <- family$df ; k <- family$k
  sum( g0(arg, df=df, k=k) ) - (if(!fixed) (n-p)*x[p] else 0)
}

Laplace.cs <- function(val1, val2, idx1, data, log.scale = TRUE)
{
  anc <- data$anc ; X <- data$X ; n <- length(anc) 
  beta <- data$coef ; sigma <- data$disp ; p <- length(beta) 
  family <- data$family
  g0 <- family$g0 ; g2 <- family$g2 ; df <- family$df ; k <- family$k
  ret <- matrix(0, nrow=length(val1), ncol=length(val2))
  for( i in seq(along=val1) )
    for( j in seq(along=val2) )
    {
      if( dim(X)[2]==1 )
      {
        arg <- ( X%*%(val1[i]-beta) + val2[j]*anc )/sigma
        ret[i,j] <-  val2[j]^(n-p) * exp( - sum(g0(arg, df=df, k=k)) ) 
      }
      else
      { 
        min.beta <- optim(par=beta[-idx1], fn=h.cs, bb=val1[i], 
                          log.s=log(val2[j]), idx=idx1, cf=beta, 
                          disp=sigma, anc=anc, X=X, family=family,
                          method="BFGS")$par
        beta.temp <- c( if(idx1!=1) min.beta[1:(idx1-1)], val1[i], 
                        if(idx1!=p) min.beta[idx1:(p-1)] )
        arg <- ( X %*% (beta.temp-beta) + val2[j] * anc ) / sigma
        X.temp <- X[,-idx1]
        G2 <- as.vector(g2(arg, df=df, k=k))
        G2 <- ifelse( is.na(G2), 0, G2 )
        h.2 <- t(X.temp) %*% diag(G2) %*% X.temp
        ret[i,j] <-  val2[j]^(n-p) * exp( -sum(g0(arg, df=df, k=k)) ) *
                                     1/sqrt(det(h.2))
      }
    }
  if(!log.scale)
       ret <- ret / matrix(rep(val2,length(val1)), ncol=length(val2), 
                           byrow=TRUE)
  ret[is.na(ret)] <- 0
  ret[is.infinite(ret)] <- 0
  temp <- ret[-1,]+ret[-length(val1),]
  temp <- (temp[,-1]+temp[,length(val2)])/4              
  supp <- as.vector(abs(diff(val1))) %*%               
            t( as.vector(abs(diff(if(log.scale) log(val2) 
                                  else val2))) )
  cnorm <- sum(temp * supp)
  ret <- ret/cnorm
  ret
}

h.cs <- function(x, bb, log.s, idx, cf, disp, anc, X, family)  
{
  p <- length(cf)
  cf.temp <- c( if(idx!=1) x[1:(idx-1)], bb, if(idx!=p) x[idx:(p-1)] )
  arg <- ( X %*% (cf.temp-cf) + exp(log.s)*anc )/ disp
  g0 <- family$g0 ; df <- family$df ; k <- family$k
  sum( g0(arg, df=df, k=k) )
}

Laplace.cc <- function(val1, val2, idx1, idx2, data, fixed)
{
  anc <- data$anc ; X <- data$X ; n <- length(anc) 
  beta <- data$coef ; sigma <- data$disp ; p <- length(beta) 
  family <- data$family
  g0 <- family$g0 ; g2 <- family$g2 ; df <- family$df ; k <- family$k
  ret <- matrix(0, nrow=length(val1), ncol=length(val2))
  for( i in seq(along=val1) )
    for( j in seq(along=val2) )
    { 
      if( fixed & (dim(X)[2]==2) )
      {
        val.x <- beta
        val.x[idx1] <- val1[i] ; val.x[idx2] <- val2[j]
        arg <- ( X%*%(val.x-beta) + sigma*anc )/sigma
        ret[i,j] <-  exp( - sum(g0(arg, df=df, k=k)) ) 
      }
      else
      {   
        p.start <- beta[-c(idx1,idx2)]
        if(!fixed) p.start <- c(p.start,log(sigma))
        min.par <- optim(par=p.start, fn=h.cc, bb=c(val1[i], val2[j]),
                         idx1=idx1, idx2=idx2, cf=beta, disp=sigma, 
                         anc=anc, X=X, family=family, fixed=fixed,
                         method="BFGS")$par
        min.beta <- c( if(idx1!=1) min.par[1:(idx1-1)], val1[i],
                       if(idx2!=(idx1+1)) min.par[idx1:(idx2-2)],
                       val2[j], if(idx2!=p) min.par[(idx2-1):(p-2)])
        min.sigma <- if(!fixed) exp(min.par[p-1]) else sigma
        arg <- ( X%*% (min.beta-beta) + min.sigma*anc) / sigma
        X.t <- X[,-c(idx1,idx2)] 
        if(!fixed) X.t <- cbind(X.t, min.sigma*anc)
        G2 <- as.vector(g2(arg, df=df, k=k))
        G2 <- ifelse( is.na(G2), 0, G2 )
        h.2 <- t(X.t) %*% diag(G2) %*% X.t
        if(!fixed) h.2[p-1,p-1] <- h.2[p-1,p-1] + (n-p)*sigma^2
        ret[i,j] <-  exp( - sum(g0(arg, df=df ,k=k)) ) * 
                        ( if(!fixed) min.sigma^(n-p) else 1 ) /
                        1/sqrt(det(h.2))
      }
    }
  ret[is.na(ret)] <- 0
  ret[is.infinite(ret)] <- 0
  temp <- ret[-1,]+ret[-length(val1),]
  temp <- (temp[,-1]+temp[,length(val2)])/4              
  supp <- as.vector(abs(diff(val1))) %*% 
              t( as.vector(abs(diff(val2))) ) 
  cnorm <- sum(temp * supp)
  ret <- ret/cnorm
  ret
}

h.cc <- function(x, bb, idx1, idx2, cf, disp, anc, X, family, fixed)
{
  p <- length(cf) ; n <- length(anc)
  cf.temp <- c( if(idx1!=1) x[1:(idx1-1)], bb[1], 
                if(idx2!=(idx1+1)) x[idx1:(idx2-2)],
                bb[2], if(idx2!=p) x[(idx2-1):(p-2)])
  arg <- ( X %*% (cf.temp-cf) + ( if(fixed) disp 
                                  else exp(x[p-1]) )*anc) / disp
  g0 <- family$g0 ; df <- family$df ; k <- family$k
  sum( g0(arg, df=df, k=k) ) - (if(!fixed) (n-p)*x[p-1] else 0)
}

plot.Lapl.spl <- function(x, ...)
    invisible(plot(x$dens, xlab = switch(names(x)[1], 
                                           "b" = "beta", 
                                           "s" = "scale", 
                                           "ls" = "log scale"),
                             ylab = "marginal density", ...))

plot.Lapl.cont <- function(x, ...)
    invisible(contour(x[[1]], x[[2]], x$dens, 
            xlab = if( names(x)[1]=="b" ) "beta" else "beta 1",
            ylab = switch(names(x)[2], 
                          "s" = "scale", "ls" = "log scale", 
                          "b2" = "beta 2"), ...))

rmt <- function(n, df = stop("\'df\' argument is missing, with no default"), 
                mm = rep(0, mult), cov = diag(rep(1, mult)), mult, 
                is.chol = FALSE)
{
  if( !missing(cov) && 
      (dim(as.matrix(cov))[1] != dim(as.matrix(cov))[2]) )
    stop("covariance matrix not square")
  if( !missing(mm) && !missing(cov) && 
      (length(mm) != dim(as.matrix(cov))[2]) )
    stop("size mismatch for mean vector and covariance matrix")
  if(missing(mult)) 
  {
    mult <- if(!missing(mm)) length(mm) else 
              if(!missing(cov)) dim(as.matrix(cov))[2] else 1
  }
  else 
  {
    if(!missing(mm) && (length(mm) != mult))
      stop("size mismatch")
    if(!missing(cov) && (dim(as.matrix(cov))[2] != mult))
      stop("size mismatch")
  }
  if(mult == 1)
    return(sqrt(cov) * rt(n, df = df) + mm)
  S <- if(!is.chol) t(chol(cov)) else cov
  x <- matrix(rnorm(mult * n), nrow = mult, ncol = n, byrow = FALSE)
  y <- rchisq(n, df)
  a <- if(length(n) > 1) 
	 S %*% x %*% diag(as.vector(sqrt(df/y))) + 
	   matrix(c(rep(mm, n)), nrow = mult, ncol = n, byrow = FALSE)
       else
	 S %*% x * sqrt(df/y) + 
           matrix(c(rep(mm, n)), nrow = mult, ncol = n, byrow = FALSE)
  names(a) <- c()
  a
}

dmt <-  function(x, df = stop("\'df\' argument is missing, with no default"),
                 mm = rep(0,length(x)), cov = diag(rep(1,length(x)))) 
{
  m <- length(x)
  d <- gamma((df+m)/2)/gamma(df/2)/(pi*df)^(m/2)/
         sqrt(det(cov))*
           (1+t((x-mm))%*%solve(qr(cov))%*%(x-mm)/df)^(-(df+m)/2)
  d
}
