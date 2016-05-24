#  file sn/R/sn-funct.R  (various functions)
#  This file is a component of the package 'sn' for R 
#  copyright (C) 1997-2015 Adelchi Azzalini
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#---------
dsn <- function(x, xi=0, omega=1, alpha=0, tau=0, dp=NULL, log=FALSE)
{
  if(!is.null(dp)) {
    if(!missing(alpha)) 
      stop("You cannot set both 'dp' and component parameters") 
    xi <- dp[1]
    omega <- dp[2]
    alpha <- dp[3]
    tau <- if(length(dp)>3) dp[4] else 0
    }
  z <- (x-xi)/omega
  logN <- (-log(sqrt(2*pi)) -logb(omega) -z^2/2)
  if(abs(alpha) < Inf)   
    logS <- pnorm(tau * sqrt(1+alpha^2) + alpha*z, log.p=TRUE)
  else
    logS  <- log(as.numeric(sign(alpha)*z + tau > 0)) 
  logPDF <- as.numeric(logN + logS - pnorm(tau, log.p=TRUE))
  replace(logPDF, omega<= 0, NaN)
  if(log) logPDF else exp(logPDF)
}

psn <- function(x, xi=0, omega=1, alpha=0, tau=0, dp=NULL, engine, ...)
{
  if(!is.null(dp)) {
    if(!missing(alpha)) 
      stop("You cannot set both 'dp' and component parameters")
    xi <- dp[1]
    omega <- dp[2]
    alpha <- dp[3]
    tau <- if(length(dp)>3) dp[4] else 0L
   }
  z <- as.numeric((x-xi)/omega)
  nz <- length(z)
  na <- length(alpha)
  if(missing(engine)) engine <- 
    if(na == 1 & nz > 3 & all(alpha*z > -5) & (tau == 0L)) 
      "T.Owen" else "biv.nt.prob"
  if(engine == "T.Owen") {
    if(tau != 0 | na > 1) 
      stop("engine='T.Owen' not compatible with other arguments")
    p <- pnorm(z) - 2 * T.Owen(z, alpha, ...)
    }
  else{ #  engine="biv.nt.prob"
    p <- numeric(nz)
    alpha <- cbind(z, alpha)[,2]
    delta <- delta.etc(alpha)
    p.tau <- pnorm(tau) 
    for(k in seq_len(nz)) {
      if(abs(alpha[k]) == Inf){
       p[k] <- if(alpha[k] > 0)
             (pnorm(pmax(z[k],-tau)) - pnorm(-tau))/p.tau
           else
             1- (pnorm(tau) - pnorm(pmin(z[k], tau)))/p.tau
      }
    else { # SNbook: formula (2.48), p.40
      R <- matrix(c(1, -delta[k], -delta[k], 1), 2, 2)
      p[k]<- mnormt::biv.nt.prob(0, rep(-Inf,2), c(z[k], tau), c(0, 0), R)/p.tau
      }
    }}
  p <- pmin(1, pmax(0, as.numeric(p)))
  replace(p, omega <= 0, NaN)
}

#
qsn <- function(p, xi = 0, omega = 1, alpha = 0, tau=0, dp=NULL, tol = 1e-08, 
         solver="NR", ...)            
{ if(!is.null(dp)) {
    if(!missing(alpha)) 
        stop("You cannot set both 'dp' and component parameters")
      xi <- dp[1]
      omega <- dp[2]
      alpha <- dp[3]
      tau <- if(length(dp) > 3) dp[4] else 0
      }
  max.q <- sqrt(qchisq(p,1)) + tau
  min.q <- -sqrt(qchisq(1-p,1)) + tau
  if(tau == 0) {
    if(alpha == Inf)  return(as.numeric(xi + omega * max.q))
    if(alpha == -Inf) return(as.numeric(xi + omega * min.q))
    }
  na <- is.na(p) | (p < 0) | (p > 1)
  zero <- (p == 0)
  one <- (p == 1)
  p <- replace(p, (na | zero | one), 0.5)
  dp0 <- c(0, 1, alpha, tau)
  if(solver == "NR") {
    dp0 <- c(0, 1, alpha, tau)
    cum <- sn.cumulants(dp=dp0, n=4)
    g1 <- cum[3]/cum[2]^(3/2)
    g2 <- cum[4]/cum[2]^2
    x <- qnorm(p)
    x <- (x + (x^2 - 1) * g1/6 + x * (x^2 - 3) * g2/24 -
          x * (2 * x^2 - 5) * g1^2/36)
    x <- cum[1] + sqrt(cum[2]) * x
    px <- psn(x, dp=dp0, ...)
    max.err <- 1
    while (max.err > tol) { # cat("qsn:", x, "\n")
      # cat('x, px:', format(c(x,px)),"\n")
      x1 <- x - (px - p)/dsn(x, dp=dp0)
      # x1 <- pmin(x1,max.q)
      # x1 <- pmax(x1,min.q)
      x <- x1
      px <- psn(x, dp=dp0, ...)
      max.err <- max(abs(px-p))
      if(is.na(max.err)) stop('failed convergence, try with solver="RFB"')
    }
    x <- replace(x, na, NA)
    x <- replace(x, zero, -Inf)
    x <- replace(x, one, Inf)
    q <- as.numeric(xi + omega * x)
  } else { if(solver == "RFB") {
	  abs.alpha <- abs(alpha)
	  if(alpha < 0) p <- (1-p)
	  x <- xa <- xb <- xc <- fa <- fb <- fc <- rep(NA, length(p))
	  nc <- rep(TRUE, length(p)) # not converged (yet)
	  nc[(na| zero| one)] <- FALSE
	  fc[!nc] <- 0
	  xa[nc] <- qnorm(p[nc])
	  xb[nc] <- sqrt(qchisq(p[nc], 1)) + abs(tau) 
	  fa[nc] <- psn(xa[nc], 0, 1, abs.alpha, tau, ...) - p[nc]
	  fb[nc] <- psn(xb[nc], 0, 1, abs.alpha, tau, ...) - p[nc]
	  regula.falsi <- FALSE
	  while (sum(nc) > 0) { # alternate regula falsi/bisection
		xc[nc] <- if(regula.falsi) 
		   xb[nc] - fb[nc] * (xb[nc] - xa[nc])/(fb[nc] - fa[nc])    else
		   (xb[nc] + xa[nc])/2
		fc[nc] <- psn(xc[nc], 0, 1, abs.alpha, tau, ...) - p[nc]
		pos <- (fc[nc] > 0)
		xa[nc][!pos] <- xc[nc][!pos]
		fa[nc][!pos] <- fc[nc][!pos]
		xb[nc][pos] <- xc[nc][pos]
		fb[nc][pos] <- fc[nc][pos]
		x[nc] <- xc[nc]
		nc[(abs(fc) < tol)] <- FALSE
		regula.falsi <- !regula.falsi 
		}
	  # x <- replace(x, na, NA)
	  x <- replace(x, zero, -Inf)
	  x <- replace(x, one, Inf)
	  Sign <- function(x) sign(x) + as.numeric(x==0)
	  q <- as.numeric(xi + omega * Sign(alpha)* x)
  } else stop("unknown solver")}
  names(q) <- names(p)
  return(q)
}
#
rsn <- function(n=1, xi=0, omega=1, alpha=0, tau=0, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(alpha)) 
      stop("You cannot set both 'dp' and component parameters")        
    xi <- dp[1]
    omega <- dp[2]
    alpha <- dp[3]
    tau <- if(length(dp)>3) dp[4] else 0
    }
  if(tau == 0) {
    u1 <- rnorm(n)
    u2 <- rnorm(n)
    id <- (u2 > alpha*u1)
    u1[id] <- (-u1[id])
    z <- u1
    }
  else { # for ESN use transformation method
    delta <- alpha/sqrt(1+alpha^2)
    truncN <- qnorm(runif(n, min=pnorm(-tau), max=1))
    z <- delta * truncN + sqrt(1-delta^2) * rnorm(n)
    }
  y <- xi+omega*z
  attr(y, "family") <- "SN"
  attr(y, "parameters") <- c(xi,omega,alpha,tau)
  return(y)
  }

dmsn <- function(x, xi=rep(0,length(alpha)), Omega, alpha,
                 tau=0, dp=NULL, log=FALSE)
{
    if(!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
       stop("You cannot set both component parameters and dp")
    if(!is.null(dp)){
      if(length(dp) < 3) stop("wrong length of non-null 'dp'")
      xi <- drop(dp[[1]])
      Omega <- dp[[2]]
      alpha <- dp[[3]]
      tau <- if(length(dp) == 4) dp[[4]] else 0
    }
    if(any(abs(alpha) == Inf)) stop("Inf's in alpha are not allowed")
    d <- length(alpha)
    Omega <- matrix(Omega,d,d)
    invOmega <- pd.solve(Omega, silent=TRUE, log.det=TRUE)
    if (is.null(invOmega))  stop("Omega matrix is not positive definite")
    logDet <- attr(invOmega, "log.det")
    x <- if(is.vector(x)) matrix(x, 1, d) else data.matrix(x) 
    if (is.vector(xi)) xi <- outer(rep(1, nrow(x)), as.vector(matrix(xi,1,d)))
    if(tau == 0){
      log.const <- logb(2)
      alpha0 <- 0
      }
    else {
      log.const <- -pnorm(tau, log.p=TRUE)
      O.alpha <- cov2cor(Omega) %*% alpha
      alpha0 <- tau*sqrt(1+sum(alpha* O.alpha))
      }
    X <- t(x - xi)
    # Q <- apply((invOmega %*% X) * X, 2, sum)
    Q <- colSums((invOmega %*% X) * X)
    L <- alpha0 + as.vector(t(X/sqrt(diag(Omega))) %*% as.matrix(alpha))
    logPDF <- (log.const - 0.5 * Q + pnorm(L, log.p = TRUE)
               - 0.5 * (d * logb(2 * pi) + logDet))
    if (log) logPDF
    else exp(logPDF)
}

pmsn <- function(x, xi=rep(0,length(alpha)), Omega, alpha, tau=0, 
                 dp=NULL, ...)
{
  if(!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
       stop("You cannot set both component parameters and dp")
  if(!is.null(dp)){
    xi <- dp$xi
    Omega <- dp$Omega
    alpha <- dp$alpha
    tau <- if(is.null(dp$tau)) 0 else dp$tau
    }
  if(any(abs(alpha) == Inf)) stop("Inf's in alpha are not allowed")
  d <- length(alpha)
  Omega <- matrix(Omega, d, d) 
  omega <- sqrt(diag(Omega))
  delta_etc <- delta.etc(alpha, Omega)
  delta <- delta_etc$delta
  Ocor <- delta_etc$Omega.cor
  Obig <- matrix(rbind(c(1,-delta), cbind(-delta,Ocor)), d+1, d+1)
  x <- if (is.vector(x)) matrix(x, 1, d) else data.matrix(x)
  if (is.vector(xi)) xi <- outer(rep(1, nrow(x)), as.vector(matrix(xi,1,d)))
  z0 <- cbind(tau, t(t(x - xi))/omega) 
  mnormt::pmnorm(z0, mean=rep(0,d+1), varcov=Obig, ...)/pnorm(tau) 
}

rmsn <- function(n=1, xi=rep(0,length(alpha)), Omega, alpha, tau=0, dp=NULL)
{# generates SN_d(..) variates using the additive (=transformation) method
  # if(!(missing(alpha) & missing(Omega) & !is.null(dp)))
  #     stop("You cannot set both component parameters and dp")
  if(!is.null(dp)) {  
     dp0 <- dp  
     dp0$nu <- NULL
     if(is.null(dp0$tau)) dp0$tau <- 0 
     if(names(dp)[1] == "beta") {
        dp0[[1]] <- as.vector(dp[[1]])
        names(dp0)[1] <- "xi"
        } 
     }
  else dp0 <- list(xi=xi, Omega=Omega, alpha=alpha, tau=tau)
  if(any(abs(dp0$alpha) == Inf)) stop("Inf's in alpha are not allowed")
  lot <- dp2cpMv(dp=dp0, family="SN", aux=TRUE)
  d <- length(dp0$alpha)
  y <- matrix(rnorm(n*d), n, d) %*% chol(lot$aux$Psi) # each row is N_d(0,Psi)
  if(dp0$tau == 0)    
    truncN <- abs(rnorm(n))  
  else 
    truncN <- qnorm(runif(n, min=pnorm(-dp0$tau), max=1))
  truncN <- matrix(rep(truncN, d), ncol=d)
  delta  <- lot$aux$delta
  z <- delta * t(truncN) + sqrt(1-delta^2) * t(y)
  y <- t(dp0$xi + lot$aux$omega * z)
  attr(y, "family") <- "SN"
  attr(y, "parameters") <- dp0
  return(y)
}

#

dst <-  function (x, xi=0, omega=1, alpha=0, nu=Inf, dp=NULL, log=FALSE)
{ 
  if(!is.null(dp)) {
     if(!missing(alpha)) 
        stop("You cannot set both component parameters and dp")
     xi <- dp[1]
     omega <- dp[2]
     alpha <- dp[3]
     nu <- dp[4]
    }
  if (nu == Inf) return(dsn(x, xi, omega, alpha, log=log))
  if (nu == 1) return(dsc(x, xi, omega, alpha, log=log))
  z   <- (x - xi)/omega
  pdf <- dt(z, df=nu, log=log)
  cdf <- pt(alpha*z*sqrt((nu+1)/(z^2+nu)), df=nu+1, log.p=log)
  if(log)
    logb(2) + pdf + cdf -logb(omega)
  else
    2 * pdf * cdf / omega
}


rst <- function (n=1, xi = 0, omega = 1, alpha = 0, nu=Inf, dp=NULL)
{ 
    if(!is.null(dp)) {
     if(!missing(alpha)) 
        stop("You cannot set both component parameters and dp")
     xi <- dp[1]
     omega <- dp[2]
     alpha <- dp[3]
     nu <- dp[4]
    }
  z <- rsn(n, 0, omega, alpha)
  if(nu < Inf) {  
    v <- rchisq(n,nu)/nu
    y <- z/sqrt(v) + xi
    }
    else y <- z+xi
  attr(y, "family") <- "ST"
  attr(y, "parameters") <- c(xi,omega,alpha,nu)
  return(y)
}

pst <- function (x, xi=0, omega=1, alpha=0, nu=Inf, dp=NULL, method=0, ...) 
{     
  if(!is.null(dp)) {
    if(!missing(alpha)) 
      stop("You cannot set both component parameters and dp")
    xi <- dp[1]
    omega <- dp[2]
    alpha <- dp[3]
    nu <- dp[4]
   }
  if(length(alpha) > 1) stop("'alpha' must be a single value")  
  if(length(nu) > 1) stop("'nu' must be a single value")  
  if (nu <= 0) stop("nu must be non-negative")
  if (nu == Inf) return(psn(x, xi, omega, alpha))
  if (nu == 1) return(psc(x, xi, omega, alpha))
  int.nu <- (round(nu) == nu)
  if((method == 1 | method ==4) & !int.nu) 
    stop("selected method does not work for non-integer nu")
  ok <- !(is.na(x) | (x==Inf) | (x==-Inf))
  z <- ((x-xi)/omega)[ok]
  if(abs(alpha) == Inf) {
    z0 <- replace(z, alpha*z < 0, 0)
    p <- pf(z0^2, 1, nu)
    return(if(alpha>0) p else (1-p))
    }  
  fp <- function(v, alpha, nu, t.value) 
          psn(sqrt(v) * t.value, 0, 1, alpha) * dchisq(v * nu, nu) * nu  
  if(method == 4 || (method ==0  && int.nu && 
     (nu < (8.2 + 3.55* log(log(length(z)+1))))))
    p <- pst_int(z, 0, 1, alpha, nu)  # "method 4"
  else  {
    p <- numeric(length(z))
    for (i in seq_len(length(z))) {
     if(abs(z[i]) == Inf)
       p[i] <- (1+sign(z[i]))/2
    else {      
      if(method==1 | method == 0) 
         p[i] <- pmst(z[i], 0, matrix(1,1,1), alpha, nu, ...) # method 1
    else {
      # upper <- if(absalpha> 1) 5/absalpha + 25/(absalpha*nu) else 5+25/nu
      upper <- 10 + 50/nu
      if(method==2 || (method==0 & (z[i] < upper) ))
          p[i] <- integrate(dst, -Inf, z[i], dp=c(0,1,alpha, nu), ...)$value
          # method 2
        else 
          p[i] <- integrate(fp, 0, Inf, alpha, nu, z[i], ...)$value
          # method 3 
   }}
  }}
  pr <- rep(NA, length(x))
  pr[x==Inf] <- 1
  pr[x==-Inf] <- 0
  pr[ok] <- p
  return(pmax(0,pmin(1,pr)))
}


pst_int <- function (x, xi=0, omega=1, alpha=0, nu=Inf) 
{# Jamalizadeh, A. and Khosravi, M. and Balakrishnan, N. (2009)
    if(nu != round(nu) | nu < 1) stop("nu not integer or not positive")
  z <- (x-xi)/omega
  if(nu == 1) 
    atan(z)/pi + acos(alpha/sqrt((1+alpha^2)*(1+z^2)))/pi
    else { if(nu==2)
      0.5 - atan(alpha)/pi + (0.5 + atan(z*alpha/sqrt(2+z^2))/pi)*z/sqrt(2+z^2)
    else
      (pst_int(sqrt((nu-2)/nu)*z, 0, 1, alpha, nu-2) + 
        pst_int(sqrt(nu-1)*alpha*z/sqrt(nu+z^2), 0, 1, 0, nu-1) * z *
        exp(lgamma((nu-1)/2) +(nu/2-1)*log(nu)-0.5*log(pi)-lgamma(nu/2)
        -0.5*(nu-1)*log(nu+z^2)))
    } 
}


qst <- function (p, xi = 0, omega = 1, alpha = 0, nu=Inf, tol = 1e-8, 
                 dp = NULL, method=0, ...)
{
  if(!is.null(dp)) {
    if(!missing(alpha)) 
      stop("You cannot set both component parameters and dp")
    xi <- dp[1]
    omega <- dp[2]
    alpha <- dp[3]
    nu <- dp[4]
    }
  if(length(alpha) > 1) stop("'alpha' must be a single value")  
  if(length(nu) > 1) stop("'nu' must be a single value")  
  if (nu <= 0) stop("nu must be non-negative")  
  if (nu > 1e4) return(qsn(p, xi, omega, alpha))
  if (nu == 1) return(qsc(p, xi, omega, alpha))
  if (alpha == Inf) 
      return(xi + omega * sqrt(qf(p, 1, nu)))
  if (alpha == -Inf) 
      return(xi - omega * sqrt(qf(1 - p, 1, nu)))
  na <- is.na(p) | (p < 0) | (p > 1)
  abs.alpha <- abs(alpha)
  if(alpha < 0) p <- (1-p)
  zero <- (p == 0)
  one <- (p == 1)
  x <- xa <- xb <- xc <- fa <- fb <- fc <- rep(NA, length(p))
  nc <- rep(TRUE, length(p)) # not converged (yet)
  nc[(na| zero| one)] <- FALSE
  fc[!nc] <- 0
  xa[nc] <- qt(p[nc], nu)
  xb[nc] <- sqrt(qf(p[nc], 1, nu)) 
  fa[nc] <- pst(xa[nc], 0, 1, abs.alpha, nu, method=method, ...) - p[nc]
  fb[nc] <- pst(xb[nc], 0, 1, abs.alpha, nu, method=method, ...) - p[nc]
  regula.falsi <- FALSE
  while (sum(nc) > 0) { # alternate regula falsi/bisection
    xc[nc] <- if(regula.falsi) 
       xb[nc] - fb[nc] * (xb[nc] - xa[nc])/(fb[nc] - fa[nc])    else
       (xb[nc] + xa[nc])/2
    fc[nc] <- pst(xc[nc], 0, 1, abs.alpha, nu, method=method, ...) - p[nc]
    pos <- (fc[nc] > 0)
    xa[nc][!pos] <- xc[nc][!pos]
    fa[nc][!pos] <- fc[nc][!pos]
    xb[nc][pos] <- xc[nc][pos]
    fb[nc][pos] <- fc[nc][pos]
    x[nc] <- xc[nc]
    nc[(abs(fc) < tol)] <- FALSE
    regula.falsi <- !regula.falsi 
    }
  # x <- replace(x, na, NA)
  x <- replace(x, zero, -Inf)
  x <- replace(x, one, Inf)
  Sign <- function(x) sign(x) + as.numeric(x==0)
  q <- as.numeric(xi + omega * Sign(alpha)* x)
  names(q) <- names(p)
  return(q)
}
 
dmst <- function(x, xi=rep(0,length(alpha)), Omega, alpha, nu=Inf, dp=NULL,
                  log = FALSE) 
{
    if(!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
       stop("You cannot set both component parameters and dp")
    if(!is.null(dp)) {
      if(length(dp) != 4) stop("wrong length of non-null 'dp'")
      xi <- drop(dp[[1]])
      Omega <- dp[[2]]
      alpha <- dp[[3]]
      nu <- dp[[4]]
      }  
    if(any(abs(alpha) == Inf)) stop("Inf's in alpha are not allowed")
    if (nu == Inf) return(dmsn(x, xi, Omega, alpha, log = log))
    d <- length(alpha)
    Omega <- matrix(Omega, d, d)
    if(!all(Omega - t(Omega) == 0)) return(NA)
      # stop("Omega not a symmetric matrix")
    invOmega <- pd.solve(Omega, silent=TRUE, log.det=TRUE)
    if(is.null(invOmega))  return(NA)
      # stop("Omega matrix is not positive definite")
    logDet <- attr(invOmega, "log.det")
    x <- if(is.vector(x)) matrix(x, 1, d) else data.matrix(x)
    if (is.vector(xi)) xi <- outer(rep(1, nrow(x)), as.vector(matrix(xi,1,d)))
    X <- t(x - xi)
    # Q <- apply((invOmega %*% X) * X, 2, sum)
    Q <- colSums((invOmega %*% X) * X)
    L <- as.vector(t(X/sqrt(diag(Omega))) %*% as.matrix(alpha))
    if(nu < 1e4) {
      log.const <- lgamma((nu + d)/2)- lgamma(nu/2)-0.5*d*logb(nu)
      log1Q <- logb(1+Q/nu) 
      }
    else {
      log.const <- (-0.5*d*logb(2)+ log1p((d/2)*(d/2-1)/nu))
      log1Q <- log1p(Q/nu)
      }
    log.dmt <- log.const - 0.5*(d * logb(pi) + logDet + (nu + d)* log1Q) 
    log.pt <- pt(L * sqrt((nu + d)/(Q + nu)), df = nu + d, log.p = TRUE)
    logPDF <-  logb(2) + log.dmt + log.pt
    if (log) logPDF else exp(logPDF)
}

rmst <- function(n=1, xi=rep(0,length(alpha)), Omega, alpha, nu=Inf, dp=NULL)
{ 
  if(!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
       stop("You cannot set both component parameters and dp")
  if(!is.null(dp)){
      if(!is.null(dp$xi)) xi <- dp$xi
        else
      if(!is.null(dp$beta)) xi <- as.vector(dp$beta)
      Omega <- dp$Omega
      alpha <- dp$alpha
      nu <- dp$nu
     }  
  if(any(abs(alpha) == Inf)) stop("Inf's in alpha are not allowed")
  d <- length(alpha)
  x <- if(nu==Inf) 1 else rchisq(n,nu)/nu
  z <- rmsn(n, rep(0,d), Omega, alpha)
  y <- t(xi+ t(z/sqrt(x)))
  attr(y, "family") <- "ST"
  attr(y, "parameters") <- list(xi=xi, Omega=Omega, alpha=alpha, nu=nu)
  return(y)
}

pmst <- function(x, xi=rep(0,length(alpha)), Omega, alpha, nu=Inf, dp=NULL, ...)
{
  if(!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
       stop("You cannot set both component parameters and dp")
  if(!is.null(dp)){
    if(!is.null(dp$xi)) xi <- dp$xi     else
      if(!is.null(dp$beta)) xi <- as.vector(dp$beta)
      Omega <- dp$Omega
      alpha <- dp$alpha
      nu <- dp$nu
      }  
  if(!is.vector(x)) stop("x must be a vector")
  if(any(abs(alpha) == Inf)) stop("Inf's in alpha are not allowed")
  if(nu == Inf) return(pmsn(x, xi, Omega, alpha))
  d <- length(alpha)
  Omega<- matrix(Omega,d,d) 
  omega<- sqrt(diag(Omega))
  Ocor <- cov2cor(Omega)
  O.alpha <- as.vector(Ocor %*% alpha)
  delta <- O.alpha/sqrt(1 + sum(alpha*O.alpha))
  Obig <- matrix(rbind(c(1, -delta), cbind(-delta, Ocor)), d+1, d+1)
  if(nu == as.integer(nu)) {
    z0 <- c(0,(x-xi)/omega)
    if(nu < .Machine$integer.max)  
      p <- 2 * mnormt::pmt(z0, mean=rep(0,d+1), S=Obig, df=nu, ...)
    else 
      p <- 2 * mnormt::pmnorm(z0, mean=rep(0,d+1), varcov=Obig, ...)    
    }
  else {# for fractional nu, use formula in Azzalini & Capitanio (2003),
        # full-length paper, last paragraph of Section 4.2[Distr.function]) 
    z <- (x-xi)/omega
    fp <- function(v, Ocor, alpha, nu, t.value) {
            pv <-  numeric(length(v))
            for(k in seq_len(length(v))) pv[k] <- (dchisq(v[k] * nu, nu) * nu *
                 pmsn(sqrt(v[k]) * t.value, rep(0,d), Ocor, alpha) )
            pv}
    p <- integrate(fp, 0, Inf, Ocor, alpha, nu, z, ...)$value
    }
  p
}

  
dmsc <- function(x, xi=rep(0,length(alpha)), Omega, alpha, dp=NULL, 
                log = FALSE) 
{
  if(is.null(dp))  
     dp <- list(xi=xi, Omega=Omega, alpha=alpha, nu=1)
  else
     dp$nu <- 1
  dmst(x, dp=dp, log = log) 
}
  
  
pmsc <- function(x, xi=rep(0,length(alpha)), Omega, alpha, dp=NULL, ...)
{
  if(is.null(dp))  
     dp <- list(xi=xi, Omega=Omega, alpha=alpha, nu=1)
  else
     dp$nu <- 1
  pmst(x, dp=dp, ...) 
}
  
  
rmsc <- function(n=1, xi=rep(0,length(alpha)), Omega, alpha, dp=NULL)
{
  if(is.null(dp))  
     dp <- list(xi=xi, Omega=Omega, alpha=alpha, nu=1)
  else
     dp$nu <- 1
  y <- rmst(n, dp=dp) 
  attr(y, "family") <- "SC"
  attr(y, "parameters") <- dp[-4]
  return(y) 
}

dsc <- function(x, xi=0, omega=1, alpha=0, dp=NULL, log = FALSE) {
  # log.pt2 <- function(x) log1p(x/sqrt(2+x^2)) - log(2)
  if(!is.null(dp)){
     if(!missing(alpha)) 
       stop("You cannot set both 'dp' and component parameters")
    xi <- dp[1]
    omega <- dp[2]
    alpha <- dp[3]
    }
  z <- (x-xi)/omega
  logPDF <- (dcauchy(x, xi, omega, log=TRUE)
             + log1p(alpha*z/sqrt(1+z^2*(1+alpha^2))))
  if(log) logPDF else exp(logPDF)
}
 
psc <- function(x, xi=0, omega=1, alpha=0, dp=NULL) 
{# Behboodian et al. / Stat. & Prob. Letters 76 (2006) p.1490, line 2
  if(!is.null(dp)){
    if(!missing(alpha)) 
      stop("You cannot set both 'dp' and component parameters")    
    xi <- dp[1]
    omega <- dp[2]
    alpha <- dp[3]
    }
  z <- (x-xi)/omega
  delta <- if(abs(alpha)==Inf) sign(alpha) else alpha/sqrt(1+alpha^2)
  atan(z)/pi + acos(delta/sqrt(1+z^2))/pi
  }
   
qsc <- function(p, xi=0, omega=1, alpha=0, dp=NULL) 
{# Behboodian et al. / Stat. & Prob. Letters 76 (2006) p.1490, formula (4)
  if(!is.null(dp)){
    if(!missing(alpha)) 
      stop("You cannot set both 'dp' and component parameters")
    xi<- dp[1]
    omega <- dp[2]
    alpha <- dp[3]
    }
  na <- is.na(p) | (p < 0) | (p > 1)
  zero <- (p == 0)
  one <- (p == 1)
  p <- replace(p, (na | zero | one), 0.5)
  u <- (p - 0.5) * pi
  delta <- if(abs(alpha) == Inf) sign(alpha) else alpha/sqrt(1+alpha^2)
  z <- delta/cos(u) + tan(u)
  z <- replace(z, na, NA)
  z <- replace(z, zero, -Inf)
  z <- replace(z, one, Inf)
  q <- (xi + omega*z)
  names(q) <- names(p)
  return(q)
  }
  
rsc <- function(n=1, xi=0, omega=1, alpha=0, dp=NULL) {
  if(!is.null(dp)){
     if(!missing(alpha)) 
       stop("You cannot set both 'dp' and component parameters")
    xi <- dp[1]
    omega <- dp[2]
    alpha <- dp[3]
  }
  y <- xi + rsn(n, 0, omega, alpha)/abs(rnorm(n))
  attr(y, "family") <- "SC"
  attr(y, "parameters") <- c(xi, omega, alpha)
  return(y) 
}

sn.cumulants <- function(xi = 0, omega = 1, alpha = 0, tau=0, 
                         dp=NULL, n=4)
 {
   cumulants.half.norm <- function(n=4){
     n <- max(n,2)
     n <- as.integer(2*ceiling(n/2))
     half.n  <-  as.integer(n/2)
     m <- 0:(half.n-1)
     a <- sqrt(2/pi)/(gamma(m+1)*2^m*(2*m+1))
     signs <- rep(c(1, -1), half.n)[seq_len(half.n)]
     a <- as.vector(rbind(signs*a, rep(0,half.n)))
     coeff <- rep(a[1],n)
     for (k in 2:n) {
        ind <- seq_len(k-1)
        coeff[k] <- a[k] - sum(ind*coeff[ind]*a[rev(ind)]/k)
        }
     kappa <- coeff*gamma(seq_len(n)+1)
     kappa[2] <- 1 + kappa[2]
     return(kappa)
    }
  if(!is.null(dp)) {
    if(!missing(alpha)) 
      stop("You cannot set both component parameters and dp")
    dp <- c(dp,0)[1:4]
    dp <- matrix(dp, 1, ncol=length(dp))
    }
  else  dp <- cbind(xi,omega,alpha,tau)
  delta <- ifelse(abs(dp[,3])<Inf, dp[,3]/sqrt(1+dp[,3]^2), sign(dp[,3]))
  tau <- dp[,4]
  if(all(tau==0)) {
    kv <- cumulants.half.norm(n)
    if(length(kv)>n) kv <- kv[-(n+1)]
    kv[2] <- kv[2] - 1
    kappa <- outer(delta,1:n,"^") * matrix(rep(kv,nrow(dp)),ncol=n,byrow=TRUE)
    }
  else{ # ESN
    if(n>4){
       warning("n>4 not allowed with ESN distribution")
       n <- min(n, 4)
       }
    kappa <- matrix(0, nrow=length(delta), ncol=0)
    for (k in 1:n) kappa <- cbind(kappa, zeta(k,tau)*delta^k)
    }
  kappa[,2] <- kappa[,2] + 1 
  kappa <- kappa * outer(dp[,2],(1:n),"^")
  kappa[,1] <- kappa[,1] + dp[,1]
  kappa[,,drop=TRUE]
} 


zeta <- function(k, x)
{ # k integer in (0,5)
  if(k<0 | k>5 | k != round(k)) return(NULL)
  na <- is.na(x)
  x  <- replace(x,na,0)
  x2 <- x^2
  z <- switch(k+1,
            pnorm(x, log.p=TRUE) + log(2),           
            ifelse(x>(-50), exp(dnorm(x, log=TRUE) - pnorm(x, log.p=TRUE)),
                            -x/(1 -1/(x2+2) +1/((x2+2)*(x2+4)) 
                              -5/((x2+2)*(x2+4)*(x2+6))
                              +9/((x2+2)*(x2+4)*(x2+6)*(x2+8)) 
                              -129/((x2+2)*(x2+4)*(x2+6)*(x2+8)*(x2+10)) )), 
            (-zeta(1,x)*(x+zeta(1,x))),
            (-zeta(2,x)*(x+zeta(1,x)) - zeta(1,x)*(1+zeta(2,x))),
            (-zeta(3,x)*(x+2*zeta(1,x)) - 2*zeta(2,x)*(1+zeta(2,x))),
            (-zeta(4,x)*(x+2*zeta(1,x)) -zeta(3,x)*(3+4*zeta(2,x))
                 -2*zeta(2,x)*zeta(3,x)),
            NULL)
  neg.inf <- (x == -Inf)
  if(any(neg.inf))
    z <- switch(k+1,
                z,
                replace(z, neg.inf, Inf),
                replace(z, neg.inf, -1),
                replace(z, neg.inf, 0),
                replace(z, neg.inf, 0),
                replace(z, neg.inf, 0),
                NULL)
  if(k>1) z<- replace(z, x==Inf, 0)
  replace(z, na, NA)
}

st.cumulants <- function(xi=0, omega=1, alpha=0, nu=Inf, dp=NULL, n=4)
{
  if(!is.null(dp)) {
      if(!missing(alpha)) 
        stop("You cannot set both component parameters and dp")
      xi <- dp[1]
      omega <- dp[2]
      alpha <- dp[3]
      nu <- dp[4]
      }
  if(nu == Inf) return(sn.cumulants(xi, omega, alpha, n=n))
  n <- min(as.integer(n),4)
  # if(nu <= n) stop("need nu>n")
  par <- cbind(xi,omega,alpha)
  delta <- par[,3]/sqrt(1+par[,3]^2)
  cum<- matrix(NA, nrow=nrow(par), ncol=n)
  cum[,1]<- mu <- b(nu)*delta
  # if(n>1) cum[,2] <- nu/(nu-2) - mu^2
  # if(n>2) cum[,3] <- mu*(nu*(3-delta^2)/(nu-3) - 3*nu/(nu-2)+2*mu^2)
  # if(n>3) cum[,4] <- (3*nu^2/((nu-2)*(nu-4)) - 4*mu^2*nu*(3-delta^2)/(nu-3)
  #                    + 6*mu^2*nu/(nu-2)-3*mu^4)- 3*cum[,2]^2
  r <- function(nu, k1, k2) 1/(1-k2/nu) - k1/(nu-k2) # (nu-k1)/(nu-k2)
  if(n>1 & nu>2) cum[,2] <- r(nu,0,2) - mu^2
  if(n>2 & nu>3) cum[,3] <- mu*((3-delta^2)*r(nu,0,3) - 3*r(nu,0,2) + 2*mu^2)
  if(n>3 & nu>4) cum[,4] <- (3*r(nu,0,2)*r(nu,0,4) - 4*mu^2*(3-delta^2)*r(nu,0,3)
                      + 6*mu^2*r(nu,0,2)-3*mu^4) - 3*cum[,2]^2
  cum <- cum*outer(par[,2],1:n,"^")
  cum[,1] <- cum[,1]+par[,1]
  cum[,,drop=TRUE]
}
 

T.Owen <- function(h, a, jmax=50, cut.point=8)
{
 T.int <-function(h, a, jmax, cut.point)
   {
     fui <- function(h,i) (h^(2*i))/((2^i)*gamma(i+1)) 
     seriesL <- seriesH <- NULL
     i  <- 0:jmax
     low<- (h <= cut.point)
     hL <- h[low]
     hH <- h[!low]
     L  <- length(hL)
     if (L > 0) {
       b    <- outer(hL, i, fui)
       cumb <- apply(b, 1, cumsum)
       b1   <- exp(-0.5*hL^2) * t(cumb)
       matr <- matrix(1, jmax+1, L) - t(b1)
       jk   <- rep(c(1,-1), jmax)[1:(jmax+1)]/(2*i+1)
       matr <- t(matr*jk) %*%  a^(2*i+1)
       seriesL  <- (atan(a) - as.vector(matr))/(2*pi)
     }
     if (length(hH) > 0)  seriesH <- 
          atan(a)*exp(-0.5*(hH^2)*a/atan(a)) * (1+0.00868*(hH*a)^4)/(2*pi)
     series <- c(seriesL, seriesH)
     id <- c((1:length(h))[low],(1:length(h))[!low]) 
     series[id] <- series  # re-sets in original order
     series
  }
  if(!is.vector(a) | length(a)>1) stop("'a' must be a vector of length 1")
  if(!is.vector(h)) stop("'h' must be a vector")
  aa <- abs(a)    
  ah <- abs(h)
  if(is.na(aa)) stop("parameter 'a' is NA") 
  if(aa==Inf) return(sign(a)*0.5*pnorm(-ah)) # sign(a): 16.07.2007
  if(aa==0)   return(rep(0,length(h)))
  na  <- is.na(h)
  inf <- (ah == Inf)
  ah  <- replace(ah,(na|inf),0)
  if(aa <= 1)
    owen <- T.int(ah,aa,jmax,cut.point)
  else
    owen<- (0.5*pnorm(ah) + pnorm(aa*ah)*(0.5-pnorm(ah)) 
               - T.int(aa*ah,(1/aa),jmax,cut.point))
  owen <- replace(owen,na,NA)
  owen <- replace(owen,inf,0)
  return(owen*sign(a))
}
          
#=========================================================================

makeSECdistr <- function(dp, family, name, compNames)
{
  ndp <- switch(tolower(family), "sn" = 3, "esn" = 4, "st" = 4, "sc" = 3, NULL)
  if(is.null(ndp)) stop(gettextf("unknown family '%s'", family))
  family <- toupper(family)
  if(length(dp) != ndp) 
    stop(gettextf("wrong number of dp components for family '%s'", family))
  if(family == "ST") {
    nu <- as.numeric(dp[4])
    if(nu <= 0) stop("'nu' for ST family must be positive")
    if(nu == Inf) {
      warning("ST family with 'nu==Inf' is changed to SN family")
      family <- "SN"
      dp <- dp[-4]
    }}

  if(is.numeric(dp)){ # univariate distribution
    if(dp[2] <= 0) stop("omega parameter must be positive") 
    fourth <- switch(family, "SN"=NULL, "ESN"="tau", "SC"=NULL, "ST"="nu")
    names(dp) <- c("xi","omega","alpha",fourth)
    name <- if(!missing(name)) as.character(name)[1]  else 
      paste("Unnamed-", toupper(family), sep="")
    obj <- new("SECdistrUv", dp=dp, family=family, name=name)
    }
  else {if(is.list(dp)) {# multivariate distribution
    names(dp) <- rep(NULL,ndp)
    d <- length(dp[[3]])
    if(any(abs(dp[[3]]) == Inf)) stop("Inf in alpha not allowed") 
    if(length(dp[[1]]) != d) stop("mismatch of parameters size")
    Omega <- matrix(dp[[2]],d,d)  
    if(any(Omega != t(Omega))) stop("Omega matrix must be symmetric")
    if(min(eigen(Omega, symmetric=TRUE, only.values=TRUE)$values) <= 0)
      stop("Omega matrix must be positive definite")
    dp0 <- list(xi=as.vector(dp[[1]]), Omega=Omega, alpha=dp[[3]])
    name <- if(!missing(name)) as.character(name)[1]  else 
      paste("Unnamed-", toupper(family), "[d=", as.character(d), "]", sep="")
    if(family=="ST") dp0$nu <- nu
    if(family=="ESN") dp0$tau <- dp[[4]]
    if(d == 1)  warning(paste(
      "A multivariate distribution with dimension=1 is a near-oxymoron.",
      "\nConsider using a 'dp' vector to define a univariate distribution.",
      "\nHowever, I still build a multivariate distribution for you."))
    if(missing(compNames)) { compNames <-
      if(length(names(dp[[1]])) == d) names(dp[[1]]) else
        as.vector(outer("V",as.character(1:d),paste,sep=""))
      }
    else {
      if(length(compNames) != d) stop("Wrong length of 'compNames'")
      compNames <- as.character(as.vector(compNames))
      }
    names(dp0$alpha) <- names(dp0$xi) <- compNames
    dimnames(dp0$Omega) <- list(compNames, compNames)  
    obj <- new("SECdistrMv", dp=dp0, family=family, name=name, 
               compNames=compNames) } 
     else stop("'dp' must be either a numeric vector or a list")}
  obj
}

summary.SECdistrUv <- function(object, cp.type="auto", probs)
{
  cp.type <- match.arg(tolower(cp.type), c("proper", "pseudo", "auto"))
  family <- slot(object,"family")
  lc.family <- tolower(family)
  name <- slot(object,"name")
  dp  <- slot(object,"dp") 
  # op <- dp2op(dp, family)
  if(family=="ST" || family=="SC") { if(cp.type=="auto") 
    cp.type <- if(family == "SC" | dp[4] <= 4) "pseudo" else "proper" }  
  if(family=="SN" || family=="ESN") cp.type <- "proper"  
  cp <- dp2cpUv(dp, family, cp.type)
  if(is.null(cp)) stop('Stop. Consider using cp.type=="pseudo"')
  if(missing(probs)) probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  if(lc.family == "esn") lc.family <- "sn"
  q.fn <- get(paste("q",lc.family, sep=""), inherits = TRUE)
  q <- q.fn(probs, dp=dp)
  names(q) <- format(probs)
  cum <- switch(lc.family,
           "sn" = sn.cumulants(dp=dp, n=4),
           "st" = st.cumulants(dp=dp, n=4),
           rep(NA,4)
           )
  std.cum <- c(gamma1=cum[3]/cum[2]^1.5, gamma2=cum[4]/cum[2]^2)
  oct <- q.fn(p=(1:7)/8, dp=dp)
  mode <- modeSECdistrUv(dp, family)
  alpha<- as.numeric(dp[3])
  delta <- delta.etc(alpha)
  q.measures <- c(bowley=(oct[6]-2*oct[4]+oct[2])/(oct[6]-oct[2]),
                  moors=(oct[7]-oct[5]+oct[3]-oct[1])/(oct[6]-oct[2]))
  aux <- list(delta=delta, mode=mode, quantiles=q, 
              std.cum=std.cum, q.measures=q.measures)
  new("summary.SECdistrUv", dp=dp, family=family, name=name,   
      cp=cp, cp.type=cp.type, aux=aux)
}

modeSECdistr <- function(dp, family, object=NULL) 
{
  if(!is.null(object)) {
     if(!missing(dp)) stop("you cannot set both arguments dp and obj")
    obj.class <- class(object)
    if(!(obj.class %in% c("SECdistrUv", "SECdistrMv"))) 
      stop(gettextf("wrong object class: '%s'", obj.class), domain = NA)
    family <- slot(object, "family")
    dp <- slot(object, "dp")
    }  
  else {
    if(missing(family)) stop("family required")
    family <- toupper(family)
    if(!(family %in% c("SN", "ESN", "ST","SC")))
      stop(gettextf("family '%s' is not supported", family), domain = NA)
    } 
  if(is.list(dp)) modeSECdistrMv(dp, family) else modeSECdistrUv(dp, family)
}

modeSECdistrUv <- function(dp, family)
{
  if(abs(dp[3]) < .Machine$double.eps) return(as.numeric(dp[1]))
  cp <- dp2cpUv(dp, family, cp.type="auto", upto=1)
  lc.family <- tolower(family)
  if(lc.family == "esn") lc.family <- "sn"
  d.fn <- get(paste("d", lc.family, sep=""), inherits = TRUE)
  int <- c(dp[1], cp[1])
  if(abs(diff(int)) < .Machine$double.eps) return(mean(int))
  opt <- optimize(d.fn, lower=min(int), upper=max(int), maximum=TRUE, dp=dp)
  as.numeric(opt$maximum)
}


modeSECdistrMv <- function(dp, family)
{
  Omega <- dp[[2]]
  alpha <- dp[[3]]
  delta_etc <- delta.etc(alpha, Omega)
  alpha.star <- delta_etc$alpha.star
  if(alpha.star <  .Machine$double.eps) return(dp[[1]])
  lc.family <- tolower(family)
  if(lc.family == "esn") lc.family <- "sn"
  direct <- sqrt(diag(Omega)) * (delta_etc$delta/delta_etc$delta.star)
  if(lc.family == "sn") {# case SN: book (5.49), +ESN
    dp1 <- c(xi=0, omega=1, alpha=alpha.star, dp$tau) 
    mode.can <- modeSECdistrUv(dp1, family)
    mode <- as.numeric(dp[[1]] + mode.can * direct)  
  } else { #  ST, SC: book Prop. 6.2
    d.fn <- get(paste("dm", lc.family, sep=""), inherits = TRUE)
    f <- function(u, dp, direct) -d.fn(dp[[1]]+ u*direct, dp=dp, log=TRUE)
    maxM <- max(dp2cpMv(dp, family, "auto", upto=1)[[1]] - dp[[1]]/direct)
    opt <- optimize(f, lower=0, upper=maxM, dp=dp, direct=direct)
    mode <- as.numeric(dp[[1]]+ opt$minimum * direct)
  }
  return(mode)
}


summary.SECdistrMv <- function(object, cp.type="auto")
{
  cp.type <- match.arg(tolower(cp.type), c("proper", "pseudo", "auto")) 
  family <- slot(object,"family")
  name <- slot(object,"name")
  dp <- slot(object,"dp")
  # op <- dp2op(dp, family)
  if(family == "SN" || family == "ESN") cp.type <- "proper"
  if(family=="ST" || family=="SC") { if(cp.type=="auto") 
    cp.type <- if(family == "SC" || dp$nu <= 4) "pseudo" else "proper"}
  cp <- dp2cpMv(dp, family, cp.type, aux=TRUE)
  aux <- cp$aux
  if(family=="SN" | family=="SC") cp <- cp[1:3] 
  cp[["aux"]] <- NULL
  mode <- modeSECdistrMv(dp, family)
  aux0 <- list(mode=mode, delta=aux$delta,  
     alpha.star=aux$alpha.star, delta.star=aux$delta.star, mardia=aux$mardia)
  new("summary.SECdistrMv", dp=dp, family=family, name=object@name, 
      compNames=object@compNames,  cp=cp, cp.type=cp.type, aux=aux0)
}

dp2cp <- function(dp, family, object=NULL, cp.type="proper", upto=NULL)
{
  if(!is.null(object)){
    if(!missing(dp)) stop("you cannot set both arguments dp and object")
    obj.class <- class(object)
    if(!(obj.class %in% c("SECdistrUv", "SECdistrMv"))) 
      stop(gettextf("wrong object class: '%s'", obj.class), domain = NA)     
    family <- slot(object, "family")
    dp <- slot(object,"dp")
    multiv <- (obj.class == "SECdistrMv")
    }
  else{
    if(missing(family)) stop("family required")
    family <- toupper(family)
    if(!(family %in% c("SN", "ESN", "ST","SC")))
      stop(gettextf("family '%s' is not supported", family), domain = NA)
    multiv <- is.list(dp)
    }
  if(!is.null(upto)) if(upto<0 | upto>4 | upto != round(upto)) { 
      warning("unsuitable value of argument 'upto', reset to NULL")
      upto <- NULL} 
  if(multiv)
    dp2cpMv(dp, family, cp.type, upto=upto)
  else
    dp2cpUv(dp, family, cp.type, upto=upto)
}
 
dp2cpUv <- function(dp, family, cp.type="proper", upto=NULL) 
{ # internal function; works also with regression parameters included
  cp.type <- match.arg(tolower(cp.type), c("proper", "pseudo", "auto"))
  family <- toupper(family)
  if(!(family %in% c("SN", "ESN", "ST", "SC")))
    stop(gettextf("family = '%s' is not supported", family), domain = NA)
  if(family %in% c("SN","ESN")){
    if(cp.type == "pseudo") 
      warning("'cp.type=pseudo' makes no sense for SN and ESN families")
    p <- length(dp)-2-as.numeric(family=="ESN")
    omega <- dp[p+1]
    if(omega <= 0) stop("scale parameter 'omega' must be positive")
    alpha <- dp[p+2]
    tau <- if(family=="ESN") as.numeric(dp[p+3]) else 0
    delta <- if(abs(alpha) < Inf) alpha/sqrt(1+alpha^2) else sign(alpha)
    mu.Z  <- zeta(1,tau)*delta
    s.Z   <- sqrt(1+zeta(2,tau)*delta^2)
    gamma1 <- zeta(3,tau)*(delta/s.Z)^3
    sigma <- omega*s.Z
    mu    <- dp[1:p]
    mu[1] <- dp[1]+sigma*mu.Z/s.Z
    beta1 <- if(p>1) mu[2:p] else NULL
    cp    <- c(mu, sigma, gamma1, if(family=="ESN") tau else NULL)
    names(cp) <- param.names("CP", family, p, x.names=names(beta1))
    if(!is.null(upto)) cp <- cp[1:(upto+p-1)]
    }
  if(family=="ST" || family=="SC") { if(cp.type=="auto") 
    cp.type <- if(family == "SC" || dp[4] <= 4) "pseudo" else "proper" }
  if(family %in%  c("SC", "ST")) {
    fixed.nu <- if(family=="SC") 1 else NULL
    cp <- st.dp2cp(dp, cp.type, fixed.nu, jacobian=FALSE, upto=upto)
    if(is.null(cp)) {warning("no CP could be found"); return(invisible())}
    # param.type <- switch(cp.type, proper="CP", pseudo="pseudo-CP")
    # names(cp) <- param.names(param.type, family)
    } 
   return(cp)
}

dp2cpMv <- 
function(dp, family, cp.type="proper", fixed.nu=NULL, aux=FALSE, upto=NULL) 
{# internal. NB: name of cp[1] must change according to dp[1]
  cp.type <- match.arg(cp.type, c("proper", "pseudo", "auto"))
  family <- toupper(family)
  if(!(family %in% c("SN", "ESN", "ST","SC")))
    stop(gettextf("family '%s' is not supported", family), domain = NA)
  if(family %in% c("SN","ESN")){  
    if(cp.type == "pseudo") 
      warning("'cp.type=pseudo' makes no sense for SN and ESN families")
    cp <- msn.dp2cp(dp, aux=aux)
    if(!is.null(upto)) cp <- cp[1:upto]
    }
  if(family %in% c("SC","ST")){
    if(cp.type=="auto") cp.type <- 
      if(family == "SC" || dp$nu <= 4) "pseudo" else "proper"
    if(family == "SC") fixed.nu <- 1
    cp <- mst.dp2cp(dp, cp.type=cp.type, fixed.nu=fixed.nu, aux=aux, upto=upto)
    if(is.null(cp)) {warning("no CP could be found"); return(invisible())}
    }
  return(cp)
}
  
msn.dp2cp <- function(dp, aux=FALSE)
{# dp2cp for multivariate SN and ESN 
  alpha <- dp$alpha
  d <- length(alpha)
  Omega <- matrix(dp$Omega, d, d)  
  omega <- sqrt(diag(Omega))
  lot <- delta.etc(alpha, Omega)
  delta <- lot$delta
  delta.star <- lot$delta.star
  alpha.star <- lot$alpha.star
  names(delta) <- names(dp$alpha)
  tau <- if(is.null(dp$tau)) 0 else dp$tau
  mu.z  <- zeta(1, tau) * delta
  sd.z  <- sqrt(1 + zeta(2, tau) * delta^2)
  Sigma <- Omega + zeta(2,tau) * outer(omega*delta, omega*delta)
  gamma1 <- zeta(3, tau) * (delta/sd.z)^3
  if(is.vector(dp[[1]])) { 
    cp <- list(mean=dp[[1]] + mu.z*omega, var.cov=Sigma, gamma1=gamma1)
    }
  else {
    beta <- dp[[1]]  
    beta[1,] <- beta[1,] + mu.z*omega
    cp <- list(beta=beta, var.cov=Sigma, gamma1=gamma1)
  }
  if(!is.null(dp$tau)) cp$tau <- tau
  if(aux){
    lambda <- delta/sqrt(1-delta^2)
    D <- diag(sqrt(1+lambda^2), d, d)
    Ocor <- lot$Omega.cor
    Psi <- D %*% (Ocor-outer(delta,delta)) %*% D
    Psi <- (Psi + t(Psi))/2
    O.inv <- pd.solve(Omega)
    O.pcor <- -cov2cor(O.inv) 
    O.pcor[cbind(1:d, 1:d)] <- 1
    R <- force.symmetry(Ocor + zeta(2,tau)*outer(delta,delta))
    ratio2 <- delta.star^2/(1+zeta(2,tau)*delta.star^2)
    mardia <- c(gamma1M=zeta(3,tau)^2*ratio2^3, gamma2M=zeta(4,tau)*ratio2^2)
    # book: (5.74), (5.75) on p.153
    cp$aux <- list(omega=omega, cor=R, Omega.inv=O.inv, Omega.cor=Ocor, 
      Omega.pcor=O.pcor, lambda=lambda, Psi=Psi, delta=delta, lambda=lambda,
      delta.star=delta.star, alpha.star=alpha.star, mardia=mardia)
    }
  return(cp)  
}

mst.dp2cp <- function(dp, cp.type="proper", fixed.nu=NULL, symmetr=FALSE, 
   aux=FALSE, upto=NULL)
{# dp2cp for multivariate ST, returns NULL if CP not found (implicitly silent)
  nu <- if(is.null(fixed.nu)) dp$nu else fixed.nu
  if(is.null(upto)) upto <- 4L
  if((round(upto) != upto)||(upto < 1)) stop("'upto' must be positive integer")
  if(nu <= upto && (cp.type =="proper")) return(NULL)
  if(cp.type == "proper")  {
    if(nu <= upto) 
      # stop(gettextf("d.f. '%s' too small, CP is undefined", nu), domain = NA)
      return(NULL)
      a <- rep(0, upto) 
      tilde <- NULL
    } else {
      a <- (1:upto) 
      tilde <- rep("~", upto)
    }
  Omega <- dp$Omega 
  d <- ncol(Omega)
  comp.names <- colnames(dp$Omega)
  alpha <- if(symmetr) rep(0, d) else dp$alpha 
  omega <- sqrt(diag(Omega))
  lot <- delta.etc(alpha, Omega)
  delta <- lot$delta
  delta.star <- lot$delta.star
  alpha.star <- lot$alpha.star
  comp.names <- colnames(dp$Omega)
  names(delta) <- comp.names  
  mu0 <- b(nu+a[1]) * delta * omega
  names(mu0) <- comp.names
  mu.2 <- b(nu+a[2]) * delta * omega
  if(is.vector(dp[[1]])) cp <- list(mean=dp[[1]] + mu0)  else {
    beta <- dp[[1]]  
    beta[1,] <- beta[1,] + mu0
    cp <- list(beta=beta)  }
  if(upto > 1) {
    Sigma <- Omega * (nu+a[2])/(nu+a[2]-2) - outer(mu.2, mu.2)
    dimnames(Sigma) <- list(comp.names, comp.names)
    cp$var.cov <- Sigma
    }
  cp$gamma1 <- if(upto > 2 & !symmetr) st.gamma1(delta, nu+a[3]) else NULL
  cp$gamma2M <- if(upto > 3 & is.null(fixed.nu))  
      mst.gamma2M(delta.star^2, nu+a[4], d) else NULL
  names(cp) <- paste(names(cp), tilde[1:length(cp)], sep="")  
  # cp <- cp[1:length(dp1)]
  if(aux){
    mardia <- mst.mardia(delta.star^2, nu, d)
    cp$aux <- list(fixed.nu=fixed.nu, 
                omega=omega, Omega.cor=lot$Omega.cor, delta=delta,
                delta.star=delta.star, alpha.star=alpha.star, mardia=mardia)
    }
  return(cp)  
}


mst.gamma2M <- function(delta.sq, nu, d)
{ # Mardia's index of kurtosis gamma_2 for ST-d
  if(delta.sq < 0 | delta.sq >1 )  stop("delta.sq not in (0,1)")
  ifelse(nu>4, 
    {R <- b(nu)^2 * delta.sq * (nu-2)/nu
     R1R <- R/(1-R)
     (2*d*(d+2)/(nu-4) + (R/(1-R)^2)*8/((nu-3)*(nu-4))
      +2*R1R^2*(-(nu^2-4*nu+1)/((nu-3)*(nu-4))+2*(nu/((nu-3)*b(nu)^2)-1))
      +4*d*R1R/((nu-3)*(nu-4))) },
    Inf)
}

mst.mardia <- function(delta.sq, nu, d) 
{# Mardia's gamma1 and gamam2 for MST; book: (6.31), (6.32), p.178
  if(delta.sq < 0 | delta.sq > 1)  stop("delta.sq not in (0,1)")
  if(d < 1) stop("d < 1") 
  cum <- st.cumulants(0, 1, sqrt(delta.sq/(1-delta.sq)), nu)
  mu <- cum[1]
  sigma <- sqrt(cum[2])
  gamma1 <- cum[3]/sigma^3
  gamma2 <- cum[4]/sigma^4
  gamma1M <- if(nu > 3) (gamma1^2 + 3*(d-1)*mu^2/((nu-3)*sigma^2)) else Inf
  r <- function(nu, k1, k2) 1/(1 - k2/nu) - k1/(nu - k2) # (nu-k1)/(nu-k2)
  gamma2M <- if(nu > 4) (gamma2 + 3 +(d^2-1)*r(nu,2,4) +2*(d-1)*(r(nu,0,4)
                -mu^2*r(nu,1,3))/sigma^2 - d*(d+2)) else Inf
  return(c(gamma1M=gamma1M, gamma2M=gamma2M))            
}    
	
cp2dp <- function(cp, family){
  family <- toupper(family)
  if(!(family %in% c("SN", "ESN", "ST","SC")))
      stop(gettextf("family '%s' is not supported", family), domain = NA)
  dp <- if(is.list(cp))  cp2dpMv(cp, family)  else cp2dpUv(cp, family)
  if(anyNA(dp)) dp <- NULL
  return(dp)
}
 
cp2dpUv <- function(cp, family, silent=FALSE, tol=1e-8) 
{ # internal function; works also with regression parameters included
   family <- toupper(family)
   if(family=="ESN") stop("cp2dp for ESN not yet implemented")
   if(family == "SN") {     
     p <- length(cp)-2-as.numeric(family=="ESN")
     beta1 <- if (p>1) cp[2:p] else NULL
     b <- sqrt(2/pi) 
     sigma  <- cp[p+1]
     excess <- max(0, -sigma)
     gamma1 <- cp[p+2]
     tau <- if(family=="ESN") as.numeric(cp[p+3]) else 0
     max.gamma1 <- 0.5*(4-pi)*(2/(pi-2))^1.5
     if (abs(gamma1) >= max.gamma1) {
       if (silent) excess <- excess + (abs(gamma1) - max.gamma1) else 
         {message("gamma1 outside admissible range"); return(invisible())}}
     if(excess > 0) {
       out <- NA
       attr(out, "excess") <- excess
       return(out)
       }
     r  <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^(1/3)
     delta <- r/(b*sqrt(1+r^2))
     alpha <- delta/sqrt(1-delta^2)
     mu.z <- b*delta
     sd.z <- sqrt(1-mu.z^2)
     beta <- cp[1:p]
     omega <- cp[p+1]/sd.z
     beta[1] <- cp[1] - omega*mu.z
     dp <- as.numeric(c(beta, omega, alpha))
     names(dp) <- param.names("DP", family, p, x.names=names(beta1))
     return(dp)
     }
  if(family == "ST") return(st.cp2dp(cp, silent=silent, tol=tol))
  if(family == "SC") stop("this makes no sense for SC family")
  warning(gettextf("family = '%s' is not supported", family), domain = NA)
  invisible(NULL)
}

cp2dpMv <- function(cp, family, silent=FALSE, tol=1e-8) 
{ # internal function
  if(family == "SN")  dp <- msn.cp2dp(cp, silent)
  else if(family == "ESN") stop("cp2dp for ESN not yet implemented")
  else if(family == "ST") dp <- mst.cp2dp(cp, silent, tol=tol)
  else if(family == "SC") stop("this makes no sense for SC family")
  else warning(gettextf("family = '%s' is not supported", family), domain = NA)
  return(dp)
}


msn.cp2dp <- function(cp, silent=FALSE) {
  beta <- cp[[1]]
  Sigma <- cp[[2]]
  gamma1 <- cp[[3]]
  d <- length(gamma1)
  b <- sqrt(2/pi)  
  max.gamma1 <- 0.5*(4-pi)*(2/(pi-2))^1.5
  if(any(abs(gamma1) >= max.gamma1))  
    {if(silent) return(NULL) else stop("non-admissible CP")}
  R <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^(1/3)
  delta <-  R/(b*sqrt(1+R^2))
  mu.z <- b*delta
  omega <- sqrt(diag(Sigma)/(1-mu.z^2))
  Omega <- Sigma + outer(mu.z*omega, mu.z*omega)
  Omega.bar <- cov2cor(Omega)
  Obar.inv <- pd.solve(Omega.bar, silent=silent)
  if(is.null(Obar.inv))  
    {if(silent) return(NULL) else stop("non-admissible CP")}
  Obar.inv.delta <- as.vector(Obar.inv %*% delta)
  delta.sq <- sum(delta * Obar.inv.delta)
  if(delta.sq >= 1) 
    {if(silent) return(NULL) else stop("non-admissible CP")}
  alpha <- Obar.inv.delta/sqrt(1-delta.sq)
  if(is.vector(beta)) {
    beta <- beta - omega*mu.z
    dp <- list(beta=beta, Omega=Omega, alpha=alpha)
    }
  else {
    beta[1,] <- beta[1,] - omega*mu.z
    dp <- list(beta=beta, Omega=Omega, alpha=alpha)  
    }
  attr(dp, "delta.star") <- sqrt(delta.sq)
  return(dp)
  }

st.dp2cp <- 
function(dp, cp.type="proper", fixed.nu=NULL, jacobian=FALSE, upto=NULL) 
{
  if(any(is.na(dp))) stop("NA's in argument 'dp'")
  if(!(cp.type %in% c("proper", "pseudo"))) stop("invalid cp.type") 
  nu <- if(is.null(fixed.nu)) dp[length(dp)] else fixed.nu
  if(is.null(upto)) upto <- 4L
  if((round(upto) != upto)||(upto < 1)) stop("'upto' must be positive integer")
  if(nu <= upto && (cp.type =="proper")) return(NULL)
  p <- length(dp) - 2 - is.null(fixed.nu)
  beta1 <- if(p>1) dp[2:p] else  NULL   
  dp <- c(dp[1], dp[p+1], dp[p+2], nu)
  a <- if(cp.type == "proper") rep(0,upto) else (1:upto) 
  omega <- dp[2]
  alpha <- dp[3]
  delta <- delta.etc(alpha)
  mu.z <- function(delta, nu) delta*b(nu)
  mu <- dp[1] + dp[2]* mu.z(delta, nu+a[1])
  cp <- c(mu, beta1)
  if(upto > 1) {
    kappa2 <- function(delta,nu) nu/(nu-2) - mu.z(delta,nu)^2
    sigma <- omega * sqrt(kappa2(delta, nu+a[2]))
    cp <- c(cp, sigma)
    }
  if(upto > 2) {
    g1 <- st.gamma1(delta, nu+a[3])
    cp <- c(cp, g1)
    }
  if(upto > 3) { 
    g2 <- st.gamma2(delta, nu+a[4])
    cp <- c(cp, g2)}
  rv.comp <- c(rep(TRUE,upto-1), rep(FALSE, 4-upto))
  param.type <- switch(cp.type, proper="CP", pseudo="pseudo-CP")
  names(cp) <- param.names(param.type, "ST", p, x.names=names(beta1), rv.comp)
  if(!is.null(fixed.nu) && upto==4) cp <- cp[-length(cp)]
  if(jacobian && (nu+a[3] > 3)) {
    u <- function(nu) 0.5*(1/nu + digamma((nu-1)/2) - digamma(nu/2)) 
    Ddelta <- 1/(1+alpha^2)^1.5
    Dkappa2.nu <- function(delta,nu) 
      (-2)*(1/(nu-2)^2 + mu.z(delta,nu)^2 * u(nu))
    Dg1.delta <- function(delta,nu) { # derivative di gamma1 wrt delta
      k2 <- kappa2(delta,nu)
      tmp <- nu/(nu-2)-delta^2*(nu-2*b(nu)^2*(nu-2))    
      (3*b(nu) *nu *tmp)/(k2^2.5 * (nu-2)*(nu-3))
      }
    Dg1.nu <-  function(delta,nu) {# derivative di gamma1 wrt nu
      k1 <- mu.z(delta,nu)
      k2 <- kappa2(delta,nu)
      Dk2.nu <- Dkappa2.nu(delta,nu)
      (g1*u(nu)
       + k1/k2^1.5*(-3*(3-delta^2)/(nu-3)^2 + 6/(nu-2)^2 + 4*k1^2*u(nu))
       -3*g1*Dk2.nu/(2*k2))
       }
    Dg2.delta <- function(delta,nu) {# derivative di gamma2 wrt delta
      k1 <- mu.z(delta, nu)
      k2 <- kappa2(delta,nu)
      4*b(nu)^2*delta/k2 * (g2 + 3 -(2*(3-2*delta^2)*nu/(nu-3)
                 -3*nu/(nu-2)+3*k1^2)/k2)
      }
    Dg2.nu <- function (delta, nu) {# derivative di gamma2 wrt nu
      k1 <- mu.z(delta, nu)
      k2 <- kappa2(delta,nu)
      b. <- b(nu)
      u. <- u(nu)
      k4 <- (3 * nu^2/((nu - 2) * (nu - 4))
              -6*(delta*b.)^2 * nu*(nu-1)/((nu-2)*(nu-3))
              + delta^4 * b.^2* (4*nu/(nu-3)-3*b.^2))
      Dk4.nu <- (-6*nu*(3*nu-8)/((nu-2)*(nu-4))^2
               -4*k1^2*(3-delta^2)*((2*u.*nu+1)*(nu-3)-nu)/(nu-3)^2
               +6*k1^2*((2*u(nu)*nu+1)*(nu-2)-nu)/(nu-2)^2
               -12*k1^4*u.)
      Dk2.nu <- Dkappa2.nu(delta,nu)
      Dk4.nu/k2^2 - 2*k4*Dk2.nu/k2^3
      }
    Dcp.dp <- if(is.null(fixed.nu)) diag(1, p+3) else  diag(1, p+2) 
    Dcp.dp[1, p+1] <- mu.z(delta, nu+a[1])
    Dcp.dp[1, p+2] <- omega * Ddelta * b(nu+a[1])
    sigma.z <- sqrt(kappa2(delta, nu+a[2]))
    Dcp.dp[p+1,p+1] <- sigma.z
    Dcp.dp[p+1,p+2] <- -omega *delta *b(nu+a[2])^2 *Ddelta/sigma.z
    Dcp.dp[p+2,p+2] <- Dg1.delta(delta, nu+a[3]) * Ddelta
    if(is.null(fixed.nu) && (nu+a[4] > 4)) {
      Dcp.dp[1, p+3] <- omega * mu.z(delta, nu+a[1]) * u(nu+a[1])
      Dcp.dp[p+1,p+3] <- omega * Dkappa2.nu(delta, nu+a[2])/(2 * sigma.z)
      Dcp.dp[p+2,p+3] <- Dg1.nu(delta, nu+a[3])
      Dcp.dp[p+3,p+2] <- Dg2.delta(delta, nu+a[4]) * Ddelta 
      Dcp.dp[p+3,p+3] <- Dg2.nu(delta, nu+a[4])
      }
    attr(cp, "jacobian") <- Dcp.dp
    }
  return(cp)
}

# b <- function (nu)  ifelse(nu>1, ifelse(nu < 1e8, 
#        sqrt(nu/pi)*exp(lgamma((nu-1)/2)-lgamma(nu/2)), sqrt(2/pi)), NA)

b <- function(nu){ 
   out <- rep(NA, length(nu))
   big.nu <- 1e4
   big <- (nu > big.nu)
   ok  <- ((nu > 1) & (!big) & (!is.na(nu)))  
   out[big] <- sqrt(2/pi) * (1 + 0.75/nu[big] + 0.78125/nu[big]^2)
   out[ok] <-  sqrt(nu[ok]/pi) * exp(lgamma((nu[ok]-1)/2) - lgamma(nu[ok]/2))
   out}
#
st.gamma1 <- function(delta, nu)
{# this function is vectorized for delta, works with a single value of nu
  if(nu > 1e6) {  
    mu <- delta*sqrt(2/pi)
    return(0.5*(4-pi)*mu^3/(1-mu^2)^1.5)
    }
  if(nu > 3) {
    mu <- delta*b(nu)
    k2 <-  nu/(nu-2)- mu^2
    k3 <- mu * (nu * (3 - delta^2)/(nu-3) -3 * nu/(nu - 2) + 2 * mu^2)
    gamma1 <- k3/sqrt(k2)^3 } 
  else 
    gamma1<- Inf*sign(delta)
  gamma1
}
#     
st.gamma2 <- function(delta, nu) 
{# this function is vectorized for delta, works a single value of nu
 #
  if(nu > 1e6) {  
    mu <- delta*sqrt(2/pi)
    return(2*(pi-3)*mu^4/(1-mu^2)^2)
    }
  if(nu > 4) {
    mu <- delta*b(nu)
    k2 <- nu/(nu-2)- mu^2
    k4 <- (3 * nu^2/((nu - 2) * (nu - 4))
           - 4 * mu^2 * nu * (3 - delta^2)/(nu - 3)
           + 6 * mu^2 * nu/(nu - 2) -3*mu^4)                  
    gamma2 <- k4/k2^2 - 3 }   
  else  
    gamma2 <- Inf
  gamma2
  }
#
st.cp2dp <- function(cp, silent=FALSE, tol=1e-8, trace=FALSE) 
{
  fn0 <- function(log.nu, g1) st.gamma1(1, exp(log.nu)) - g1
  if(any(is.na(cp))) stop("NA's in argument 'cp'")
  p <- length(cp)-3
  x.names <- if(p>1) names(cp[2:p]) else NULL
  gamma1 <- cp[p+2]
  abs.g1 <- abs(gamma1)
  gamma2 <- cp[p+3]
  if(abs.g1 <=  0.5*(4-pi)*(2/(pi-2))^1.5) {
    # book: (2.29)+(3.20)
    feasible <- (gamma2 > 2*(pi-3)*(2*abs.g1/(4-pi))^(4/3))
    excess <- max(0, 2*(pi-3)*(2*abs.g1/(4-pi))^(4/3) - gamma2)
    } 
  else {
    if(abs.g1 >= 4) { feasible <- FALSE; excess <- Inf } else {
      r0 <- uniroot(fn0, interval=c(log(4), 1000), tol=tol, g1=abs.g1)
      nu0 <- exp(r0$root) 
      feasible <- (gamma2 >= st.gamma2(1,nu0))
      excess <- max(0, st.gamma2(1,nu0) - gamma2)
      }
    }
  if(!feasible) {
    if(silent) {
      out <- NA
      attr(out, "excess") <- excess
      return(out)} 
    else stop("CP outside feasible region")}
  delta <- 0.75 * sign(gamma1)
  old <- c(delta, Inf)
  step <- Inf
  fn1 <- function(delta, g1, nu) st.gamma1(delta, nu) - g1
  fn2 <- function(log.nu, g2, delta) st.gamma2(delta, exp(log.nu)) - g2
  out <- NULL
  while(step > tol){
    fn21 <- fn2(log(4 + sqrt(.Machine$double.eps)), gamma2, delta)
    fn22 <- fn2(log(1e4), gamma2, delta)
    if(any(is.na(c(fn21, fn22)))) stop("parameter inversion failed")  
    if(fn21 * fn22 > 0) {
      out <- NA
      attr(out, "excess") <- fn21*fn22
      break}
    r2 <- uniroot(fn2, interval=c(log(4 + sqrt(.Machine$double.eps)), 100), 
           tol=tol, g2=gamma2, delta=delta)
    nu <- exp(r2$root)
    if(fn1(-1, gamma1, nu) * fn1(1, gamma1, nu)> 0) {
      out <- NA
      attr(out, "excess") <- fn1(-1, gamma1, nu) * fn1(1, gamma1, nu)
      break}
    r1 <- uniroot(fn1, interval=c(-1,1), tol=tol, g1=gamma1, nu=nu)
    delta <- r1$root
    new <- c(delta, nu)
    step <- abs(old-new)[1] + abs(log(old[2])- log(new[2]))
    if(trace) 
      cat("delta, nu, log(step):", format(c(delta, nu, log(step))),"\n")
    old <- new
    }
  if(anyNA(out)) return(out)
  mu.z <- delta * b(nu)
  omega <- cp[p+1]/sqrt(nu/(nu-2) - mu.z^2)
  if(omega < 0) {
    if(silent) {
      out <- NA
      attr(out, "excess") <- abs(omega)
      return(out)} 
    else stop("CP outside feasible region")}
  alpha <- delta/sqrt(1-delta^2)
  dp <- c(cp[1] - omega*mu.z, if(p>1) cp[2:p] else NULL, omega, alpha, nu)
  names(dp) <- param.names("DP", "ST", p, x.names=x.names)
  return(dp)
}

mst.cp2dp <- function(cp, silent=FALSE, tol=1e-8, trace=FALSE) 
{
  mu <- drop(cp[[1]])
  Sigma <- cp[[2]]
  gamma1 <- cp[[3]]
  gamma2M <- cp[[4]]
  d <- length(gamma1)
  # fn1 <- function(delta, g1, nu) st.gamma1(delta, nu) - g1
  # fn2 <- function(log.nu, g2, delta.sq, d)
  #                mst.gamma2M(delta.sq, exp(log.nu), d) - g2
  if(any(abs(gamma1) >= 4)) 
    {if(silent) return(NULL) else stop("cp$gamma1 not admissible")}
  dp.marg <- matrix(NA, d, 4)
  for(j in 1:d) {  
     dp <- st.cp2dp(c(0,1,gamma1[j], gamma2M), silent=silent)
     if(is.null(dp)) 
       {if(silent) return(NULL) else stop("no CP could be found")}
     dp.marg[j,] <- dp
  }
  if(trace) {cat("starting dp:\n"); print(dp.marg)}
  fn <- function(par, Sigma, gamma1, gamma2M, trace=FALSE){
    if(trace)  cat("[mst.cp2dp[fn]] par:", format(par), "\n")
    nu <- exp(par[1])+4
    delta <- par[-1]/sqrt(1+par[-1]^2)
    d <- length(delta)
    mu.z <- delta*b(nu)
    omega <- sqrt(diag(Sigma)/(nu/(nu-2)-mu.z^2))
    Omega.bar <- (diag(1/omega, d, d) %*% Sigma %*% diag(1/omega, d, d)
                   + outer(mu.z, mu.z)) * (nu-2)/nu
    Obar.inv <- pd.solve(force.symmetry(Omega.bar))
    delta.sq <- sum(delta * as.vector(Obar.inv %*% delta))
    if(delta.sq >= 1) return(delta.sq*10^10)
    L1 <- sum((st.gamma1(delta, nu) - gamma1)^2)
    L2 <- (mst.gamma2M(delta.sq, nu, d) - gamma2M)^2
    # if(trace){  ecat(c(nu,delta,L1,L2))} # ; readline("<cr>")}
    L1 + L2
    }
  nu <- min(dp.marg[,4])
  par <- c(log(nu-4), dp.marg[,3])
  if(trace) cat("[mst.cp2dp] par:", format(par), "\n")
  opt <- nlminb(par, fn, Sigma=Sigma, gamma1=gamma1, gamma2M=gamma2M,
                trace=trace)
  if(trace) cat("[mst.cp2dp]\nopt$convergence:", opt$convergence, 
                "\nopt$message", opt$message, "\n")
  if(opt$convergence != 0) 
    { if(silent) return(NULL) else stop ("no CP could be found") }
  par <- opt$par
  nu <- exp(par[1])+4
  delta <- par[-1]/sqrt(1+par[-1]^2)
  if(trace) {
    cat("[mst.cp2dp]min opt$fn:", format(opt$obj),"\n")
    print(c(nu,delta))
    }
  mu.z <- delta*b(nu)
  omega<- sqrt(diag(Sigma)/(nu/(nu-2)-mu.z^2))
  Omega.bar <- (diag(1/omega, d, d) %*% Sigma %*% diag(1/omega, d, d)
                   + outer(mu.z,mu.z)) * (nu-2)/nu
  Obar.inv <- pd.solve(Omega.bar)
  delta.sq <- sum(delta * as.vector(Obar.inv %*% delta))
  alpha <- as.vector(Obar.inv %*% delta)/sqrt(1-delta.sq)
  if(is.matrix(mu)) {
     xi <- mu
     xi[1,] <- mu[1,] - omega*mu.z }
  else xi <- mu - omega*mu.z
  Omega <- diag(omega) %*% Omega.bar %*% diag(omega)
  return(list(xi=xi, Omega=Omega, alpha=alpha, nu=nu))
}
 

affineTransSECdistr <- function(object, a, A, name, compNames, drop=TRUE)
{# object is of class SECdistrMv
 # computes distribution of affine transformation of SEC variable T=a+t(A)Y
  if(class(object) != "SECdistrMv") stop("wrong object class")
  dp <- slot(object, "dp")
  alpha <- dp$alpha
  d <- length(alpha)
  if(!is.matrix(A) || nrow(A) != d) stop("A is not a matrix or wrong nrow(A)")
  h <- ncol(A)
  if(length(a) != h) stop("size mismatch of arguments 'a' and 'A'")
  if(missing(name)) name<- paste(deparse(substitute(a)), " + t(",  
    deparse(substitute(A)), ") %*% (", deparse(substitute(object)),")", sep="")
  else name <- as.character(name)[1]
  compNames <- if(missing(compNames)) 
    as.vector(outer("V",as.character(1:h),paste,sep=""))
    else as.character(as.vector(compNames)[1:h])
  family <- object@family
  xi.X  <- as.vector(a + t(A) %*% matrix(dp$xi, ncol=1))
  Omega <- dp$Omega
  omega <- sqrt(diag(Omega))
  Omega.X <- as.matrix(t(A) %*% Omega %*% A) 
  invOmega.X <- pd.solve(Omega.X, silent=TRUE)
  if (is.null(invOmega.X)) stop("not full-rank transformation") 
  omega.X <- sqrt(diag(Omega.X))
  omega.delta <- omega * delta.etc(alpha, Omega)$delta
  m <- as.vector(invOmega.X %*% t(A) %*% matrix(omega.delta, ncol=1))
  u <- sum(omega.delta * as.vector(A %*% matrix(m, ncol=1)))
  alpha.X <- (omega.X * m)/sqrt(1 - u)
  dp.X <- list(xi=xi.X, Omega=Omega.X, alpha=alpha.X)
  if(family == "ESN") dp.X$tau <- dp$tau
  if(family == "ST") dp.X$nu <- dp$nu
  if(h==1 & drop) {
     dp1 <- unlist(dp.X)
     dp1[2] <- sqrt(dp1[2])
     names(dp1) <- names(dp.X) 
     names(dp1)[2] <- tolower(names(dp)[2])
     # new.obj <- new("SECdistrUv", dp=dp1, family=family, name=name) #??
     new.obj <- makeSECdistr(dp=dp1, family=family, name=name)
     } else 
  new.obj <- makeSECdistr(dp.X, family, name, compNames) 
  # new.obj <- new("SECdistrMv", dp.X, family, name, compNames) #??
  return(new.obj)
}
  
                       
marginalSECdistr <- function(object, comp, name, drop=TRUE)   
{# marginals of SECdistrMv obj; 2nd version, computing marginal delta's
  family <- slot(object,"family")
  if(missing(name)) {
     basename <- if(object@name != "") object@name 
                 else deparse(substitute(object))
     name<- paste(basename, ".components=(",
                paste(as.character(comp),collapse=","), ")", sep="")
     }
  else name <- as.character(name)[1]
  dp <- slot(object,"dp")
  xi    <- dp$xi
  Omega <- dp$Omega
  alpha <- dp$alpha
  compNames <- slot(object,"compNames")
  d <- length(alpha)
  comp <- as.integer(comp)
  Omega11 <- Omega[comp,comp,drop=FALSE]
  if(length(comp) < d){
    if(any(comp>d | comp<1)) stop("comp makes no sense")
    delta_etc <- delta.etc(alpha, Omega)
    delta1 <- delta_etc$delta[comp]
    R11 <- delta_etc$Omega.cor[comp, comp, drop=FALSE]
    iR11.delta1 <- as.vector(pd.solve(R11, silent=TRUE) %*% delta1)
    diRd <- sum(delta1*iR11.delta1)
    alpha1_2 <- if(diRd < 1) iR11.delta1/sqrt(1 - diRd) else sign(delta1)*Inf
    dp0 <- list(xi=xi[comp], Omega=Omega11, alpha=alpha1_2)
  }
  else {
    if(any(sort(comp) != (1:d))) stop("comp makes no sense")
    dp0 <- list(xi=xi[comp], Omega=Omega11, alpha=alpha[comp])
  }
  if(family=="ESN") dp0$tau <- dp$tau
  if(family=="ST") dp0$nu <- dp$nu
  new.obj <- new("SECdistrMv", dp=dp0, family=family, name=name, 
                 compNames=compNames[comp])
  if(length(comp)==1 & drop) 
    {# new.obj <- as(new.obj, "SECdistrUv") # non va..
     dp <- unlist(dp0)
     names(dp) <- names(dp0)
     dp[2] <- sqrt(dp[2])
     names(dp)[2] <- "omega"
     new.obj <- new("SECdistrUv", dp=dp, family=family, name=compNames[comp])
     }
  new.obj
}
                     
conditionalSECdistr <- 
function(object, fixed.comp, fixed.values, name, drop=TRUE)
{ # conditional distribution of SN/ESN object 
  family <- slot(object,"family")
  if(!(family %in% c("SN", "ESN"))) stop("family must be either SN or ESN")
  dp <- slot(object,"dp")
  xi    <- dp$xi
  Omega <- dp$Omega
  alpha <- dp$alpha
  tau   <- if(family=="SN") 0 else dp$tau
  d <- length(alpha)
  fix <- fixed.comp
  h <- length(fix)
  if(any(fix != round(fix)) | !all(fix %in% 1:d) | h == d) 
    stop("fixed.comp makes no sense")
  if(length(fixed.values) != h) 
    stop("length(fixed.comp) != lenght(fixed.values)")
  compNames <- slot(object,"compNames")
  if(missing(name)) {
     basename <- if(object@name != "") object@name 
                 else deparse(substitute(object))
     name<- paste(basename,"|(",
                paste(compNames[fix],collapse=","), ")=(", 
                paste(format(fixed.values),collapse=","), ")",
                sep="")
     }
  else name <- as.character(name)[1]
  # free.fix <- setdiff(1:d, fix)
  omega <- sqrt(diag(Omega))
  omega1 <- omega[fix]
  omega2 <- omega[-fix]
  R   <- cov2cor(Omega)
  R11 <- R[fix,fix, drop=FALSE]
  R12 <- R[fix,-fix, drop=FALSE]
  R21 <- R[-fix,fix, drop=FALSE]
  R22 <- R[-fix,-fix, drop=FALSE]
  alpha1 <- matrix(alpha[fix], ncol=1)
  alpha2 <- matrix(alpha[-fix], ncol=1)
  iR11  <- pd.solve(R11)
  R22.1 <- R22 - R21 %*% iR11 %*% R12
  a.sum <- as.vector(t(alpha2) %*% R22.1 %*% alpha2)
  alpha1_2 <- as.vector(alpha1 + iR11 %*% R12 %*% alpha2)/sqrt(1+a.sum)
  tau2.1 <- (tau * sqrt(1 + sum(alpha1_2 * as.vector(iR11 %*% alpha1_2)))
             + sum(alpha1_2 * (fixed.values-xi[fix])/omega1))
  O11 <- Omega[fix,fix, drop=FALSE]
  O12 <- Omega[fix,-fix, drop=FALSE]
  O21 <- Omega[-fix,fix, drop=FALSE]
  O22 <- Omega[-fix,-fix, drop=FALSE]
  iO11<- (1/omega1) * iR11 * rep(1/omega1, each=h)  # solve(O11)
  reg <- O21 %*% iO11
  xi2.1 <- as.vector(xi[-fix]+ reg %*% (fixed.values - xi[fix]))
  O22.1 <- O22 - reg %*% O12
  omega22.1 <- sqrt(diag(O22.1))
  alpha2.1 <- as.vector((omega22.1/omega2)*alpha2)
  dp2.1 <- list(xi=xi2.1, Omega=O22.1, alpha=alpha2.1, tau=tau2.1)
  obj <- if((d-h)==1 & drop) {
    dp2.1 <- unlist(dp2.1)
    dp2.1[2] <- sqrt(dp2.1[2])
    names(dp2.1) <- c("xi","omega","alpha","tau")
    new("SECdistrUv", dp=dp2.1, family="ESN", name=name)
    } else new("SECdistrMv", dp=dp2.1, family="ESN", name=name, 
               compNames=compNames[-fix])   
  return(obj)
}


delta.etc <- function(alpha, Omega=NULL) 
{ 
  inf <- which(abs(alpha) == Inf)
  if(is.null(Omega)){ # case d=1
    delta <- alpha/sqrt(1+alpha^2)
    delta[inf] <- sign(alpha[inf])
    return(delta)
    }
  else { # d>1
    if(any(dim(Omega) != rep(length(alpha),2))) stop("dimension mismatch")
    Ocor <- cov2cor(Omega)
    if(length(inf) == 0) { # d>1, standard case
      Ocor.alpha <- as.vector(Ocor %*% alpha)
      alpha.sq <- sum(alpha * Ocor.alpha)
      delta <- Ocor.alpha/sqrt(1+alpha.sq)
      alpha. <- sqrt(alpha.sq)
      delta. <- sqrt(alpha.sq/(1+alpha.sq))
      }
     else { # d>1, case with some abs(alpha)=Inf
       if(length(inf) > 1) 
         warning("Several abs(alpha)==Inf, I handle them as 'equal-rate Inf'") 
       k <- rep(0,length(alpha))
       k[inf] <- sign(alpha[inf])
       Ocor.k <- as.vector(Ocor %*% k) 
       delta <- Ocor.k/sqrt(sum(k * Ocor.k))
       delta. <- 1
       alpha. <- Inf
       }
  return(
    list(delta=delta, alpha.star=alpha., delta.star=delta., Omega.cor=Ocor))
  }
}

selm <- function (formula, family="SN", data, weights, subset, na.action, 
    start=NULL, fixed.param=list(), method="MLE",  penalty=NULL, offset, 
    model=TRUE, x = FALSE, y = FALSE, ...) 
{
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    formula <- as.formula(formula)
    if (length(formula) < 3)  stop("formula must be a two-sided formula")
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    method <- toupper(method)
    if(!(method %in% c("MLE", "MPLE"))) {
      warning(gettextf("method = '%s' is not supported, replaced by 'MLE'", 
         method), domain = NA)
      method <- "MLE"}
    penalty.name <- if(method == "MPLE") {
      if(is.null(penalty)) "Qpenalty" else penalty }
      else   NULL  
    contr <- list(penalty=penalty.name, trace=FALSE,  info.type="observed", 
                  opt.method="nlminb", opt.control=list())
    control <- list(...)
    contr[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% names(contr)])) warning(
       "unknown names in control: ", paste(noNms, collapse = ", "))
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    if(is.null(w))  w <- rep(1, NROW(y))      
    if(any(w != round(w)) | all(w == 0))
      stop("weights must be non-negative integers (=frequencies), not all 0")  
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
      if (length(offset) == 1) 
        offset <- rep(offset, NROW(y))
      else if (length(offset) != NROW(y)) 
        stop(gettextf(
         "number of offsets is %d, should equal %d (number of observations)", 
         length(offset), NROW(y)), domain = NA)
      }          
    if(length(fixed.param) > 0) {
      if(!all(names(fixed.param)  %in%  c("nu", "alpha")))
        stop("Not admissible component of 'fixed.param'")
      if(!is.null(fixed.param$alpha)) { 
        if(fixed.param$alpha != 0) stop("'alpha' can only be fixed at 0")
        if(method == "MPLE") stop('method MPLE not allowed when alpha=0')      
        }
      }  
    if (is.empty.model(mt)) stop("empty model") else
    {
      x <- model.matrix(mt, mf, contrasts)
      xt <- pd.solve(t(x) %*% (w*x), silent=TRUE)
      if(is.null(xt)) stop("design matrix appears to be of non-full rank")
      z <- selm.fit(x, y, family=family, start, w=w, fixed.param=fixed.param, 
             offset=offset, selm.control=contr)
    }
    class(z) <- c(if (is.matrix(y)) "mselm", "selm")
    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    input <- list()
    if (model) input$model <- mf
    if (ret.x) input$x <- x
    if (ret.y) input$y <- y
    # input$weights <- as.vector(model.weights(mf))
    # input$offset <- as.vector(model.offset(mf))
    # cl.obj <- if(is.matrix(y)) "mselm" else "selm"
    obj <- new(class(z), call=cl, family=toupper(family), logL=z$logL, 
               method=c(method, contr$penalty),  param=z$param,
               param.var=z$param.var, size=z$size,  
               residuals.dp=z$resid.dp, fitted.values.dp=z$fitted.dp,
               control=control, input=input, opt.method=z$opt.method)
    return(obj)
}
#
#selm.control <- function(method="MLE", info.type="observed",  
#   trace=FALSE, algorithm="nlminb", opt.control=list()) 
#{     
#  if(algorithm !="nlminb") stop("only algorithm='nlminb' handled so far")
#  if(info.type !="observed") stop("only info.type='observed' handled so far")
#  list(method=method, info.type=info.type,  trace=trace, 
#    algorithm=algorithm, opt.control=opt.control)
#}


#------------------------------------------------------
selm.fit <- function (x, y, family="SN", start=NULL, w, fixed.param=list(), 
                 offset = NULL, selm.control) 
{
    if (!(toupper(family) %in% c("SN", "ST", "SC")))
        stop(gettextf("I do not know family '%s'", family), domain = NA)
    family <- toupper(family)    
    if (is.null(n <- nrow(x))) stop("'x' must be a matrix")
    if (n == 0L) stop("0 (non-NA) cases")
    if(NROW(y) != n) stop("'x' and 'y' have non-compatible dimensions")
    p <- ncol(x)
    if ((p == 0L) || !(all(data.matrix(x)[,1] == 1))) 
      stop("first column of model matrix is not all 1's")
    y <- drop(y)
    d <- NCOL(y)
    if(d>1 && is.null(colnames(y))) colnames(y) <- paste("V", 1:d, sep="") 
    if(is.null(colnames(x))) colnames(x) <- paste("x", 0L:(p-1), sep=".")
    if (!is.null(offset))  y <- (y - offset)
    if (NROW(y) != n)  stop("incompatible dimensions")
    if (missing(w) || is.null(w)) w <- rep(1, n)
    nw <- sum(w)
    n.obs <- NROW(y)
    contr <- list(method="MLE", penalty=NULL, trace=FALSE,  
                 info.type="observed", opt.method="nlminb", opt.control=list())
    control <- selm.control
    contr[(namc <- names(control))] <- control   
    symmetr <- FALSE   
    if(length(fixed.param) > 0) {
      if(!all(names(fixed.param)  %in%  c("nu", "alpha")))
        stop("Not admissible component of 'fixed.param'")
      if(!is.null(fixed.param$alpha)) {      
        if( fixed.param$alpha != 0 ) stop("'alpha' can only be fixed at 0") 
        else symmetr <- TRUE }
      }    
    zero.weights <- any(w == 0)
    if(zero.weights) {
      save.r <- y
      save.f <- y
      save.w <- w
      ok <- (w != 0)
      nok <- !ok
      w <- w[ok]
      x0 <- x[!ok, , drop = FALSE]
      x <- x[ok, , drop = FALSE]
      n <- nrow(x)
      y0 <- if (d > 1L) y[!ok, , drop = FALSE] else y[!ok]
      y <- if (d > 1L) y[ok, , drop = FALSE] else y[ok]
      }
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    info.type <- contr$info.type # so far, only "observed"
    yInfo <- if(contr$info.type == "observed") y else NULL
    penalty <- contr$penalty  # either NULL or a char string 
    penalty.fn <- if(is.null(penalty)) NULL else get(penalty, inherits=TRUE) 
    trace <- contr$trace
    if(d == 1) {
      y <- as.vector(y) 
      if(family == "SN") {
        npar <- p + 2 - as.numeric(symmetr)
        if(symmetr) { # SN with alpha=0 is Gaussian case
          ls <- lm.wfit(x, y, w) # note: offset already subtracted if any
          res <- residuals(ls)
          s2 <- sum(w*res^2)/nw
          param <- c(coef(ls), sqrt(s2))
          j <- rbind(cbind(t(x) %*% (w*x)/s2, 0), c(rep(0,p), 2*nw/s2))
          j.inv <- solve(j)
          se <- sqrt(diag(j.inv))
          info <- list(dp=param, cp=param, info.dp=j, info.cp=j, 
                   asyvar.dp=j.inv, asyvar.cp=j.inv, se.dp=se, se.cp=se,
                   aux=NULL)
          logL <- (-0.5*nw)*(log(2*pi*s2) +1)
          fit <- list(cp=param, dp=param, dp.complete=c(param,0), 
                      opt.method=list(ls$qr), logL=logL)  
          boundary <- FALSE
          fit$opt.method <- list(method="least_squares", called.by= "lm.wfit")
          mu0 <- 0
          }
        else { # proper SN case 
        cp <- if(is.null(start)) NULL else dp2cpUv(start, "SN")
        fit <- sn.mple(x, y, cp, w, penalty, trace, contr$opt.method, 
                 contr$control)
        fit$dp <- cp2dpUv(cp=fit$cp, family="SN")
        boundary <- fit$boundary
        mu0 <- fit$cp[1] - fit$dp[1]
        info <- if(boundary) NULL else 
          sn.infoUv(dp=fit$dp, x=x, y=yInfo, w=w, penalty=penalty)
        }}
      if(family == "ST") {
        fixed.nu <- fixed.param$nu  
        npar <- p + 2 + as.numeric(is.null(fixed.nu)) - as.numeric(symmetr)
        fit <- st.mple(x, y, dp=start, w, fixed.nu, symmetr, penalty, trace,
           contr$opt.method, contr$control)
        dp <- fit$dp
        cp <- st.dp2cp(dp, cp.type="proper", fixed.nu=fixed.nu, 
                 upto=4-as.numeric(!is.null(fixed.nu)))
        p_cp<- st.dp2cp(dp, cp.type="pseudo", fixed.nu=fixed.nu, jacobian=TRUE)
        fit$cp <- cp[1:npar]     
        fit$p_cp <- p_cp[1:npar]
        Dpseudocp.dp <- attr(p_cp, "jacobian")[1:npar, 1:npar]
        attr(p_cp, "jacobian") <- NULL
        boundary <- fit$boundary
        nu <- if(is.null(fixed.nu)) dp[npar] else fixed.nu
        mu0 <- if(nu <= 1) NA else 
          st.dp2cp(dp, fixed.nu=fixed.nu, upto=1)[1] - dp[1] 
        info <- if(boundary)  NULL  else 
          st.infoUv(dp=fit$dp, NULL, x, yInfo, w, fixed.nu, symmetr, penalty) 
        }
      if(family == "SC") {
        npar <- p + 2  - as.numeric(symmetr)
        fit <- st.mple(x, y, dp=start, w, fixed.nu=1, symmetr, penalty, trace,
                       contr$opt.method, contr$control)
        fit$cp <- NULL
        p_cp0 <- st.dp2cp(fit$dp, cp.type="pseudo", fixed.nu=1, jacobian=TRUE) 
        fit$p_cp <- p_cp0[1:npar]
        Dpseudocp.dp <- attr(p_cp0, "jacobian")[1:npar, 1:npar]
        attr(p_cp0, "jacobian") <- NULL
        boundary <- fit$boundary
        mu0 <- NA
        info <- if(boundary) NULL  else
          st.infoUv(dp=fit$dp, x=x, y=yInfo, w=w, fixed.nu=1, symmetr=symmetr) 
        }
      if(!boundary && family %in% c("ST","SC"))  info$asyvar.p_cp <- 
          Dpseudocp.dp %*% info$asyvar.dp %*% t(Dpseudocp.dp)
      beta.dp <- fit$dp[1:p]
      dp <- fit$dp
      cp <- fit$cp
      }
    else { # d>1
      npar0 <- p*d + d*(d+1)/2
      if(family == "SN") {
        if(symmetr) { # SN with alpha=0 is Gaussian case
          npar <-  npar0
          ls <- lm.wfit(x, y, w) # note: offset already subtracted if any
          beta <- coef(ls)
          res <- residuals(ls)
          s2 <- t(res) %*% (w*res)/nw
          dp <- dp. <- list(beta=beta, Omega=s2)
          dp.$alpha <- rep(0,d)
          param <- c(beta, vech(s2))
          conc <- solve(s2)
          betaBlock <- conc %x% (t(x) %*% (w*x))
          D <- duplicationMatrix(d)
          varBlock <- (n/2) * t(D) %*% (conc %x% conc) %*% D
          m0 <- matrix(0, p*d, d*(d+1)/2)
          j <- rbind(cbind(betaBlock, m0), cbind(t(m0), varBlock)) 
          # use (10) in section 15.8 of Magnus & Neudecker (1988/1999, p.321)
          j.inv <- rbind(cbind(solve(betaBlock), m0), 
                         cbind(t(m0), solve(varBlock))) 
          diags.dp <- sqrt(diag(j.inv))
          se.beta <- matrix(diags.dp[1:(p*d)], p, d)
          se.diagOmega <- diags.dp[p*d + d*(d+1)/2 +1 -rev(cumsum(1:d))]
          se <- list(beta=se.beta, diagOmega=se.diagOmega)
          info <- list(dp=param, cp=param, info.dp=j, info.cp=j, 
                   asyvar.dp=j.inv, asyvar.cp=j.inv, se.dp=se, se.cp=se,           
                   aux=NULL)
          logL <- (-0.5*nw)*(determinant(2*pi*s2, logarithm=TRUE)$modulus + d)
          # see (6.2.7) of Mardia, Kent & Bibby (1979)
          fit <- list(dp=dp, cp=dp, dp.complete=dp., logL=logL)  
          fit$opt.method <- list(method="lm.wfit")
          boundary <- FALSE
          mu0 <- rep(0, d)
          }
        else { # proper SN case 
        npar <-  npar0 + d
        if(is.null(penalty)) { # MLE
          fit <- msn.mle(x, y, start, w, trace=trace, 
                 opt.method=contr$opt.method, control=contr$opt.control)
          boundary <- ((1 - fit$aux$delta.star) < .Machine$double.eps^(1/4))
          if(!boundary) info <- 
             sn.infoMv(fit$dp, x=x, y=yInfo, w=w)
          } else { # MPLE
          fit <- msn.mple(x, y, start, w, penalty, trace=trace, 
                   opt.method=contr$opt.method, control=contr$opt.control)
          boundary <- FALSE
          info <- sn.infoMv(fit$dp, x=x, y=y, w=w, penalty=penalty)
          }
        fit$cp <- msn.dp2cp(fit$dp)
        mu0 <- as.vector(fit$cp[[1]][1,] - fit$dp[[1]][1,])
        }}
      if(family == "ST"){
        fixed.nu <- fixed.param$nu 
        npar <- npar0 + d*as.numeric(!symmetr) + as.numeric(is.null(fixed.nu))
        fit <- mst.mple(x, y, start, w, fixed.nu=fixed.nu, symmetr=symmetr,
                  penalty=penalty, trace=trace, opt.method=contr$opt.method, 
                  control=contr$opt.control)
        fit$opt.method$called.by <- "mst.mple"
        boundary <- fit$boundary
        dp <- fit$dp
        nu <- if(is.null(fixed.nu)) dp$nu else fixed.nu
        mu0 <- if(nu <= 1) NA else as.vector(mst.dp2cp(dp, fixed.nu=fixed.nu,
          symmetr=symmetr, upto=1)[[1]][1,] - dp[[1]][1,])
        fit$cp <- mst.dp2cp(dp, cp.type="proper", fixed.nu, symmetr)
        fit$p_cp <- mst.dp2cp(dp, cp.type="pseudo", fixed.nu, symmetr)
        if(!boundary) info <- 
           st.infoMv(dp, x=x, y=yInfo, w, fixed.nu, symmetr, penalty)
        }
      if(family == "SC") {
        npar <- npar0 + d*as.numeric(!symmetr)
        if(is.null(start)) {
          fit.sn <- msn.mle(x, y, NULL, w, control=list(rel.tol=1e-4))
          start <- fit.sn$dp  
          }
        fit <- mst.mple(x, y, start, w,  fixed.nu=1, symmetr=symmetr,
                  penalty=penalty, trace=trace,
                  opt.method=contr$opt.method, control=contr$opt.control)
        fit$opt.method$called.by <- "mst.mple"
        npar <- p*d + d*(d+1)/2 + d*as.numeric(!symmetr)
        boundary <- fit$boundary 
        mu0 <- NA
        fit$cp <- NULL
        fit$p_cp <- mst.dp2cp(fit$dp, "pseudo", fixed.nu=1)   
        if(!boundary)  info <-
          st.infoMv(fit$dp, x=x, y=yInfo, w, fixed.nu=1, symmetr, penalty)
        }
      beta.dp <- fit$dp[[1]]
      }
    param <- list(dp=fit$dp, cp=fit$cp, "pseudo-cp"=fit$p_cp, 
                boundary=boundary, mu0=mu0)
    if(!boundary && !is.null(info)) {
      asyvar.dp <- info$asyvar.dp[1:npar, 1:npar] 
      asyvar.cp <- info$asyvar.cp[1:npar, 1:npar]
      asyvar.p_cp <- info$asyvar.p_cp[1:npar, 1:npar]
      param.var <- list(info.type=info.type, dp=asyvar.dp, cp=asyvar.cp, 
        "pseudo-cp"=asyvar.p_cp) 
      } 
    else  param.var <- list()
    dn <- colnames(x)  
    fv <- drop(x %*% beta.dp)
    if(is.matrix(fv)) colnames(fv) <- colnames(y)
    size <- c(d=d, p=p, n.param=npar, n.obs=n.obs, nw.obs=sum(w)) 
    z <- list(call=match.call(), logL=fit$logL, param=param, 
            param.var=param.var, fitted.dp=fv, resid.dp=y-fv, size=size,
            selm.control=contr, opt.method=fit$opt.method)
    r1 <- y - z$resid.dp 
    z$weights <- w
    if (zero.weights) {
        # coef[is.na(coef)] <- 0
        f0 <- x0 %*% beta.dp
        if (d > 1) {
            save.r[ok, ] <- z$resid.dp
            save.r[nok, ] <- y0 - f0
            save.f[ok, ] <- z$fitted.dp
            save.f[nok, ] <- f0
        }
        else {
            save.r[ok] <- z$resid.dp
            save.r[nok] <- y0 - f0
            save.f[ok] <- z$fitted.dp
            save.f[nok] <- f0
        }
        z$resid.dp <- save.r
        z$fitted.dp <- save.f
        z$weights <- save.w
    }
    if(!is.null(offset)) {
      z$fitted.dp <- z$fitted.dp + offset
      r1 <- r1 + offset
      }
    # z$fitted.dp <- r1
    if(length(fixed.param) > 0)  {
        z$param$fixed <- fixed.param 
        z$param$dp.complete <- fit$dp.complete } 
      else z$param$fixed <- z$param$dp.complete<- list() 
    return(z)
}

#---------------------------------------------------

summary.selm <- function(object, param.type="CP", cov=FALSE, cor=FALSE)
{
  fixed <- slot(object, "param")$fixed
  if(length(fixed$alpha==0)>0 && fixed$alpha==0) {
    param.type <- "DP"
    note <- "param.type=DP has been set because of constraint alpha=0"
    } else note <- ""
  lc.param.type <- tolower(param.type) 
  if(!(lc.param.type %in% c("cp", "op", "dp", "pseudo-cp")))
     stop(gettextf("unknown param.type '%s'", param.type), domain = NA)     
  param.type <- switch(lc.param.type, 
     "dp"="DP", "op"="OP", "cp"="CP", "pseudo-cp"="pseudo-CP")
  family <- slot(object,"family")
  if(param.type=="pseudo-CP" && !(family %in% c("ST", "SC"))) 
    stop("pseudo-CP makes sense only for ST and SC families")
  if (!(family %in% c("SN","ST","SC"))) 
     stop(gettextf("family '%s' not (yet) handled", family), domain = NA)
  param <- slot(object, "param")[[lc.param.type]]
  if(param.type=="CP" && is.null(param)) { 
    if(family %in% c("ST", "SC")) {
      {message("CP does not exist. Consider param.type='DP' or 'pseudo-CP'") 
      return(invisible())}}}
  param.var <- slot(object, "param.var")[[lc.param.type]]
  if(is.null(param.var)) param.var <- diag(NA, length(param))
  se <- sqrt(diag(param.var))
  z <- param/se
  param.table <- cbind(param, se, z, 2*pnorm(-abs(z)))
  dimnames(param.table) <- list(names(param), 
    c("estimate","std.err","z-ratio", "Pr{>|z|}"))
  resid <- residuals(object, lc.param.type)
  aux <- list()
  aux$param.cov <- if(cov) param.var else NULL
  aux$param.cor <- if(cor) cov2cor(param.var) else NULL
  out <- new("summary.selm", call=slot(object,"call"), 
           family = slot(object, "family"), 
           logL = slot(object, "logL"),
           method=slot(object, "method"),
           resid = resid, 
           param.type = param.type,
           param.table = param.table,
           param.fixed = fixed,
           control = slot(object, "control"),
           aux = aux,
           boundary=slot(object, "param")$boundary,
           size=object@size,
           note=note)
   out        
}


residuals.selm <- function(object, param.type="CP", ...){
  param.type <- tolower(param.type) 
  if(!(param.type %in% c("cp", "dp", "pseudo-cp"))) 
     stop("param.type must be either 'CP' or 'DP' or 'pseudo-CP'")
  # param <- slot(object, "param")[[param.type]]
  p <- object@size["p"]
  n <- object@size["n.obs"]
  r <- slot(object, "residuals.dp") 
  dp <- slot(object, "param")$dp
  pseudo.mu0 <- (slot(object, "param")$"pseudo-cp"[1] - dp[1])
  resid <- switch(param.type, 
     'dp' = r, 
     'cp' = r - rep(slot(object,"param")$mu0, n),
     'pseudo-cp' = r - rep(pseudo.mu0, n))
  # resid <- resid/param[p+1] # AA: standardize resid?
  w <- slot(object,"input")$weights
  if(!is.null(w)) attr(resid,"weights") <- w
  return(resid)
  }


fitted.selm <- function(object, param.type="CP", ...) {
  param.type <- tolower(param.type) 
  if(!(param.type %in% c("cp", "dp", "pseudo-cp")))
   stop("param.type must be either 'CP' or 'DP' or 'pseudo-CP'")
  # param <- slot(object, "param")[[param.type]]
  n <- object@size["n.obs"]
  dp <- slot(object, "param")$dp
  fit.dp <- slot(object,"fitted.values.dp")
  pseudo.mu0 <- (slot(object, "param")$"pseudo-cp"[1] - dp[1])
  fitted <- switch(param.type,
    'dp' = fit.dp,
    'cp' = fit.dp + rep(slot(object,"param")$mu0, n),
    'pseudo-cp' = fit.dp + rep(pseudo.mu0, n))
  w <- slot(object, "input")$weights
  if(!is.null(w)) attr(fitted,"weights") <- w
  return(fitted)
  }
  
weights.selm <- function(object, ...) slot(object, "input")$weights

summary.mselm <- function(object, param.type="CP", cov=FALSE, cor=FALSE) 
{
  fixed <- slot(object, "param")$fixed
  if(length(fixed$alpha==0)>0 && fixed$alpha==0) {
    param.type <- "DP"
    note <- "param.type=DP has been set because of constraint alpha=0"
    } else note <- ""
  lc.param.type <- tolower(param.type) 
  if(!(lc.param.type %in% c("dp", "op", "cp", "pseudo-cp")))
     stop(gettextf("unknown param.type '%s'", param.type), domain = NA)
  param.type <- switch(lc.param.type, 
     "dp"="DP", "op"="DP", "cp"="CP", "pseudo-cp"="pseudo-CP")
  # OP not yet implemented, so far re-directed to DP     
  family <- slot(object, "family")
  method <- slot(object, "method")
  if(param.type=="pseudo-CP" & !(family %in% c("ST","SC"))) 
    stop("pseudo-CP makes sense only for ST and SC families")
  p <- object@size["p"]
  d <- object@size["d"]
  npar <- object@size["n.param"]
  param <- object@param[[lc.param.type]]
  if(is.null(param) && family %in% c("ST", "SC")) {
    message("CP does not exist. Consider param.type='DP' or 'pseudo-CP'")
    return(invisible())}
  beta <- param[[1]]
  param.var <- slot(object, "param.var")[[lc.param.type]]
  if(object@param$boundary | is.null(param.var)) 
    param.var <- matrix(NA, npar, npar)
  coef.tables <- list()
  par.names <- param.names(param.type, family, p, x.names=rownames(beta)[-1])
  for(j in 1:d) {
    beta.j <- beta[,j]
    var.j <- param.var[((j-1)*p+1):(j*p), ((j-1)*p+1):(j*p), drop=FALSE]
    se.j <- sqrt(diag(var.j))
    z <- beta.j/se.j
    coef.table <- cbind(beta.j, se.j, z, 2*pnorm(-abs(z)))
    dimnames(coef.table) <- list(par.names[1:p], 
      c("estimate","std.err","z-ratio", "Pr{>|z|}"))
    coef.tables[[j]] <- coef.table
    }
  scatter <- list(matrix=param[[2]], name=names(param)[2])
  resid <- residuals.mselm(object, param.type)
  # resid <- t(t(resid)/sqrt(diag(scatter$matrix))) # for normalized/std resid
  if(is.null(fixed$alpha)) {
    se.slant <- sqrt(diag(param.var)[(p*d+d*(d+1)/2+1):(p*d+d*(d+1)/2+d)])
    slant <- list(param=param[[3]], se=se.slant, name=names(param)[3])} 
    else { if(fixed$alpha == 0) slant <- list() else 
       stop('cannot have fixed alpha at non-zero value, please report')} 
  tail <- if(family== "ST" & is.null(fixed$nu) )
             list(param=param[[length(param)]],     
              se=sqrt(diag(param.var)[npar]),  name=names(param)[length(param)]) 
          else list()
  aux <- list()
  aux$param.cov <- if(cov) param.var else NULL
  aux$param.cor <- if(cor) cov2cor(param.var) else NULL
  out <- new("summary.mselm", call=slot(object,"call"), 
           family = family, 
           logL = slot(object, "logL"),
           method=slot(object, "method"),
           resid = resid,
           param.type=param.type,
           coef.tables = coef.tables,
           param.fixed = fixed,
           scatter = scatter,
           slant = slant,
           tail = tail,
           control = slot(object, "control"),
           aux = aux,
           boundary=slot(object, "param")$boundary,
           size=slot(object, "size"))
   out        
}

residuals.mselm <- function(object, param.type="CP", ...){
  param.type <- tolower(param.type) 
  if(!(param.type %in% c("cp", "dp", "pseudo-cp"))) 
     stop("param.type must be either 'CP' or 'DP' or 'pseudo-CP'")
  # param <- slot(object, "param")[[param.type]]
  # beta <- param[[1]]
  n <- object@size["n.obs"]
  r <- slot(object,"residuals.dp")
  param <- slot(object, "param")
  pseudo.mu0 <- as.vector(param$"pseudo-cp"[[1]][1,] - param$dp[[1]][1, ])
  resid <- switch(param.type, 
    'dp' = r, 
    'cp' = r - outer(rep(1,n), param$mu0),
    'pseudo-cp' = r  - outer(rep(1,n), pseudo.mu0))
  w <- slot(object, "input")$weights
  if(!is.null(w)) attr(resid,"weights") <- w
  return(resid)
  }

fitted.mselm <- function(object, param.type="CP", ...) {
  param.type <- tolower(param.type) 
  if(!(param.type %in% c("cp", "dp", "pseudo-cp"))) 
     stop("param.type must be either 'CP' or 'DP' or 'pseudo-CP'")
  n <- object@size["n.obs"]
  fit.dp <- slot(object, "fitted.values.dp")
  param <- slot(object, "param")
  pseudo.mu0 <- as.vector(param$"pseudo-cp"[[1]][1,] - param$dp[[1]][1, ])
  fitted <- switch(param.type, 
    'dp' = fit.dp, 
    'cp' = fit.dp + outer(rep(1,n), param$mu0),
    'pseudo-cp' = fit.dp + outer(rep(1,n), pseudo.mu0))
  w <- slot(object, "input")$weights
  if(!is.null(w)) attr(fitted,"weights") <- w
  return(fitted)
  }

weights.mselm <- function(object, ...) slot(object, "input")$weights

#------------------------------------------------------------
# 
# sn.info<- function(dp=NULL, cp=NULL, x=NULL, y=NULL, w, penalty=NULL, 
#              type="observed", norm2.tol=1e-6) { 
# if(any(is.list(dp), is.list(cp))) {
#   if(is.null(dp)) stop("in the multivariate case, 'dp' must be non-NULL")
#   info <-  sn.infoMv(dp=dp, x=x, y=y, w=w, type=type, norm2.tol=norm2.tol)
#   } else {
#   if(any(is.numeric(dp), is.numeric(cp)))
#   info <- sn.infoUv(dp=dp, cp=cp, x=x, y=y, w=w, penalty=penalty, 
#     type=type, norm2.tol = norm2.tol)
#   else stop("invalid input")
#   }
# return(info)
# } 
 
sn.infoUv <- function(dp=NULL, cp=NULL, x=NULL, y, w, penalty=NULL,   
                      norm2.tol=1e-6)
{# computes observed/expected Fisher information for univariate SN variates
  if(missing(y)) {y <- NULL; type <- "expected"} else type <- "observed"
  if(type == "observed") {if(!is.numeric(y)) stop("y is non-numeric")} 
  if(is.null(dp) & is.null(cp)) stop("either dp or cp must be set")
  if(!is.null(dp) & !is.null(cp)) stop("cannot set both dp and cp")
  if(missing(w)) w <- rep(1, max(NROW(cbind(x,y)),1)) 
  if(any(w != round(w)) | any(w<0))
    stop("weights must be non-negative integers")
  n <- length(w)
  nw <- sum(w)
  if(is.null(x)) {
    p <- 1
    wx <- w
    xx <- sum.x <- nw
    x <- matrix(1, nrow=n, ncol=1)
    }
  else { 
    p <- NCOL(x)
    # x <- matrix(x, n, p)
    wx <- w*x
    xx <- t(x) %*% (wx)
    sum.x <- matrix(colSums(wx))
    }
  x.names <- if(length(colnames(x)) == p) colnames(x)[2:p]  else
               { if(p==1) NULL else paste("x", 1L:(p-1), sep=".")}
  if(is.null(cp)) {
    if(length(dp) != (p+2)) stop("length(dp) must be equal to ncol(x)+2")
    if(is.null(names(dp))) names(dp) <- param.names("DP", "SN", p, x.names)
    cp <- dp2cpUv(dp, "SN")
    }
  if(is.null(dp)) {
    if(length(cp) != (p+2)) stop("length(cp) must be equal to ncol(x)+2")
    if(is.null(names(cp))) names(cp) <- param.names("CP", "SN", p, x.names)
    dp <- cp2dpUv(cp, "SN")
    }       
  penalty.fn <- if(is.null(penalty)) NULL else get(penalty, inherits=TRUE)
  omega <- dp[p+1]
  alpha <- dp[p+2]
  mu.z   <- sqrt(2/pi)*alpha/sqrt(1+alpha^2)
  sd.z   <- sqrt(1-mu.z^2)
  sigma  <- cp[p+1]
  gamma1 <- cp[p+2]
  R <- mu.z/sd.z
  T <- sqrt(2/pi-(1-2/pi)*R^2)
  Da.Dg <- 2*(T/(T*R)^2+(1-2/pi)/T^3)/(3*(4-pi))
  Dmu.z <- sqrt(2/pi)/(1+alpha^2)^1.5
  Dsd.z <- (-mu.z/sd.z)*Dmu.z
  Ddp.cp <- diag(p+2)
  Ddp.cp[1,p+1] <- (-R)
  Ddp.cp[1,p+2] <- (-sigma*R)/(3*gamma1)
  Ddp.cp[p+1,p+1] <- 1/sd.z
  Ddp.cp[p+1,p+2] <- (-sigma)* Dsd.z* Da.Dg/sd.z^2
  Ddp.cp[p+2,p+2] <- Da.Dg
  I.dp <- I.cp  <- matrix(NA,p+2,p+2)
  if(type == "observed"){
    score <- sn.pdev.gh(cp, x, y, w, penalty.fn, trace=FALSE, hessian=TRUE)/(-2)
    I.cp <- attr(score, "hessian")/2
    attr(score,"hessian") <- NULL
    Dcp.dp <- solve(Ddp.cp)
    I.dp <- force.symmetry(t(Dcp.dp) %*% I.cp %*% Dcp.dp)
    a.coef <- NULL
    asyvar.cp <- pd.solve(I.cp, silent=TRUE)
    if(is.null(asyvar.cp)) {
      asyvar.dp <- NULL
      not.mle <- TRUE} 
    else {
      not.mle <- (abs(sum(score * as.vector(asyvar.cp %*% score))) > norm2.tol)
      asyvar.dp <- pd.solve(I.dp, silent=TRUE)
      }
    if(not.mle) warning("parameters do not seem at MLE")  
    #--Iinfo.dp 2nd form 
    I2 <- matrix(NA,p+2,p+2)
    z <- (y - as.vector(x%*% dp[1:p]))/omega
    z1 <- zeta(1, alpha*z)
    z2 <- zeta(2, alpha*z)
    I2[1:p,1:p] <- t(wx) %*% ((1 - alpha^2*z2)*x)/omega^2
    I2[1:p,p+1] <- t(wx) %*% (2*z - alpha*z1 - alpha^2*z2*z)/omega^2
    I2[p+1,1:p] <- t(I2[1:p,p+1])
    I2[1:p,p+2] <- t(wx) %*% (z1 + alpha*z2*z)/omega
    I2[p+2,1:p] <- t(I2[1:p,p+2])
    I2[p+1,p+1] <- (-nw + 3*sum(w*z^2) -2*alpha*sum(w*z1*z)
                    -alpha^2*sum(w*z2*z^2))/omega^2 
    I2[p+1,p+2] <- I2[p+2,p+1] <- (sum(w*z*z1) + alpha*sum(w*z2*z^2))/omega
    I2[p+2,p+2] <- sum(-w*z2*z^2)
   }
  else { # type == "expected"
    I2 <- NULL
    if(abs(alpha) < 200) {
      f.a <- function(x, alpha, k) x^k * dsn(x,0,1,alpha) * zeta(1,alpha*x)^2
      err <- .Machine$double.eps^0.5
      a0 <- integrate(f.a, -Inf, Inf, alpha=alpha, k=0, rel.tol=err)$value
      a1 <- integrate(f.a, -Inf, Inf, alpha=alpha, k=1, rel.tol=err)$value
      a2 <- integrate(f.a, -Inf, Inf, alpha=alpha, k=2, rel.tol=err)$value
      }
    else {# approx of Bayes & Branco (2007) with multiplicative adjustment
      u <- 1 + 8*(alpha/pi)^2
      b <- sqrt(2/pi) 
      a0 <- 1.019149098 * b^2/sqrt(u)
      a1 <- 1.020466516 * (-alpha * b^3/sqrt(u^3*(1+alpha^2/u)))
      a2 <- 1.009258704 * b^2/sqrt(u)^3
      }
    a.coef <- c(a0, a1, a2)
    I.dp[1:p,1:p] <- xx * (1+alpha^2*a0)/omega^2  
    I.dp[p+1,p+1] <- nw * (2+alpha^2*a2)/omega^2
    I.dp[p+2,p+2] <- nw * a2
    I.dp[1:p,p+1] <- sum.x * (mu.z*(1+mu.z^2*pi/2)+alpha^2*a1)/omega^2
    I.dp[p+1,1:p] <- t(I.dp[1:p,p+1])
    I.dp[1:p,p+2] <- sum.x * (sqrt(2/pi)/(1+alpha^2)^1.5-alpha*a1)/omega
    I.dp[p+2,1:p] <- t(I.dp[1:p,p+2])
    I.dp[p+1,p+2] <- I.dp[p+2,p+1] <- nw*(-alpha*a2)/omega 
    eps <- 0.005
    if(abs(alpha) >  eps) 
      I.cp  <- force.symmetry(t(Ddp.cp) %*% I.dp %*% Ddp.cp)
    else{ 
      if(alpha == 0) 
        I.cp <- diag(c(1/omega^2, 2/omega^2, 1/6))
      else {       
        add <- c(rep(0,p+1), 3*eps)
        i1 <- sn.infoUv(dp=dp+add, x=x, w=w)
        i2 <- sn.infoUv(dp=dp-add, x=x, w=w)
        I.cp <- (i1$info.cp + i2$info.cp)/2
        }
      }
    score <- NULL
    asyvar.dp <- pd.solve(I.dp, silent=TRUE)
    asyvar.cp <- pd.solve(I.cp, silent=TRUE)
    }
  dimnames(I.dp) <- list(names(dp), names(dp))
  if(!is.null(I.cp)) dimnames(I.cp) <- list(names(cp), names(cp))
  aux <- list(Ddp.cp=Ddp.cp, a.coef=a.coef, score.cp=score)
  list(dp=dp, cp=cp, type=type, info.dp=I.dp, info.cp=I.cp, 
       asyvar.dp=asyvar.dp, asyvar.cp=asyvar.cp, aux=aux)
}

sn.infoMv <- function(dp, x=NULL, y, w, penalty=NULL, norm2.tol=1e-6)
{# computes observed/expected Fisher information matrix for multiv.SN variates
 # using results in Arellano-Valle & Azzalini (JMVA, 2008+erratum)
  type <- if(missing(y)) "expected" else "observed"
  if(type == "observed") {if(!is.matrix(y)) stop("y is not a matrix")}
  cp <- dp2cpMv(dp, "SN")
  d <- length(dp$alpha)
  d2 <- d*(d+1)/2
  if(missing(w)) w <- rep(1, max(NROW(cbind(x,y)),1))
  if(any(w != round(w)) | any(w<0))
    stop("weights must be non-negative integers")
  n <- length(w)
  nw <- sum(w)
  if(is.null(x)) {
    p <- 1
    xx <- sum.x <- nw
    x <- matrix(1, nrow=n, ncol=1)
    }
  else { 
    p <- NCOL(x)
    # x <- matrix(x, n, p)
    xx <- drop(t(x) %*% (w*x))
    sum.x <- drop(matrix(colSums(w*x)))
    }
  beta <- as.matrix(dp[[1]],p,d)
  Omega <- dp$Omega
  omega <- sqrt(diag(Omega))
  alpha <- dp$alpha
  eta   <- alpha/omega
  # vOmega <- Omega[lower.tri(Omega,TRUE)]
  Obar <- cov2cor(Omega)
  Obar.alpha <-  as.vector(Obar %*% alpha)
  alpha.star <- sqrt(sum(alpha * Obar.alpha)) 
  if(alpha.star < 1e-4) {
    warning("information matrix of multivariate SN not computed near alpha=0")
    return(NULL)
    }
  # delta.star <- alpha.star/sqrt(1+alpha.star^2)
  c1 <- sqrt(2/pi)/sqrt(1+alpha.star^2)
  c2 <- 1/(pi*sqrt(1+2*alpha.star^2))
  # theta <- c(beta,vOmega,eta)
  D <- duplicationMatrix(d)
  i1 <- 1:prod(dim(beta))
  i2 <- max(i1) + 1:(d*(d+1)/2)
  i3 <- max(i2) + 1:d
  # ind <- list(i1=i1, i2=i2, i3=i3)
  O.inv <- pd.solve(Omega, silent=TRUE)
  if(type == "observed"){ 
    y0 <- y - x %*% beta
    S0 <- t(y0) %*% (w*y0) / nw
    y0.eta <- as.vector(y0 %*% eta)
    z1 <- zeta(1, y0.eta) * w
    z2 <- (-zeta(2, y0.eta) * w)
    # Z2 <- diag(z2, n)
    S1 <- (O.inv %x% t(x)) %*% as.vector(w*y0)- (eta %x% t(x)) %*% z1
    S2 <- (nw/2) * t(D) %*% ((O.inv %x% O.inv) %*% as.vector(S0-Omega))
    S3 <- t(y0) %*% z1
    score <- c(S1,S2,S3)
    u  <- t(x) %*% z1
    U  <- t(x) %*% (z2 * y0)
    V  <- O.inv %*% (2*S0-Omega) %*% O.inv
    # terms as given in the last but one matrix of p.16 
    j11 <- O.inv %x% xx + outer(eta,eta) %x% (t(x) %*% (z2 *x) )
    j12 <- (O.inv %x% (t(x) %*% (w*y0) %*% O.inv))  %*% D
    j13 <- diag(d) %x% u - eta %x% U
    j22 <- (nw/2) * t(D) %*% (O.inv %x% V) %*% D
    j23 <- matrix(0, d*(d+1)/2, d)
    j33 <- t(y0) %*% (z2 * y0)            
    uaA.coef <- NULL
    }
  else { # expected information
    Omega.eta <- omega * Obar.alpha
    mu.c <- Omega.eta/alpha.star^2 
    Omega.c <- Omega - outer(Omega.eta, Omega.eta)/alpha.star^2 
    alpha.bar <- alpha.star/sqrt(1+2*alpha.star^2)
    ginvMills <- function(x, m=0, s=1)  
        # generalized inverse Mills ratio: \phi(x; m, s^2)/\Phi(x)
        exp(-0.5*((x-m)^2/s^2-x^2)+log(zeta(1,x))-log(s))
    fn.u <- function(x, sd, k) x^k * ginvMills(x,0,sd) 
    if(alpha.bar > 0) {
      err<- .Machine$double.eps^0.5
      u0 <- integrate(fn.u, -Inf, Inf, sd=alpha.bar, k=0, rel.tol=err)$value
      u1 <- integrate(fn.u, -Inf, Inf, sd=alpha.bar, k=1, rel.tol=err)$value
      u2 <- integrate(fn.u, -Inf, Inf, sd=alpha.bar, k=2, rel.tol=err)$value }
    else {u0 <- 2; u1<- u2 <- 0}
    a0 <- u0
    a1 <- u1 * mu.c
    A2 <- u2 * outer(mu.c, mu.c) + u0 * Omega.c                    # cfr (19)
    A1 <- (c1*(diag(d)-outer(eta,eta) %*% Omega/(1+alpha.star^2))
           - c2*outer(eta, a1))   # cfr line after (12)
    # terms as given in the last matrix of p.16
    j11 <- (O.inv + c2*a0*outer(eta,eta)) %x% xx
    j12 <- c1*(O.inv %x% outer(sum.x, eta)) %*% D
    j13 <- A1 %x% sum.x
    j22 <- 0.5*nw *t(D) %*% (O.inv %x% O.inv) %*% D
    j23 <- matrix(0, d*(d+1)/2, d)
    j33 <- nw *c2 * A2
    uaA.coef <- list(u0=u0, u1=u1, u2=u2, a1=a1, A1=A1, A2=A2)
    score <- NULL
    }
  I.theta <-rbind(cbind( j11,    j12,   j13),
                  cbind(t(j12),  j22,   j23),
                  cbind(t(j13), t(j23), j33))
  if(!is.null(penalty)) { 
    # penalization depends on blocks (2,3) of the parameter set only
    penalty.fn <- if(is.null(penalty)) NULL else get(penalty, inherits=TRUE) 
    penalty.theta <- function(theta23, penalty, d) {
      vOmega <- theta23[1:(d*(d+1)/2)]
      eta <- theta23[(d*(d+1)/2) + (1:d)]
      Omega <- vech2mat(vOmega)
      alpha <- eta *sqrt(diag(Omega))
      penalty(list(alpha=alpha, Omega=Omega))
      } 
    i23 <- c(i2,i3)
    theta23 <- c(Omega[lower.tri(Omega,TRUE)], eta) # beta does not enter here
    score[i23] <- (score[i23] - 
      numDeriv::grad(penalty.theta, theta23, penalty=penalty.fn, d=d))
    jQ <- numDeriv::hessian(penalty.theta, theta23, penalty=penalty.fn, d=d)
    I.theta[i23, i23] <- I.theta[i23, i23] + jQ
    }                 
  I.theta <- force.symmetry(I.theta, tol=1e3)
  inv_I.theta <- pd.solve(I.theta, silent=TRUE)
  if(is.null(inv_I.theta)) {
     warning("numerically unstable information matrix")
     return(NULL)
     }
  if(type == "observed" ) {
    score.norm2 <- sum(score * as.vector(inv_I.theta %*% score))
    if(score.norm2/d > norm2.tol) stop("'dp' does not seem to be at MLE")
    }
  D32 <- matrix(0,d, d2)
  tmp32 <- matrix(0,d^2,d^2)
  for(i in 1:d){
    Eii <- matrix(0,d,d)
    Eii[i,i] <- 1
    tmp32 <- tmp32 + Eii %x% Eii
    }
  D32 <- (-0.5)* (t(eta) %x% diag(1/omega^2, d,d)) %*% tmp32 %*% D
  # here we use the expression given in the notes, not in the paper
  Dlow <- cbind(matrix(0,d,d*p), D32, diag(1/omega,d,d))
  Dtheta.dp <- rbind(cbind(diag(d*p+d2), matrix(0,d*p+d2,d)), Dlow)
  I.dp <- t(Dtheta.dp) %*% I.theta %*% Dtheta.dp                     # cfr (14)
  I.dp <- force.symmetry(I.dp, tol=1e3)
  #
  # psi<- c(mu, vSigma, mu0)
  Sigma <- cp$var.cov
  sigma <- sqrt(diag(Sigma))
  Sigma.inv <- pd.solve(Sigma)
  mu0 <- c1* omega * Obar.alpha
  beta0.sq <- as.vector(t(mu0) %*% Sigma.inv %*% mu0)
  beta0 <- sqrt(beta0.sq)
  q1 <- 1/(c1*(1+beta0.sq))
  q2 <- 0.5*q1*(2*c1-q1)
  Dplus <- pd.solve(t(D) %*% D) %*% t(D)
  D23 <- Dplus %*% (diag(d) %x% mu0 + mu0 %x% diag(d))
  a <- as.vector(Sigma.inv %*% mu0)
  D32 <- t(-a) %x% (q1 * Sigma.inv - q1*q2*outer(a,a)) %*% D
  D33 <- q1 * Sigma.inv - 2*q1*q2*outer(a,a)
  one00 <- c(1,rep(0,p-1))
  Dtheta.psi <- rbind(
        cbind(diag(p*d),  matrix(0,p*d,d2), -diag(d) %x% one00),
        cbind(matrix(0,d2,p*d),  diag(d2),   D23),
        cbind(matrix(0,d,p*d),    D32,       D33))                # cfr (22a)
  mu0. <- mu0/(sigma*beta0)  # \bar{\mu}_0
  D32. <- matrix(0, d, d2)   # \tilde{D}_{32}
  for(i in 1:d)  {
    Eii <- matrix(0,d,d)
    Eii[i,i] <- 1
    D32. <- D32. + (1/sigma[i])*((t(mu0.) %*% Eii) %x% Eii) %*% D
    }
  D32. <- 0.5* beta0 * D32.
  D33. <- (2/(4-pi)) * diag(sigma/mu0.^2, d, d)/(3*beta0.sq)
  Dpsi.cp <- rbind(cbind(diag(p*d+d2), matrix(0,p*d+d2,d)), 
                   cbind(matrix(0,d,p*d), D32., D33.))            # cfr (22b)
  jacob <- Dtheta.psi %*% Dpsi.cp
  I.cp <- t(jacob) %*% I.theta %*% jacob                          # cfr (17)
  I.cp <- if(any(is.na(I.cp))) NULL else force.symmetry(I.cp)  
  asyvar.dp <- pd.solve(I.dp, silent=TRUE)
  if(is.null(asyvar.dp))  se.dp <- list(NULL) else {
    diags.dp <- sqrt(diag(asyvar.dp))
    se.beta <- matrix(diags.dp[1:(p*d)], p, d)
    se.diagOmega <- diags.dp[p*d + d2 +1 -rev(cumsum(1:d))]
    # se.omega <- se.Omega/(2*omega)
    se.alpha <- diags.dp[p*d +d2 +(1:d)]
    se.dp <- list(beta=se.beta, diagOmega=se.diagOmega, alpha=se.alpha)
    }
  asyvar.cp <- pd.solve(I.cp, silent=TRUE)
  if(is.null(asyvar.cp))  se.cp <- list(NULL) else {
    diags.cp <- sqrt(diag(asyvar.cp))
    se.beta <- matrix(diags.cp[1:(p*d)], p, d)
    se.diagSigma <- diags.cp[p*d + d2 +1 -rev(cumsum(1:d))]
    # se.sigma <- se.Sigma/(2*sigma)
    se.gamma1 <- diags.cp[p*d + d2 +(1:d)]
    se.cp <- list(beta=se.beta, var=se.diagSigma, gamma1=se.gamma1)
    }
  aux <- list(info.theta=I.theta, score.theta=score,
              Dtheta.dp=Dtheta.dp, Dpsi.cp=Dpsi.cp, Dtheta.psi=Dtheta.psi, 
              uaA.coef=uaA.coef)
  list(dp=dp, cp=cp, type=type, info.dp=I.dp, info.cp=I.cp, 
       asyvar.dp=asyvar.dp, asyvar.cp=asyvar.cp, 
       se.dp=se.dp, se.cp=se.cp, aux=aux)
}



msn.mle <- function(x, y, start=NULL, w, trace=FALSE, 
                opt.method=c("nlminb", "Nelder-Mead", "BFGS", "CG",  "SANN"),
                control=list() )
{
  y <- data.matrix(y)
  if(missing(x)) x <- rep(1,nrow(y))
    else {if(!is.numeric(x)) stop("x must be numeric")}
  if(missing(w)) w <- rep(1,nrow(y))
  opt.method <- match.arg(opt.method)
  x <- data.matrix(x) 
  d <- ncol(y)  
  n <- sum(w)
  p <- ncol(x)
  y.names <- dimnames(y)[[2]] 
  x.names <- dimnames(x)[[2]]
  if(is.null(start)) {
     fit0  <- lm.wfit(x, y, w, method="qr")
     beta  <- as.matrix(coef(fit0))
     res   <- resid(fit0)
     a     <- msn.moment.fit(res)
     Omega <- a$Omega
     omega <- a$omega
     alpha <- a$alpha
     if(!a$admissible) alpha<-alpha/(1+max(abs(alpha)))
     beta[1,] <- beta[1,]-omega*a$delta*sqrt(2/pi)  
     }
  else{
    beta  <- start[[1]] # start$beta
    Omega <- start$Omega
    alpha <- start$alpha
    omega <- sqrt(diag(Omega)) 
    }
  eta <-alpha/omega
  if(trace){ 
    cat("Initial parameters:\n")
    print(cbind(t(beta),eta,Omega))
    }
  param <- c(beta,eta)
  dev <- msn.dev(param, x, y, w)    
  if(opt.method == "nlminb") {
    opt <- nlminb(param, msn.dev, msn.dev.grad, control=control, x=x, y=y, 
              w=w, trace=trace)      
    opt$value <- opt$objective 
    }
  else opt <- optim(param, fn=msn.dev, gr=msn.dev.grad, method=opt.method,
                  control=control, x=x, y=y, w=w, trace=trace)    
  if(trace) {
    cat("Message from function", opt.method, ":", opt$message,"\n")
    cat("Output parameters " , format(opt$par), "\n")
    }
  logL <- opt$value/(-2) 
  beta <- matrix(opt$par[1:(p*d)],p,d)
  dimnames(beta)[2] <- list(y.names)
  dimnames(beta)[1] <- list(x.names)
  eta   <- opt$par[(p*d+1):(p*d+d)]
  xi    <- x %*% beta
  Omega <- t(y-xi) %*% (w*(y-xi))/n
  omega <- sqrt(diag(Omega))
  alpha <- eta*omega
  # param <- cbind(omega,alpha)
  dimnames(Omega) <- list(y.names,y.names)
  names(alpha) <- y.names
  alpha2 <- sum(eta * as.vector(Omega %*% eta))
  delta.star <- sqrt(alpha2/(1+alpha2))
  # dimnames(param)[1] <- list(y.names)
  dp  <- list(beta=beta, Omega=Omega, alpha=alpha)
  opt$method <- opt.method
  opt$called.by <- "msn.mle"
  aux <- list(alpha.star=sqrt(alpha2), delta.star=delta.star)
  list(call=match.call(), dp=dp, logL=logL, aux=aux, opt.method=opt)
}

 
msn.dev <- function(param, x, y, w, trace=FALSE)
{
  d <- ncol(y)
  if(missing(w)) w <- rep(1,nrow(y))
  n <- sum(w)
  p <- ncol(x)
  beta <- matrix(param[1:(p*d)],p,d)
  eta <- param[(p*d+1):(p*d+d)]
  y0 <- y-x %*% beta
  Omega <- (t(y0) %*% (y0*w))/n  
  D <- diag(qr(2*pi*Omega)[[1]])
  logDet <- sum(log(abs(D)))
  dev <- n*logDet - 2*sum(zeta(0, y0 %*% eta) * w) + n*d
  if(trace) { 
    cat("\nmsn.dev:",dev,"\n","parameters:"); 
    print(rbind(beta,eta))
    }
  dev
}

msn.dev.grad <- function(param, x, y, w, trace=FALSE)
{
  d <- ncol(y)
  if(missing(w)) w <- rep(1,nrow(y))
  n <- sum(w)
  p <- ncol(x)
  beta <- matrix(param[1:(p*d)],p,d)
  eta <- param[(p*d+1):(p*d+d)]
  y0 <- y-x %*% beta
  Omega <- (t(y0) %*% (w*y0))/n
  p1 <- zeta(1,as.vector(y0 %*% eta)) * w
  Omega.inv <- pd.solve(Omega, silent=TRUE)
  if(is.null(Omega.inv)) return(rep(NA, p*d+d))
  Dbeta <- (t(x) %*% (y0*w) %*% Omega.inv - outer(as.vector(t(x) %*% p1), eta))
  Deta <- as.vector(t(y0) %*% p1)
  if(trace){
    cat("gradient:\n")
    print(rbind(Dbeta,Deta))}
  -2*c(Dbeta,Deta)
}


msn.moment.fit <- function(y)
{# 31-12-1997: simple fit of MSN distribution usign moments
  y     <- as.matrix(y)
  k     <- ncol(y)
  m.y   <- apply(y, 2, mean)  
  var.y <- var(y)
  y0    <- (t(y) - m.y)/sqrt(diag(var.y))
  gamma1<- apply(y0^3, 1, mean)
  out   <- (abs(gamma1) > 0.99527)
  gamma1[out] <- sign(gamma1[out])*0.995
  a     <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^0.33333
  delta <- sqrt(pi/2)*a/sqrt(1+a^2)
  m.z   <- delta * sqrt(2/pi) 
  omega <- sqrt(diag(var.y)/(1-m.z^2))
  Omega <- var.y + outer(omega*m.z, omega*m.z) 
  xi    <- m.y-omega*m.z
  O.cor <- cov2cor(Omega)
  O.inv <- pd.solve(O.cor)
  tmp   <- as.vector(1 - t(delta) %*% O.inv %*% delta)
  if(tmp<=0) {tmp <- 0.0001; admissible <- FALSE} 
        else admissible <- TRUE
  alpha <- as.vector(O.inv %*% delta)/sqrt(tmp)
  list(xi=xi, Omega=Omega, alpha=alpha, Omega.cor=O.cor, omega=omega, 
       delta=delta, skewness=gamma1, admissible=admissible) 
}

  
st.mple <- function(x, y, dp=NULL, w, fixed.nu=NULL, symmetr=FALSE, 
  penalty=NULL, trace=FALSE, 
  opt.method=c("nlminb", "Nelder-Mead", "BFGS", "CG", "SANN"), control=list())
{ # MLE of DP for univariate ST distribution, allowing case symmetr[ic]=TRUE
  if(missing(y)) stop("required argument y is missing")
  if(!is.vector(y) | !is.numeric(y)) stop("argument y must be a numeric vector")
  x <- if(missing(x)) matrix(rep(1, length(y)), ncol = 1) else data.matrix(x)
  if(!is.matrix(x)) stop("argument x must be a matrix")
  y.name <- deparse(substitute(y))
  x.name <- deparse(substitute(x))
  if(any(x[,1] != 1)) stop("first column of x must have all 1's")
  if(symmetr && !is.null(penalty)) 
    stop("Penalized log-likelihood not allowed with constraint alpha=0")
  n <- length(y)
  p <- ncol(x)
  if(missing(w)) w <- rep(1,n)
  nw <- sum(w)
  if(is.null(dp)) {
    ls <- lm.wfit(x, y, w)
    res <- ls$residuals
    s <- sqrt(sum(w*res^2)/nw)
    gamma1 <- sum(w*res^3)/(nw*s^3)
    gamma2 <- sum(res^4)/(nw*s^4) - 3 
    cp <- c(ls$coef, s, gamma1, gamma2)
    dp <- st.cp2dp(cp, silent=TRUE)
    if(is.null(dp)) dp <- rep(NA,length(cp))
    if(any(is.na(dp))) dp <- c(cp[1:(p+1)], 0, 10)
    if(!is.null(fixed.nu)) dp <- dp[-length(dp)]
    if(symmetr) dp <- dp[-length(dp)]
    }
  else{ 
    if(length(dp) != (p+2-as.numeric(symmetr)+as.numeric(is.null(fixed.nu))))
       stop("arg 'dp' has wrong length")}
  if(trace) cat("dp (starting values) =", format(dp), "\n")
  tiny <- (.Machine$double.eps)^(0.25) 
  low.dp <- c(rep(-Inf, p), tiny, if(symmetr) NULL else -Inf,   
              if(is.null(fixed.nu)) tiny)
  high.dp <- c(rep(Inf, length(dp)))
  opt.method <- match.arg(opt.method)
  penalty.fn <- if(is.null(penalty)) NULL else get(penalty, inherits=TRUE) 
  if(opt.method == "nlminb") {
    opt <- nlminb(dp, objective=st.pdev, gradient=st.pdev.gh, 
           # do NOT set: hessian=st.dev.hessian, 
           lower=low.dp, upper=high.dp, control=control,
           x=x, y=y, w=w, fixed.nu=fixed.nu, symmetr=symmetr, 
           penalty=penalty.fn, trace=trace)
    opt$value <-  opt$objective
    }
  else {
    opt <- optim(dp, fn=st.pdev, gr=st.pdev.gh,  method = opt.method,
             # arguments lower & upper not used to allow all opt.method
             control = control,
             x=x, y=y, w=w, fixed.nu=fixed.nu, symmetr=symmetr, 
             penalty=penalty.fn, trace=trace)   
    }               
  dp <- opt$par
  opt$method <- opt.method
  opt$called.by <- "st.mple"
  dp. <- if(is.null(fixed.nu)) dp else c(dp, fixed.nu)  
  if(symmetr) dp. <- c(dp.[1:(p+1)], 0, dp.[length(dp.)])
  rv.comp <- c(TRUE, !symmetr, is.null(fixed.nu))
  names(dp) <- param.names("DP", "ST", p=p, x.names=colnames(x)[-1], rv.comp)
  names(dp.) <- param.names("DP", "ST", p=p, x.names=colnames(x)[-1])
  logL <- (-opt$value)/2
  boundary <- FALSE
  if(!symmetr) boundary <- as.logical(abs(dp[p+2]) > 1000) 
  if(is.null(fixed.nu)) boundary <- (boundary | dp[length(dp)] > 1e3)
  # AA, must improve this rule
  if(trace) {
     cat("Message from function", opt.method, ": ", opt$message, "\n")
     cat("estimates (dp):", dp, "\n")
     cat("log-likelihood:", logL, "\n")
     }
  list(call=match.call(), dp=dp, fixed.nu=fixed.nu, logL=logL, 
      dp.complete=dp., boundary=boundary, opt.method=opt)
}


st.pdev <- function(dp, x, y, w, fixed.nu=NULL, symmetr=FALSE, penalty=NULL, 
  trace=FALSE)
{ # computes "penalized deviance"=-2*(logL-Q) for ST
  p <- ncol(x)
  xi <- as.vector(x %*% matrix(dp[1:p],p,1))
  alpha <- if(symmetr) 0 else dp[p+2]
  nu <- if(is.null(fixed.nu)) dp[p+3-as.numeric(symmetr)] else fixed.nu
  if(dp[p+1] <= 0 | nu <= 0) return(NA)
  logL <- sum(w * dst(y, xi, dp[p+1], alpha, nu, log=TRUE))
  Q <- if(is.null(penalty)) 0 else penalty(dp[p+2], nu, der=0)
  if(trace) cat("st.pdev: (dp,pdev) =", format(c(dp, -2*(logL-Q))),"\n")
  return(-2 * (logL - Q))
}

st.pdev.gh <- function(dp, x, y, w, fixed.nu=NULL, symmetr=FALSE, 
   penalty=NULL, trace=FALSE, hessian=FALSE)
{ # computes gradient and hessian of (penalized) deviance for ST 
  p  <- ncol(x)
  n  <- nrow(x)
  beta  <- dp[1:p]
  omega <- dp[p+1]
  alpha <- if(symmetr) 0 else dp[p+2]
  j.nu <- p + 2 + as.numeric(!symmetr)
  nu <- if(is.null(fixed.nu)) dp[j.nu] else fixed.nu 
  npar <- p + 1 + as.numeric(!symmetr) + as.numeric(is.null(fixed.nu))
  score <- numeric(npar)
  xi <- as.vector(x %*% beta)
  z <- (y - xi)/omega
  nuz2 <- (nu + z^2)
  loro.tau <- sqrt((nu+1)/nuz2) 
  zt <- z * loro.tau
  log.pdf <- dt(alpha*zt, nu+1, log=TRUE)
  log.cdf <- pt(alpha*zt, nu+1, log.p=TRUE)
  cdf <- exp(log.cdf)
  loro.w <- exp(log.pdf - log.cdf)
  tw <- loro.tau * loro.w
  zwz2 <- z*(z^2-1)*loro.w/loro.tau
  wi.beta  <- z*loro.tau^2 - nu*alpha*tw/(nu+z^2)
  score[1:p] <- colSums(w*x*wi.beta)/omega
  score[p+1] <- sum(w * (-1 + zt^2 -alpha*nu*z*tw/(nu+z^2)))/omega
  if(!symmetr) score[p+2] <- sum(w*z*tw)
  if(is.null(fixed.nu)){
    fun.g <- function(x, nu1) dt(x,nu1) *
              (((nu1+1)*x^2)/(nu1*(nu1+x^2)) - log1p(x^2/nu1))
    int.g <- numeric(n)
    for (i in 1:n)
      int.g[i] <- integrate(fun.g, -Inf, alpha*zt[i], nu1=nu+1)$value
    score[j.nu] <- 0.5 * sum(w * (digamma(1+nu/2) -digamma(nu/2)
       - (2*nu+1)/(nu*(nu+1)) -log1p(z^2/nu) + zt^2/nu 
       + alpha*zwz2/(nu+z^2)^2 + int.g/cdf))
    }
  if(is.null(penalty)) { 
    Q <- 0
    attr(Q, "der1") <- rep(0,2)
    attr(Q, "der2") <- matrix(rep(0,4), 2, 2) }  
    else  {
    if(symmetr) stop("Penalized logL not allowed with constraint alpha=0") 
    Q <- penalty(alpha, nu, der=1+as.numeric(hessian))  
    }
  score[(p+2):(p+3)] <-  score[(p+2):(p+3)] - attr(Q, "der1") 
  score <- score[1:npar]
  gradient <- (-2)*score
  if(hessian){ 
    info <- matrix(NA, npar, npar) 
    w.z  <- (-nu*(nu+2)*alpha^2*z*loro.w/((nu+z^2*(1+alpha^2))*nuz2)
             -nu*alpha*loro.tau*loro.w^2/nuz2)
    w.alpha <- (-(nu+2)* alpha*z^2*loro.w/(nu+z^2*(1+alpha^2)) -zt*loro.w^2)
    S.z  <- (-z*loro.tau^2 + alpha*nu*tw/nuz2) 
    S.zz <- (2*zt^2/nuz2 - loro.tau^2 -3*alpha*nu*z*tw/nuz2^2
             +alpha*nu*loro.tau*w.z/nuz2) 
    info[1:p,1:p] <- t(-S.zz *x) %*% (w*x)/omega^2
    info[1:p,p+1] <- info[p+1,1:p] <- colSums(-w*(S.zz*z + S.z)*x)/omega^2
    info[p+1,p+1] <- -sum(w*(1 + z^2*S.zz + 2*z*S.z))/omega^2
    S.za <- nu*loro.tau*(loro.w +alpha*w.alpha)/nuz2
    if(!symmetr) {
      info[1:p,p+2] <- info[p+2,1:p] <- colSums(w*S.za*x)/omega
      info[p+1,p+2] <- info[p+2,p+1] <- sum(w*z*S.za)/omega
      info[p+2,p+2] <- sum(-w*zt*w.alpha) + attr(Q,"der2")[1,1]
      }
    if(is.null(fixed.nu)) {
      w.nu <- (0.5*loro.w*((nu+2)*(alpha*z)^2/((nu+z^2*(1+alpha^2))*nuz2)
               - log1p((alpha*z)^2/nuz2) - int.g/cdf)
               - 0.5*alpha*zwz2*loro.w/nuz2^2)
      S.znu <- (z*(1-z^2)/nuz2^2 + alpha*nu*loro.tau*w.nu/nuz2
                + alpha*(nu*(3*z^2-1)+2*z^2)*loro.w/(2*loro.tau*nuz2^3))
      info[1:p,j.nu] <- info[j.nu,1:p] <- colSums(w* S.znu*x)/omega
      info[p+1,j.nu] <- info[j.nu,p+1] <- sum(w*z*S.znu)/omega
    
      fun.b <- function(x, nu1) dt(x,nu1) *
                 (((nu1+1)*x^2)/(nu1*(nu1+x^2)) - log1p(x^2/nu1))^2
      fun.d <- function(x, nu1) dt(x,nu1) *
                 x^2*((nu1-1)*x^2-2*nu1)/(nu1^2*(nu1+x^2)^2)
      int.b <- int.d <- numeric(n)
      for (i in 1:n) {
        int.b[i] <- integrate(fun.b, -Inf, alpha*zt[i], nu1=nu+1)$value
        int.d[i] <- integrate(fun.d, -Inf, alpha*zt[i], nu1=nu+1)$value
        }
      info[j.nu,j.nu] <- -sum(w*( (trigamma(nu/2+1) - trigamma(nu/2))/4 
        + (2*nu^2+2*nu+1)/(2*(nu*(nu+1))^2) + z^2/(2*nu*nuz2)
        - z^2*(nu^2+2*nu+z^2)/(2*nu^2*nuz2^2)
        - alpha*zwz2*(z^2+4*nu+3)/(4*(nu+1)*nuz2^3)
        + alpha*z*(1-loro.tau^2)*w.nu/(2*loro.tau*nuz2)
        - (int.g/(2*cdf))^2 - alpha*zwz2*int.g/(4*cdf*nuz2^2)
        + (2*int.d + int.b)/(4*cdf)
        + (alpha*zwz2/(4*nuz2^2))*
          ((nu+2)*alpha^2*z^2/((nu+1)*(nu+z^2*(1+alpha^2))) 
            - log1p((alpha*z)^2/nuz2)) ))
      info[j.nu,j.nu] <- info[j.nu,j.nu] + attr(Q,"der2")[2,2]
      if(!symmetr) { 
        info[p+2,p+3] <- info[p+3,p+2] <- -sum(w*(0.5*zwz2/nuz2^2 + zt*w.nu))
        info[p+2,p+3] <- info[p+2,p+3] + attr(Q,"der2")[1,2]
        info[p+3,p+2] <- info[p+3,p+2] + attr(Q,"der2")[2,1]
        }
      }
    attr(gradient,"hessian") <- force.symmetry(2*info)
    if(trace) cat("Hessian matrix has been computed\n")
    }
  if(trace) cat("st.pdev.gh: gradient = ", format(gradient),"\n")
  return(gradient)
}

st.pdev.hessian <- function(dp, x, y, w, fixed.nu=NULL, symmetr=FALSE,  
     penalty = NULL, trace=FALSE)
  attr(st.pdev.gh(dp, x, y, w, fixed.nu, symmetr, penalty, trace, 
     hessian=TRUE), "hessian")


st.infoUv <- function(dp=NULL, cp=NULL, x=NULL, y, w, fixed.nu=NULL, 
   symmetr=FALSE, penalty=NULL, norm2.tol=1e-06)
{# computes observed Fisher information matrix for univariate ST variates
  if(missing(y)) stop("y is missing")
  if(!is.numeric(y)) stop("y is non-numeric")
  type <- "observed"
  if(is.null(dp) & is.null(cp)) stop("either dp or cp must be set")
  if(!is.null(dp) & !is.null(cp)) stop("cannot set both dp and cp")
  # if(is.null(cp)) cp <- st.dp2cp(c(dp, fixed.nu)) # completa DP se necessario
  if(is.null(dp)) dp <- st.cp2dp(cp) # AA, CP deve essere comunque completo
  if(missing(w)) w <- rep(1, max(nrow(cbind(x, y)), 1))
  if(any(w != round(w)) | any(w<0))
    stop("weights must be non-negative integers")
  npar <- length(dp)
  n <- length(w)
  nw <- sum(w)
  nu <- if(is.null(fixed.nu)) dp[npar] else fixed.nu
  if(is.null(x)) {
    n <- if(is.null(y)) 1 else NROW(y) 
    p <- 1
    xx <- sum.x <- nw 
    x <- matrix(1, nrow=n, ncol=1)
    }
  else { 
    p <- NCOL(x)
    # x <- matrix(x, n, p)
    xx <- t(x) %*% (w * x)
    sum.x <- matrix(colSums(x))
    }
  penalty.fn <- if(is.null(penalty)) NULL else get(penalty, inherits=TRUE)   
  score <- st.pdev.gh(dp, x, y, w, fixed.nu, symmetr, penalty.fn, trace=FALSE, 
              hessian=TRUE)
  I.dp <- attr(score, "hessian")/2
  if((d2 <- sum(score * as.vector(solve(I.dp) %*% score))) > norm2.tol*npar) {
    warning("'dp' does not seem to be at MLE; score not quite 0")
    cat("score(dp): ", score, "\n")
    cat("norm(score)^2:", d2,"\n")
    }
  attr(score, "hessian") <- NULL
  dimnames(I.dp) <- list(names(dp), names(dp))
  asyvar.dp <- pd.solve(I.dp, silent=TRUE)
  aux <- list(score.dp=score)
  if(nu > 4) {
    dp0 <- c(dp[1:(p+1)], if(symmetr) 0 else dp[p+2], if(is.null(fixed.nu)) nu)
    cp <- st.dp2cp(dp=dp0, cp.type="proper", fixed.nu=fixed.nu,
            upto=if(is.null(fixed.nu)) 4 else 3, jacobian=TRUE)
    Dcp.dp <- attr(cp, "jacobian")
    attr(cp, "jacobian") <- NULL
    ind <- c(1:(p+1), if(symmetr) NULL else (p+2), if(is.null(fixed.nu)) p+3) 
    Dcp.dp <- Dcp.dp[ind, ind]
    cp <- cp[ind]
    Ddp.cp <- solve(Dcp.dp)
    I.cp <- force.symmetry(t(Ddp.cp) %*% I.dp %*% Ddp.cp)
    dimnames(I.cp) <- list(names(cp), names(cp))
    asyvar.cp <- pd.solve(I.cp)
    aux$Dcp.dp <- Dcp.dp 
    aux$Ddp.cp <- Ddp.cp
    }  
  else  {
    I.cp <- NULL
    asyvar.cp <- NULL
    aux <- NULL
    }
  list(dp=dp, cp=cp, type=type, info.dp=I.dp, info.cp=I.cp, 
       asyvar.dp=asyvar.dp, asyvar.cp=asyvar.cp, aux=aux)
}


param.names <- function(param.type, family="SN", p=1, x.names=NULL, rv.comp)
{# NB: x.names= names of covariates except intercept, having length (p-1); 
 # rv.comp=random variable components (those not part of the linear predictor)
  if(!(param.type %in% c("DP","CP","pseudo-CP"))) stop("invalid param.type")
  if(!(family %in% c("SN", "ESN", "ST", "SC")))  stop("unknown family")
  if(p > 1  && (length(x.names) < (p-1)))
    x.names <- outer("x", as.character(1L:(p-1)), paste, sep=".")
  if(param.type == "DP"){
    name0 <-  if(p > 1) "(Intercept.DP)" else "xi"
    par.names <- c(name0, x.names, "omega", "alpha")
    if(family == "ESN") par.names <- c(par.names, "tau")
    if(family == "ST") par.names <- c(par.names, "nu")
    }
  if(param.type == "CP"){
    name0 <-  if(p > 1) "(Intercept.CP)" else "mean"
    par.names <- c(name0, x.names, "s.d.", "gamma1")
    if(family == "ESN") par.names <- c(par.names, "tau")
    if(family == "ST") par.names <- c(par.names, "gamma2")
    }
  if(param.type == "pseudo-CP"){
    if(!(family %in% c("ST", "SC"))) 
      stop("pseudo-CP makes sense only for ST and SC families")
    name0 <-  if(p > 1) "(Intercept.CP~)" else "mean~"
    par.names <- c(name0, x.names, "s.d.~", "gamma1~")
    if(family == "ST") par.names <- c(par.names, "gamma2~")
    }
  if(missing(rv.comp)) rv.comp <- rep(TRUE, length(par.names)-p)
  par.names[c(rep(TRUE,p), rv.comp)]
}


mst.mple <- function (x, y, start=NULL, w, fixed.nu = NULL, symmetr=FALSE,
                penalty=NULL, trace = FALSE, 
                opt.method = c("nlminb", "Nelder-Mead", "BFGS", "CG", "SANN"),
                control = list()) 
{
  opt.method <- match.arg(opt.method)
  if(missing(y)) stop("required argument y is missing")
  if(!is.matrix(y) | !is.numeric(y)) stop("argument y must be a numeric matrix") 
  y.name <- deparse(substitute(y))
  y.names <- dimnames(y)[[2]]
  n <- nrow(y)
  x <- if (missing(x)) matrix(rep(1, n), ncol = 1) else data.matrix(x)
  if (missing(w)) w <- rep(1, n)
  nw <- sum(w)
  x.names <- dimnames(x)[[2]]
  d <- ncol(y)
  p <- ncol(x)
  if (is.null(start)) {
    ls <- lm.wfit(x, y, w, singular.ok=FALSE)
    beta <- coef(ls)
	Omega <-  var(resid(ls))
	omega <- sqrt(diag(Omega))
	alpha <- rep(0, d)
    nu <- if(is.null(fixed.nu)) 8 else fixed.nu
    if (trace) cat("mst.mple: starting dp = (",
         c(beta, Omega[!upper.tri(Omega)], alpha, nu),  ")\n")
  }
  else {
    if (!is.null(fixed.nu))  start$nu <- fixed.nu
    if (all(names(start)[2:4] == c("Omega", "alpha", "nu"))) {
        beta <- start[[1]] # was start$beta
        Omega <- start$Omega
        alpha <- start$alpha
        nu <- start$nu
      }
    else stop("argument 'start' is not in the form that I expected")
    }
  if(symmetr) alpha <- rep(0,d) 
  param <- dplist2optpar(list(beta=beta, Omega=Omega, alpha=alpha))
  if(symmetr) param <- param[-(p*d + d*(d+1)/2 + (1:d))]
  if(is.null(fixed.nu)) param <- c(param, log(nu))
  if(!is.null(penalty)) penalty <- get(penalty, inherits=TRUE)
  opt.method <- match.arg(opt.method)
  if(opt.method == "nlminb") {
    opt <- nlminb(param, objective = mst.pdev, gradient = mst.pdev.grad, 
             control = control,  x = x, y = y, w = w, fixed.nu = fixed.nu, 
             symmetr=symmetr,  penalty=penalty, trace = trace)
    # info <- num.deriv2(opt$par, FUN="mst.dev.grad", X=X, y=y,
    #           w=w, fixed.nu = fixed.nu)/2
    opt$value <-  opt$objective
    }
  else {
      opt <- optim(param, fn = mst.pdev, gr = mst.pdev.grad, 
               method = opt.method, control = control, hessian = TRUE,
               x = x, y = y, w = w,   fixed.nu = fixed.nu, 
               symmetr=symmetr, penalty=penalty, trace = trace)
      # info <- opt$hessian/2
      }
  dev   <- opt$value
  param <- opt$par
  opt$method <- opt.method
  opt$called.by <- "mst.mple"
  if (trace) {
      cat("Message from optimization routine:", opt$message, "\n")
      cat("(penalized) deviance:", dev, "\n")
  }
  par <- opt$par 
  npar0 <- (p*d + d*(d+1)/2)
  vp <- par[1:npar0]
  dp.comp <- (1:2)
  if(symmetr) vp <- c(vp, rep(0,d)) else {
    vp <- c(vp, par[npar0 + (1:d)]); dp.comp <- (1:3)}
  if(is.null(fixed.nu)) {
    vp <- c(vp, par[length(par)])
    dp.comp <- c(dp.comp,4)}
  dp.list <- optpar2dplist(vp, d, p, x.names, y.names)
  dp <- dp.complete <- dp.list$dp
  if(symmetr) dp.complete$alpha <- rep(0, d)
  if(!is.null(fixed.nu)) dp.complete$nu <- fixed.nu
  alpha2 <- sum(dp$alpha * as.vector(cov2cor(dp$Omega) %*% dp$alpha))
  delta.star <- sqrt(alpha2/(1+alpha2))
  dp <- dp[dp.comp]
  aux <- list(fixed.nu=fixed.nu, symmetr=symmetr, alpha.star=sqrt(alpha2), 
              delta.star=delta.star)
  boundary <- ((1 - delta.star) < .Machine$double.eps^(1/4))
  if(is.null(fixed.nu)) boundary <- (boundary | dp$nu > 1e3)     
  list(call=match.call(), dp=dp, dp.complete=dp.complete, logL=dev/(-2),
    boundary=boundary, aux=aux, opt.method = opt)
}


mst.pdev <- function(param, x, y, w, fixed.nu=NULL, symmetr=FALSE, 
   penalty=NULL, trace=FALSE)
{
  if(missing(w)) w <- rep(1,nrow(y))
  d <- ncol(y)
  p <- ncol(x)
  npar0 <- (p*d + d*(d+1)/2)
  param1 <- c(param[1:npar0], if(symmetr) rep(0, d) else param[npar0+(1:d)], 
    if(is.null(fixed.nu)) param[length(param)])
  dp.list <- optpar2dplist(param1, d, p)
  dp <- dp.list$dp
  nu <- if(is.null(fixed.nu)) dp$nu else fixed.nu
  logL <- sum(w * dmst(y, x %*% dp$beta, dp$Omega, dp$alpha, nu, log=TRUE))
  Q <- if(is.null(penalty)) 0 else 
       penalty(list(alpha=dp$alpha, Omega.bar=cov2cor(dp$Omega)), nu, der=0)
  pdev <- (-2) * (logL - Q)
  if(trace) cat("mst.pdev: ", pdev, "\nparam:", format(param), "\n")
  pdev
}


mst.pdev.grad <- function(param, x, y, w, fixed.nu=NULL, symmetr=FALSE, 
  penalty=NULL, trace=FALSE)
{ # based on Appendix B of Azzalini & Capitanio (2003, arXiv-0911.2342)
  # except for a few quite patent typos (transposed matrices, etc) 
  d <- ncol(y)
  p   <- ncol(x)
  beta<- matrix(param[1:(p*d)],p,d)
  D  <- exp(-2*param[(p*d+1):(p*d+d)])
  A  <- diag(d)
  i0 <- p*d + d*(d+1)/2
  if(d>1) A[!lower.tri(A,diag=TRUE)] <- param[(p*d+d+1):i0]
  eta  <- if(symmetr) rep(0,d) else param[(i0+1):(i0+d)]
  nu   <- if(is.null(fixed.nu)) exp(param[length(param)]) else fixed.nu
  Oinv <- t(A) %*% diag(D,d,d) %*% A
  u    <- y - x %*% beta
  u.w  <- u * w
  Q    <- as.vector(rowSums((u %*% Oinv) * u.w))
  L    <- as.vector(u.w %*% eta)
  sf   <- if(nu < 1e4) sqrt((nu+d)/(nu+Q)) else sqrt((1+d/nu)/(1+Q/nu))
  t.   <- L*sf                                     # t(L,Q,nu) in \S 5.1
  # dlogft<- (-0.5)*(1+d/nu)/(1+Q/nu)              # \tilde{g}_Q
  dlogft <- (-0.5)*sf^2                            # \tilde{g}_Q, again
  dt.dL <- sf                                      # \dot{t}_L
  dt.dQ <- (-0.5)*L*sf/(Q+nu)                      # \dot{t}_Q
  logT. <- pt(t., nu+d, log.p=TRUE)
  dlogT.<- exp(dt(t., nu+d, log=TRUE) - logT.)     # \tilde{T}_1
  Dbeta <- (-2* t(x) %*% (u.w*dlogft) %*% Oinv 
            - outer(as.vector(t(x) %*% (dlogT. * dt.dL* w)), eta)
            - 2* t(x) %*% (dlogT.* dt.dQ * u.w) %*% Oinv )
  Deta  <- colSums(dlogT.*sf*u.w)
  if(d>1) {
     M  <- 2*( diag(D,d,d) %*% A %*% t(u * dlogft
               + u * dlogT. * dt.dQ) %*% u.w)
     DA <- M[!lower.tri(M,diag=TRUE)]
     }
  else DA<- NULL
  M <- (A %*% t(u*dlogft + u*dlogT.*dt.dQ) %*% u.w %*% t(A))
  if(d>1) DD <- diag(M) + 0.5*sum(w)/D
     else DD <- as.vector(M + 0.5*sum(w)/D) 
  grad <- (-2) * c(Dbeta, DD*(-2*D), DA, if(!symmetr) Deta)
  if(is.null(fixed.nu)) {
    df0 <- min(nu, 1e8)
    if(df0 < 10000){
       diff.digamma <- digamma((df0+d)/2) - digamma(df0/2)
       log1Q<- log(1+Q/df0)
     }
    else
      {
       diff.digamma <- log1p(d/df0)
       log1Q <- log1p(Q/df0)
      }
    dlogft.ddf <- 0.5 * (diff.digamma - d/df0
                        + (1+d/df0)*Q/((1+Q/df0)*df0) - log1Q)
    eps   <- 1.0e-4
    df1 <- df0 + eps
    sf1 <- if(df0 < 1e4) sqrt((df1+d)/(Q+df1)) else sqrt((1+d/df1)/(1+Q/df1))
    logT.eps <- pt(L*sf1, df1+d, log.p=TRUE)
    dlogT.ddf <- (logT.eps-logT.)/eps
    Ddf   <- sum((dlogft.ddf + dlogT.ddf)*w)
    grad <- c(grad, -2*Ddf*df0)
    }
  if(!is.null(penalty)) { 
    if(symmetr) stop("penalized log-likelihood not allowed when alpha=0")
    Ainv <- backsolve(A, diag(d))
    Omega <- Ainv %*% diag(1/D,d,d) %*% t(Ainv)
    omega <- diag(Omega)
    alpha <- eta*omega
    Q <- Qpenalty(list(alpha, cov2cor(Omega)), nu, der=1)
    comp <-  1:(length(alpha)+is.null(fixed.nu))
    Qder <- attr(Q, "der1") * c(1/omega, 1)[comp] 
    # gradient for transformed variable (alpha --> eta)
    grad <- grad + 2*c(rep(0, p*d + d*(d+1)/2),  Qder)
    }
  if(trace) cat("mst.pdev.grad: norm is ", format(sqrt(sum(grad^2))), "\n")  
  return(grad)
}
 

mst.theta.jacobian <- function(theta, p, d, cp.type="proper") 
{ # jacobian matrices associated to transformations from
  # theta=c(beta, vech(Omega), eta, nu) to DP, CP and other parameterizations
  cp.type <- match.arg(cp.type, c("proper", "pseudo"))
  k1 <- p * d
  k2 <- k1 + d*(d+1)/2
  k3 <- k2 + d
  k4 <- k3 + 1
  if(length(theta) != k4) stop("mismatch in the arguments")
  block1 <- 1:k1
  block2 <- (k1+1):k2
  block3 <- (k2+1):k3
  block4 <- k4
  beta  <- matrix(theta[block1], p, d)
  Omega <- vech2mat(theta[block2]) 
  Omega.inv <- pd.solve(Omega)
  eta <- theta[block3]
  nu <- theta[block4]
  a.incr <- if(cp.type=="proper") rep(0,4) else 1:4
  omega <- sqrt(diag(Omega))
  alpha <- eta*omega
  # delta <- delta.etc(alpha, Omega)$delta
  D <- duplicationMatrix(d)
  P <- matrix(0, d^2, d^2)
  for (i in 1:d) {
    Eii <- matrix(0,d,d)
    Eii[i,i] <- 1
    P <- P + Eii %x% Eii
    }
  omega <- sqrt(diag(Omega))
  d <- length(omega)
  delta.plus <- delta.etc(alpha, Omega)
  delta <- delta.plus$delta
  delta.sq <- (delta.plus$delta.star)^2
  alpha.sq <- (delta.plus$alpha.star)^2
  a  <- function(nu) nu/(nu-2)
  u  <- function(nu) 0.5*(1/nu + digamma((nu-1)/2) - digamma(nu/2))
  c1 <- function(nu) b(nu)/sqrt(1 + alpha.sq)
  q1 <- function(nu) a(nu)/(c1(nu)*(1 + beta0.sq(nu)))
  q2 <- function(nu) q1(nu)*(2*c1(nu) - q1(nu))/(2*a(nu))
  beta0.sq <- function(nu) # beta0.sq = sum(mu0 * Sigma.inv_mu0) =
    b(nu)^2 * alpha.sq/(a(nu)+(a(nu)-b(nu)^2)*alpha.sq)
  #-- Dtheta.dp = D_{DP}\theta
  Dtheta.dp <- diag(k4)
  diag(Dtheta.dp)[block3] <- 1/omega
  Deta.vOmega <- (-0.5)* (t(eta) %x% diag(1/omega^2, d, d)) %*% P %*% D
  Dtheta.dp[block3, block2] <- Deta.vOmega
  #
  mu0 <- function(nu) omega * b(nu) * delta
  Sigma.etc <- function(nu) {
    mu0. <- mu0(nu)
    Omega.inv_mu0 <- as.vector(Omega.inv %*% mu0.)
    Sigma <- a(nu)*Omega - outer(mu0., mu0.)
    sigma <- sqrt(diag(Sigma))
    tmp <- a(nu) - sum(mu0. *Omega.inv_mu0)
    Sigma.inv_mu0 <- Omega.inv_mu0/tmp
    Sigma.inv <- (Omega.inv + outer(Omega.inv_mu0, Omega.inv_mu0)/tmp)/a(nu)
    list(Sigma=Sigma, Sigma.inv=Sigma.inv, Sigma.inv_mu0=Sigma.inv_mu0,    
         sigma=sigma)
    }
  Dq1.nu <- function(nu){
    beta0_sq <- beta0.sq(nu)
    (-2/(nu-2)^2 -a(nu)*(b(nu)^2*u(nu)+beta0_sq/((nu-2)^2*(1+beta0_sq)))
            /c1(nu)^2)/(c1(nu)*(1+beta0_sq))
    } 
  # blocks for D_{\Psi}\theta   
  Dplus <- solve(t(D)%*% D) %*% t(D)
  DvOmega.vSigma <- function(nu) diag(d*(d+1)/2)/a(nu)
  DvOmega.mu0 <- function(nu)
    Dplus %*% (diag(d) %x% mu0(nu) + mu0(nu) %x% diag(d))/a(nu) 
  DvOmega.nu <- function(nu){
    s <- Sigma.etc(nu)    
    2*vech(s$Sigma + outer(mu0(nu), mu0(nu)))/nu^2
    }
  Deta.vSigma <- function(nu) { 
    S <- Sigma.etc(nu)
    t(-S$Sigma.inv_mu0) %x%  (q1(nu)* S$Sigma.inv -
         q1(nu) * q2(nu) *outer(S$Sigma.inv_mu0, S$Sigma.inv_mu0)) %*% D
    }
  Deta.mu0 <- function(nu) {
    S <- Sigma.etc(nu)
    q1(nu) * (S$Sigma.inv - 2*q2(nu)*outer(S$Sigma.inv_mu0, S$Sigma.inv_mu0))
    } 
  Deta.nu <- function(nu) Dq1.nu(nu) * Sigma.etc(nu)$Sigma.inv_mu0   
  #-- Dtheta.phi(phi)= D_{\Psi}\theta
  one00 <- c(1,rep(0,p-1))
  Dtheta.phi <- diag(k4)
  Dtheta.phi[block1, block3] <- -diag(d) %x% one00
  Dtheta.phi[block2, block2] <- DvOmega.vSigma(nu+a.incr[2])
  Dtheta.phi[block2, block3] <- DvOmega.mu0(nu+a.incr[2])
  Dtheta.phi[block2, block4] <- DvOmega.nu(nu+a.incr[2])
  Dtheta.phi[block3, block2] <- Deta.vSigma(nu+a.incr[2])
  Dtheta.phi[block3, block3] <- Deta.mu0(nu+a.incr[2])
  Dtheta.phi[block3, block4] <- Deta.nu(nu +a.incr[2])
  #
  # blocks for D_{\Psi}CP    
  Dgamma2M.misc <- function(nu){
    beta0_sq <- beta0.sq(nu)
    s <- Sigma.etc(nu)
    nu.34 <- (nu-3)*(nu-4)
    tmp2 <- ( (d+2)/nu.34
      + beta0_sq * (2*nu/((nu-3)*b(nu)^2) - (3*(nu-3)^2-6)/nu.34 ))
    Dgamma2M.mu0 <- as.vector(8 * tmp2 *  t(s$Sigma.inv_mu0))
    Dgamma2M.vSigma <- (-4 * tmp2) * as.vector(( t(s$Sigma.inv_mu0) %x%  
                       t(s$Sigma.inv_mu0)) %*% D)
    R <- b(nu)^2*delta.sq*(nu-2)/nu
    R1R <- R/(1-R) 
    PDgamma2.nu <- (-2*d*(d+2)/(nu-4)^2 -4*((2*nu-7)/nu.34^2) *R1R*(2/(1-R)+d)
       +2*(2*((nu-3)-nu*(1+2*(nu-3)*u(nu)))/((nu-3)*b(nu))^2 
       +(3*nu^2-22*nu+41)/nu.34^2)*R1R^2)      #\ref{f:partial_gamma2.nu}  
    list(Dgamma2M.vSigma=Dgamma2M.vSigma, Dgamma2M.mu0=Dgamma2M.mu0,  
         PDgamma2.nu=PDgamma2.nu)    
  }
  Dgamma1.misc <- function(nu) {
    sigma <- Sigma.etc(nu)$sigma
    lambda <- mu0(nu)/sigma
    g.nu <- 3/(nu-3)
    h.nu <- 1 + nu*(1-1/b(nu)^2)/(nu-3)
    Q <- g.nu*diag(d) + 3*h.nu*diag(lambda^2)
    Dgamma1.vOmega <- (t(-lambda/2) %x% (Q %*% diag(1/sigma^2,d))) %*% P %*% D
    Dgamma1.mu0 <- Q %*% diag(1/sigma,d)                            # K_{33}
    Dgamma1.nu <- (-3*lambda/(nu-3)^2 + (-3*(1-1/b(nu)^2)/(nu-3)^2 
                       + 2*nu*u(nu)/((nu-3)*b(nu)^2))*lambda^3)     # K_{34}
    list(Dgamma1.vOmega=Dgamma1.vOmega, Dgamma1.mu0=Dgamma1.mu0, 
         Dgamma1.nu=Dgamma1.nu)                   
    }
  #
  #--
  # Dcp.phi(phi) = D_{\Psi}(CP) [in the notes] = D_{\phi}\bar\rho [paper]
  #
  Dcp.phi <- diag(k4)
  K3 <- Dgamma1.misc(nu+a.incr[3])
  K4 <- Dgamma2M.misc(nu+a.incr[4])
  Dcp.phi[block3,block2] <- K3$Dgamma1.vOmega 
  Dcp.phi[block3,block3] <- K3$Dgamma1.mu0 
  Dcp.phi[block3,block4] <- K3$Dgamma1.nu
  Dcp.phi[block4,block2] <- K4$Dgamma2M.vSigma
  Dcp.phi[block4,block3] <- K4$Dgamma2M.mu0 
  Dcp.phi[block4,block4] <- K4$PDgamma2.nu  
  #
  # Dtheta.cp <- Dtheta.phi %*% solve(Dcp.phi)
  list(Dtheta.dp=Dtheta.dp, Dtheta.cp= Dtheta.phi %*% solve(Dcp.phi),
       Dtheta.phi=Dtheta.phi, Dcp.phi=Dcp.phi)
  } 
# 
mst.vdp2vcp <- function(vdp, p, d, cp.type="proper") 
{ # vdp = c(betaDP, vech(Omega), alpha, nu), 
  # vcp=(betaCP, vech(Sigma), gamma1, gamma2M)
  # d=ncol(y), p=ncol(x)
  beta <- matrix(vdp[1:(p*d)], p, d)
  vOmega <- vdp[(p*d+1):(p*d+d*(d+1)/2)]
  Omega <- vech2mat(vOmega)
  # omega <- sqrt(diag(Omega))
  alpha <- vdp[(p*d+d*(d+1)/2+1):(p*d+d*(d+1)/2+d)]
  nu <- vdp[p*d+d*(d+1)/2+d+1]
  dp <- list(beta=beta, Omega=Omega, alpha=alpha, nu=nu)
  cp <- mst.dp2cp(dp, cp.type=cp.type)
  c(cp[[1]], vech(cp[[2]]), cp[[3]], cp[[4]])
}  
#
mst.logL <- function(vdp,  X, y, dp=TRUE, penalty=NULL) 
{ # calcola logL rispetto a DP (se dp=TRUE) oppure a theta (se dp=FALSE),
  # con eventuale inclusione del termine 'penalty' se presente;
  # funziona non solo per ST, ma anche per SN ponendo dp$nu=Inf
  n <- nrow(y)
  d <- ncol(y)
  if(missing(X)) X <- matrix(1,n,1)
  p <- ncol(X)
  beta <- matrix(vdp[1:(p*d)], p, d)
  vOmega <- vdp[(p*d+1):(p*d+d*(d+1)/2)]
  Omega <- vech2mat(vOmega)
  # if(any(eigen(Omega)$values <= 0)) return(NA)
  if(any(diag(Omega) <= 0)) return(-Inf)
  omega <- sqrt(diag(Omega))
  tmp <- vdp[(p*d+d*(d+1)/2+1):(p*d+d*(d+1)/2+d)]
  alpha <- if(dp) tmp else tmp*omega
  nu <- vdp[p*d+d*(d+1)/2+d+1]
  if(nu <= 0) return(-Inf)
  Q <- if(is.null(penalty)) 0 else penalty(list(alpha, cov2cor(Omega)), nu)
  sum(dmst(y, X %*% beta, Omega, alpha, nu, log=TRUE)) - Q
}


st.infoMv <- function(dp, x=NULL, y, w, fixed.nu=NULL, symmetr=FALSE, 
   penalty=NULL, norm2.tol=1e-06)
{# Computes observed Fisher information matrices for multiv.ST distribution
 # using expressions of score function of Arellano-Valle (2010, Metron),
 # followed by numerical differentiation. Expected info matrix not implemented.
 # Info matrices are computed for DP, CP and pseudo-CP
  if(missing(y)) stop("missing y")
  if(!is.matrix(y)) stop("y is not matrix")
  type <- "observed"
  d <- ncol(dp$Omega)
  d2 <- d*(d+1)/2
  if(missing(w)) w <- rep(1, nrow(cbind(x,y)))
  if(any(w != round(w)) || any(w<0))
    stop("weights must be non-negative integers")
  n <- length(w)
  nw <- sum(w)
  if(is.null(x)) {
    p <- 1
    xx <- sum.x <- nw
    x <- matrix(1, nrow=n, ncol=1)
    }
  else { 
    p <- NCOL(x)
    # x <- matrix(x, n, p)
    xx <- drop(t(x) %*% (w*x))
    sum.x <- drop(matrix(colSums(w*x)))
    }
  beta <- as.matrix(dp[[1]], p, d)
  Omega <- dp[[2]]
  omega <- sqrt(diag(Omega))
  alpha <- if(symmetr) rep(0,d) else dp$alpha
  eta   <- alpha/omega
  nu <- if(is.null(fixed.nu)) dp$nu else fixed.nu 
  dp1 <- list(beta=beta, Omega=Omega, alpha=alpha, nu=nu)
  Obar <- cov2cor(Omega)
  Obar.alpha <-  as.vector(Obar %*% alpha)
  alpha.star <- sqrt(sum(alpha * Obar.alpha)) # =\sqrt{\eta\T\Omega\eta}
  theta <- as.numeric(c(beta, vech(Omega), eta, nu))
  vdp <- as.numeric(c(beta, vech(Omega), alpha, nu))
  penalty.fn <- if(is.null(penalty)) NULL else get(penalty, inherits=TRUE) 
  H <- numDeriv::hessian(mst.logL, vdp, X=x, y=y, dp=TRUE, penalty=penalty.fn)
  J <- mst.theta.jacobian(theta, p=NCOL(x), d=NCOL(y))
  s <- 1:(length(theta) - as.numeric(!is.null(fixed.nu)))
  I.dp <- force.symmetry(-H[s,s])
  J1 <- solve(J$Dtheta.dp[s,s])
  I.theta <- force.symmetry(t(J1) %*% I.dp %*% J1)
  asyvar.dp <- pd.solve(I.dp, silent=TRUE)
  if(is.null(asyvar.dp)) { 
    warning("Condition 'information_matrix > 0' fails, DP seems not at MLE")
    se.dp <- list(NULL)  
    }
  else {
    diags.dp <- sqrt(diag(asyvar.dp))
    se.beta <- matrix(diags.dp[1:(p*d)], p, d)
    se.diagOmega <- diags.dp[p*d + d2 +1 - rev(cumsum(1:d))]
    se.alpha <- diags.dp[p*d +d2 +(1:d)]
    se.dp <- list(beta=se.beta, diagOmega=se.diagOmega, alpha=se.alpha)
    if(is.null(fixed.nu)) se.dp$nu<- diags.dp[p*d +d2 + d +1] 
    }
  if(nu>4) {
    cp <- mst.dp2cp(dp, cp.type="proper", fixed.nu=fixed.nu, symmetr=symmetr)
    I.cp <- force.symmetry(t(J$Dtheta.cp[s,s]) %*% I.theta %*% J$Dtheta.cp[s,s])
    asyvar.cp <- pd.solve(I.cp, silent=TRUE)
    if(is.null(asyvar.cp)) { 
      se.cp <- list(NULL)  
      }
    else {
      diags.cp <- sqrt(diag(asyvar.cp))
      se.beta <- matrix(diags.cp[1:(p*d)], p, d)
      se.diagSigma <- diags.cp[p*d + d2 +1 - rev(cumsum(1:d))]
      # se.sigma <- se.Sigma/(2*sigma)
      se.gamma1 <- if(!symmetr)  diags.cp[p*d + d2 +(1:d)] else NULL
      se.cp <- list(beta=se.beta, var=se.diagSigma, gamma1=se.gamma1)
      if(is.null(fixed.nu)) se.cp$gamma2 <- diags.cp[p*d +d2 + d +1] 
      }} 
  else 
    I.cp <- asyvar.cp <- se.cp <- cp <- NULL  
  if(is.null(asyvar.dp)) { 
    asyvar.pcp  <-  NULL
    se.pcp <- list(NULL)  
    Jp <- NULL
    }
  else {
    vdp1 <- as.numeric(c(dp1[[1]], vech(dp1[[2]]), dp1[[3]], dp1[[4]]))
    Jp <- numDeriv::jacobian(mst.vdp2vcp, vdp1, p=ncol(x), d=ncol(y), 
            cp.type="pseudo")
    asyvar.pcp <- (Jp[s,s]) %*% asyvar.dp %*% t(Jp[s,s])
    diags.pcp <- sqrt(diag(asyvar.pcp))
    se.beta <- matrix(diags.pcp[1:(p*d)], p, d)
    se.diagSigma <- diags.pcp[p*d + d2 +1 - rev(cumsum(1:d))]
    # se.sigma <- se.Sigma/(2*sigma)
    se.gamma1 <- if(!symmetr) diags.pcp[p*d + d2 +(1:d)] else NULL
    se.pcp <- list(beta=se.beta, var=se.diagSigma, gamma1=se.gamma1)
    if(is.null(fixed.nu)) se.pcp$gamma2 <- diags.pcp[p*d +d2 + d +1] 
    }
  aux <- list(Dpseudocp.dp=Jp[s,s]) 
  list(dp=dp, cp=cp, type=type, info.dp=I.dp, info.cp=I.cp, 
    asyvar.dp=asyvar.dp, asyvar.cp=asyvar.cp, asyvar.p_cp=asyvar.pcp,
    se.dp=se.dp, se.cp=se.cp, se.p_cp=se.pcp,  aux=aux)
}


sn.mple <- function(x, y, cp=NULL, w, penalty=NULL, trace=FALSE, 
  opt.method=c("nlminb", "Nelder-Mead", "BFGS", "CG", "SANN"), control=list()) 
{# MPLE for CP of univariate SN (not intendend for ESN)
  y <- drop(y)
  n <- length(y)
  if (missing(x))
    x <- matrix(rep(1,n), nrow=n, ncol=1)
  else
    if (is.null(n <- nrow(x)))  stop("'x' must be a matrix") 
  if (n == 0) stop("0-row design matrix cases")
  if (missing(w)) w <- rep(1,n)
  if(length(w) != n)  stop("incompatible dimensions")
  y.name <- deparse(substitute(y))
  x.name <- deparse(substitute(x)) 
  p <- ncol(x)
  max.gamma1 <- 0.5*(4-pi)*(2/(pi-2))^1.5 - (.Machine$double.eps)^(1/4)
  if(is.null(cp)) {
    qr.x <- qr(x)
    s <- sqrt(sum(qr.resid(qr.x, y)^2)/n)
    gamma1 <- sum(qr.resid(qr.x, y)^3)/(n*s^3)
    if(abs(gamma1) > max.gamma1) gamma1 <- sign(gamma1)*0.9*max.gamma1
    cp <- as.numeric(c(qr.coef(qr.x, y), s, gamma1))
    }
  else{ 
    if(length(cp)!= (p+2)) stop("ncol(x)+2 != length(cp)")}
  penalty.fn <- if(is.null(penalty)) NULL else get(penalty, inherits=TRUE)  
  opt.method <- match.arg(opt.method)  
  if(opt.method == "nlminb") {  
    opt <- nlminb(cp, objective=sn.pdev, 
           gradient=sn.pdev.gh, hessian=sn.pdev.hessian, 
           lower=c(-rep(Inf,p), sqrt(.Machine$double.eps), -max.gamma1), 
           upper=c(rep(Inf,p), Inf, max.gamma1), control=control,
           x=x, y=y, w=w, penalty=penalty.fn, trace=trace)
    opt$value <- opt$objective
    }
  else {
    opt <- optim(cp, fn=sn.pdev, gr=sn.pdev.gh,  
             method = opt.method, control = control,
             # lower & upper not used to allow all opt.method             
             x=x, y=y, w=w, penalty=penalty.fn, trace=trace)   
    } 
  cp <- opt$par
  names(cp) <- param.names("CP", "SN", p, colnames(x)[-1])
  logL <- (-opt$value)/2
  boundary <- as.logical(abs(cp[p+2]) >= max.gamma1)
  if(trace) {
    cat("Message from function", opt.method, ": ", opt$message, "\n")
    cat("estimates (cp):", cp, "\n")
    cat("(penalized) log-likelihood:", logL, "\n")
    }
  opt$method <- opt.method
  opt$called.by <- "sn.mple"  
  list(call=match.call(), cp=cp, logL=logL, boundary=boundary, opt.method=opt)
}


sn.pdev <- function(cp, x, y, w, penalty=NULL, trace=FALSE)
{ # "penalized deviance"=-2*(logL-Q) for centred parameters of SN distribution
  p <- ncol(x)
  if(abs(cp[p+2])> 0.9952717) return(Inf)
  if(missing(w)) w <- rep(1, length(y))
  if(any(w < 0)) stop("weights must be non-negative")
  dp <- cp2dpUv(cp, "SN") 
  xi <- as.vector(x %*% as.matrix(dp[1:p]))
  if(dp[p+1] <= 0) return(NA)
  logL <- sum(w * dsn(y, xi, dp[p+1], dp[p+2], log=TRUE))
  Q <- if(is.null(penalty)) 0 else penalty(dp[p+2], der=0)
  if(trace) cat("sn.pdev: (cp,pdev) =", format(c(cp, -2*(logL-Q))),"\n")
  return(-2 * (logL - Q))
}


sn.pdev.gh <- function(cp, x, y, w, penalty=NULL, trace=FALSE, hessian=FALSE)
{ # computes gradient and hessian of pdev=-2*(logL-Q) for centred parameters  
  p  <- ncol(x)
  n  <- nrow(x)
  if(abs(cp[p+2]) > 0.9952717) return(rep(NA,p+2))
  if(missing(w)) w <- rep(1,n)
  if(any(w < 0)) stop("weights must be non-negative")
  score <- rep(NA,p+2)
  info  <- matrix(NA,p+2,p+2)
  beta  <- cp[1:p]
  sigma <- cp[p+1]
  gamma1 <- cp[p+2]
  nw <- sum(w)
  dp <-  cp2dpUv(cp, "SN") 
  lambda <- dp[p+2]
  mu <- as.vector(x %*% as.matrix(beta))
  d  <- y-mu
  r  <- d/sigma
  mu.z<- lambda*sqrt(2/(pi*(1+lambda^2)))
  sd.z<- sqrt(1-mu.z^2)
  z  <- mu.z+sd.z*r
  p1 <- as.vector(zeta(1,lambda*z))
  p2 <- as.vector(zeta(2,lambda*z))
  omega<- sigma/sd.z
  af    <- lambda*p1-mu.z
  Dmu.z <- sqrt(2/pi)/(1+lambda^2)^1.5
  Dsd.z <- (-mu.z/sd.z)*Dmu.z
  Dz    <- Dmu.z + r*Dsd.z
  DDmu.z<- (-3)*mu.z/(1+lambda^2)^2
  DDsd.z<- -((Dmu.z*sd.z-mu.z*Dsd.z)*Dmu.z/sd.z^2+mu.z*DDmu.z/sd.z)
  DDz   <- DDmu.z + r*DDsd.z
  score[1:p] <- omega^(-2) * t(x) %*% as.matrix(w*(y-mu-omega*af))
  score[p+1] <- (-nw)/sigma + sd.z*sum(w*d*(z-p1*lambda))/sigma^2
  score.l <- nw*Dsd.z/sd.z - sum(w*z*Dz) + sum(w*p1*(z+lambda*Dz))
  if(!is.null(penalty)) {
    Q <- penalty(lambda, der=2)
    score.l <- (score.l - attr(Q, "der1"))
    }
  Dg.Dl <- 1.5*(4-pi)*mu.z^2 * (Dmu.z*sd.z - mu.z*Dsd.z)/sd.z^4
  R <- mu.z/sd.z
  T <- sqrt(2/pi-(1-2/pi)*R^2)
  Dl.Dg <- 2*(T/(T*R)^2+(1-2/pi)/T^3)/(3*(4-pi))
  R. <- 2/(3*R^2 * (4-pi))
  T. <- (-R)*R.*(1-2/pi)/T
  DDl.Dg <- (-2/(3*(4-pi))) * (T./(R*T)^2+2*R./(T*R^3)+3*(1-2/pi)*T./T^4)
  score[p+2] <- score.l/Dg.Dl  # convert deriv wrt lamda to gamma1 
  gradient <- (-2)*score
  if(hessian){ # info = -(second deriv of logL)
     info[1:p,1:p] <- omega^(-2) * t(x) %*% (w*(1-lambda^2*p2)*x)
     info[1:p,p+1] <- info[p+1,1:p] <- 
            sd.z* t(x) %*% as.matrix(w*(z-lambda*p1)+ w*d*(1-lambda^2*p2)*
            sd.z/sigma)/sigma^2
     info[p+1,p+1] <- (-nw)/sigma^2 + 2*sd.z*sum(w*d*(z-lambda*p1))/sigma^3 +
            sd.z^2*sum(w*d*(1-lambda^2*p2)*d)/sigma^4
     info[1:p,p+2] <- info[p+2,1:p] <- t(x) %*% (w*
            (-2*Dsd.z*d/omega+Dsd.z*af+sd.z*(p1+lambda*p2*(z+lambda*Dz)
            -Dmu.z)))/sigma 
     info[p+1,p+2] <- info[p+2,p+1] <- 
            -sum(w*d*(Dsd.z*(z-lambda*p1)+sd.z*(Dz-p1-p2*lambda*(z+lambda*Dz))
             ))/sigma^2
     info[p+2,p+2] <- (nw*(-DDsd.z*sd.z+Dsd.z^2)/sd.z^2+sum(w*(Dz^2+z*DDz)) -
            sum(w*p2*(z+lambda*Dz)^2)- sum(w*p1*(2*Dz+lambda*DDz)))
     if(!is.null(penalty)) info[p+2,p+2] <- info[p+2,p+2] + attr(Q, "der2")
     info[p+2,] <- info[p+2,]/Dg.Dl # convert info wrt lambda to gamma1 
     info[,p+2] <- info[,p+2]*Dl.Dg # an equivalent form of the above
     info[p+2,p+2] <- info[p+2,p+2] - score.l*DDl.Dg
     attr(gradient,"hessian") <- force.symmetry(2*info)
     }
  if(trace) cat("sn.pdev.gh: gradient = ", format(gradient),"\n")
  return(gradient)
}

sn.pdev.hessian <- function(cp, x, y, w, penalty=NULL, trace=FALSE)
{
  gh <- sn.pdev.gh(cp, x, y, w, penalty=penalty, trace=trace, hessian=TRUE)
  attr(gh, "hessian")
}    
 

Qpenalty <- function(alpha_etc, nu=NULL, der=0)
{# 'standard' penalty function of logL, possibly with derivatives
  e1 <- e1. <- 1/3
  e2 <- e2. <- 0.2854166
  if(!is.null(nu)) if(nu<Inf) { 
    g <- 0.57721
    e1 <- e1. * (nu+2)*(nu+3)/(nu+1)^2
    e2 <- e2. * (1 + 4/(nu+g))
    } else nu <- NULL
  c1 <- 1/(4*e2)
  c2 <- e2/e1 
  if(is.vector(alpha_etc) && length(alpha_etc)==1) {
    alpha<- alpha_etc
    Obar.alpha <- alpha
    alpha2 <- alpha^2
    }
  else {
    if(!is.list(alpha_etc)) stop("wrong argument alpha_etc")
    alpha <- alpha_etc[[1]]
    Omega.bar <- alpha_etc[[2]]
    if(any(dim(Omega.bar) != length(alpha))) stop("dimension mismatch")
    Obar.alpha <- as.vector(Omega.bar %*% alpha)
    alpha2 <- sum(alpha* Obar.alpha)
    }
  Q <- c1 * log(1 + c2* alpha2)
  if(der==0) return(Q)
  der1 <- 2*c1*c2*Obar.alpha/(1+ c2*alpha2)   
  if(!is.null(nu)) { 
    h <- (nu+g)*(nu+2)*(nu+3)
    dc1.dnu <- 1/(e2.*(nu+g+4)^2)
    tmp <- ((nu+1)^2 + 2*(nu+1)*(nu+g+4)) * h - (nu+1)^2*(nu+g+4)*(
            (nu+2)*(nu+3)+ (nu+2)*(nu+g)+(nu+3)*(nu+g))
    dc2.dnu <- 3*e2.*tmp/h^2
    der1 <- c(der1, Q*dc1.dnu/c1+ c1*alpha2*dc2.dnu/(1+c2*alpha2)) 
    }
  attr(Q, "der1") <- der1
  if(der==2) {
    attr(Q, "der2") <- if(is.null(nu)) 
      2*c1*c2*(1-c2*alpha^2)/(1+c2*alpha^2)^2   else 
    { 
      # Qdash <- function(x)  attr(Qpenalty(x[1], x[2], der=1), "der1") 
      # H <- jacobian(Qdash, c(alpha,nu))
      Q.fn <- function(x) Qpenalty(x[1], x[2], der=0)
      numDeriv::hessian(Q.fn, c(alpha, nu))
    }
  }
  return(Q)
}

MPpenalty <- function(alpha, der=0) 
{# penalty function associated to "matching prior" of Cabras et al.(SJS, 2012)
  a <- sn.infoUv(dp=c(0,1,alpha))$aux$a.coef
  a0 <- a[1]
  a1 <- a[2]
  a2 <- a[3]
  A <- 1+alpha^2
  num <- (a2*A^2*(pi*(1+a0*alpha^4) + alpha^2*(pi*(1+a0)-4))
          +2*sqrt(2*pi)*a1*alpha*A^1.5 - pi*a1^2*alpha^2*A^3 -2)
  den <- (pi*A^3*(2+alpha^2*(2*a0+a2)+ alpha^4*(a0*a2-a1^2))
         -2*(alpha+2*alpha^3)^2 
         -2*sqrt(2*pi)*a1*alpha^3*sqrt(A)*(1+3*alpha^2+2*alpha^4))
  prior <- sqrt(num/den)
  penalty <- -log(prior)
  if(der > 0) attr(penalty,"der1") <- numDeriv::grad(MPpenalty, alpha)
  if(der > 1) attr(penalty,"der2") <- numDeriv::hessian(MPpenalty, alpha)
  return(penalty)
}


msn.mple <- function(x, y, start=NULL, w, trace=FALSE, penalty=NULL,
                opt.method=c("nlminb", "Nelder-Mead", "BFGS", "CG",  "SANN"),
                control=list() )
{
  y <- data.matrix(y)
  if(missing(x)) x <- rep(1,nrow(y))
    else {if(!is.numeric(x)) stop("x must be numeric")}
  if(missing(w)) w <- rep(1,nrow(y))
  opt.method <- match.arg(opt.method)
  x <- data.matrix(x) 
  d <- ncol(y)  
  n <- sum(w)
  p <- ncol(x)
  y.names <- dimnames(y)[[2]] 
  x.names <- dimnames(x)[[2]]
  if(is.null(start))  start <- msn.mle(x, y, NULL, w)$dp
  if(trace){  
    cat("msn.mple initial parameters:\n")
    print(cbind(t(start[[1]]), start$Omega, start$alpha))
    }
  param <- dplist2optpar(start) 
  if(!is.null(penalty)) penalty <- get(penalty, inherits=TRUE)
  opt.method <- match.arg(opt.method)
  if(opt.method == "nlminb"){
    opt <- nlminb(param, msn.pdev, # msn.pdev.grad, 
                control=control, x=x, y=y, w=w, penalty=penalty, trace=trace)
    opt$value<- opt$objective 
    }
  else{
   opt <- optim(param, fn=msn.pdev, method=opt.method,
               control=control, x=x, y=y, w=w, penalty=penalty, trace=trace)
   }
  if(trace) 
    cat(paste("Message from optimization routine:", opt$message,"\n"))
  logL <- opt$value/(-2) 
  dp.list <- optpar2dplist(opt$par, d, p)
  beta <- dp.list$beta
  dimnames(beta)[2] <- list(y.names)
  dimnames(beta)[1] <- list(x.names)
  Omega <- dp.list$Omega
  alpha <- dp.list$alpha
  dimnames(Omega) <- list(y.names,y.names)
  names(alpha) <- y.names
  alpha2 <- sum(alpha * as.vector(cov2cor(Omega) %*% alpha))
  delta.star <- sqrt(alpha2/(1+alpha2))
  dp  <- list(beta=beta, Omega=Omega, alpha=alpha)
  opt$method <- opt.method
  opt$called.by <- "msn.mple"
  aux <- list(penalty=penalty, alpha.star=sqrt(alpha2), delta.star=delta.star)
  list(call=match.call(), dp=dp, logL=logL, aux=aux, opt.method=opt)
}

msn.pdev <- function(param, x, y, w, penalty=NULL, trace=FALSE)
{ # -2*(profile.logL - Q)
  d <- ncol(y)
  if(missing(w)) w <- rep(1, nrow(y))
  n <- sum(w)
  p <- ncol(x)
  dp. <- optpar2dplist(param, d=ncol(y), p=ncol(x))
  logL <- sum(w * dmsn(y, x %*% dp.$beta, dp.$Omega, dp.$alpha, log=TRUE))
  Q <- if(is.null(penalty)) 0 else penalty(list(dp.$alpha,dp.$Omega), der=0)
  pdev <- (-2)*(logL-Q)
  if(trace) cat("opt param:", format(param), "\nmsn.pdev:", format(pdev),"\n")
  return(pdev)
}
 
optpar2dplist <- function(param, d, p, x.names=NULL, y.names=NULL)
{# convert vector form of optimization parameters to DP list;
 # output includes inverse(Omega) and its log determinant 
  beta <- matrix(param[1:(p * d)], p, d)
  D <- exp(-2 * param[(p * d + 1):(p * d + d)])
  A <- diag(d)
  i0 <- p*d + d*(d+1)/2
  if(d>1)  A[!lower.tri(A,diag=TRUE)] <- param[(p*d+d+1):i0]
  eta <- param[(i0 + 1):(i0 + d)]
  nu <- if(length(param) == (i0 + d + 1)) exp(param[i0 + d + 1]) else NULL
  Oinv <- t(A) %*% diag(D,d,d) %*% A
  # Omega <- pd.solve(Oinv)
  Ainv <- backsolve(A, diag(d))
  Omega <- Ainv %*% diag(1/D,d,d) %*% t(Ainv)
  Omega <- (Omega + t(Omega))/2
  omega <- sqrt(diag(Omega))
  alpha <- eta * omega
  dimnames(beta) <- list(x.names, y.names)
  dimnames(Omega) <- list(y.names, y.names)
  if (length(y.names) > 0) names(alpha) <- y.names
  dp <- list(beta=beta, Omega=Omega, alpha=alpha)
  if(!is.null(nu)) dp$nu <- nu
  list(dp=dp, beta=beta, Omega=Omega, alpha=alpha, nu=nu, Omega.inv=Oinv,
     log.det=sum(log(D)))
}

dplist2optpar <- function(dp,  Omega.inv=NULL)
{# convert DP list to vector form of optimization parameters 
  beta <- dp[[1]]
  Omega <- dp[[2]]
  alpha <- dp[[3]]
  d <- length(alpha)
  nu <- if(is.null(dp$nu)) NULL else dp$null
  eta <- alpha/sqrt(diag(Omega))
  Oinv <- if(is.null(Omega.inv)) pd.solve(Omega) else Omega.inv
  if(is.null(Oinv)) stop("matrix Omega not symmetric positive definite")
  upper <- chol(Oinv)
  D <- diag(upper)
  A <- upper/D
  D <- D^2
  param <- if(d > 1)  c(beta, -log(D)/2, A[!lower.tri(A, diag = TRUE)], eta)
     else c(beta, -log(D)/2, eta)
  if(!is.null(dp$nu))  param <- c(param, log(dp$nu)) 
  param <- as.numeric(param)
  attr(param, 'ind') <- cumsum(c(length(beta), d, d*(d-1)/2, d, length(dp$nu)))
  return(param) 
}  


force.symmetry <- function(x, tol=10*sqrt(.Machine$double.eps)) 
{
  if(!is.matrix(x)) stop("x must be a matrix")
  # err <- abs(x-t(x))
  err <- abs(x-t(x))/(1+abs(x))
  max.err <- max(err/(1+err))
  if(max.err > tol) warning("matrix seems not symmetric")
  if(max.err > 100*tol) stop("this matrix really seems not symmetric")
  return((x + t(x))/2)
}
 
duplicationMatrix <- duplication_matrix <- function (n=1)
{# translated by AA from Octave code written by <Kurt.Hornik@wu-wien.ac.at>
  if ( (n<1) |  (round (n) != n) ) stop ("n must be a positive integer")
  d <- matrix (0, n * n, n * (n + 1) / 2)
  ## KH: It is clearly possible to make this a LOT faster!
  count = 0
  for (j in 1 : n){
    d [(j - 1) * n + j, count + j] = 1
    if(j<n) {
      for (i in (j + 1) : n){
       d [(j - 1) * n + i, count + i] = 1
       d [(i - 1) * n + j, count + i] = 1
    }}
    count = count + n - j
  }
  return(d)
}

vech <- function(A) if(is.matrix(A)) A[lower.tri(A, diag=TRUE)] else 
  stop("argument 'A' must be a matrix")

vech2mat <- function(v) 
{# inverse function of vech(A)
  if(mode(v) != "numeric" | !is.vector(v)) stop("wrong type of argument")
  n <- round((-1 + sqrt(1 + 8*length(v)))/2)
  if(length(v) != n*(n+1)/2) stop("wrong length of vector 'v'") 
  A <- matrix(0, n, n)
  A[lower.tri(A,TRUE)] <- v
  return(A + t(A) - diag(diag(A), n))
}

#-------

plot.SECdistrUv <- function(x, range, probs, main, npt=251, data=NULL, ...)
{# plot density of object "SECdistrUv"
  obj <- x
  lc.family <- tolower(slot(obj, "family"))
  if(lc.family == "esn") lc.family <- "sn"
  dp <- slot(obj, "dp")
  d.fn <- get(paste("d", lc.family, sep=""), inherits = TRUE)
  q.fn <- get(paste("q", lc.family, sep=""), inherits = TRUE)
  if(missing(probs)) probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  if(!is.numeric(data)) data <- NULL
  q <- if(is.null(probs)) NULL else q.fn(probs, dp=dp)
  if(missing(range)) {
    q0 <- q.fn(c(0.05, 0.2, 0.8, 0.95), dp=dp)
    dq <- diff(q0)
    range <- c(q0[1]- 1.5*dq[1], q0[length(q0)] +  1.5*dq[length(dq)])
    if(!is.null(q)) {
      range[1] <- min(range[1], min(q))
      range[2] <- max(range[2], max(q))
      }
    if(!is.null(data)) {
      range[1] <- min(range[1], min(data))
      range[2] <- max(range[2], max(data))
    }}
  dots <- list(...)
  nmdots <- names(dots)  
  topline <- if(obj@name == "") "" else	
    paste("Probability density of ", obj@name, "\n", sep="")
  if(missing(main))  main <- paste(topline, obj@family," distribution",
     ", dp = (",paste(format(dp), collapse=", "),")",sep="") 
  mar <- if ("mar" %in% nmdots) dots$mar else NULL
  if (is.null(mar)) {
    mar <- c(4.5, 4.5, 4, 2)
    if (is.null(main))  mar[3L] <- 2
    }     
  omar <- par()$mar
  on.exit(par(omar))
  par(mar=mar)    
  x <- seq(min(range), max(range), length=npt)
  pdf <- d.fn(x, dp=dp)
  xLab <- if("xlab" %in% nmdots) dots$xlab else slot(obj, "name")
  yLab <- if("ylab" %in% nmdots) dots$ylab else "probability density"
  yLim <- if("ylim" %in% nmdots) dots$ylim else c(0, max(pdf))
  plot(x, pdf, type="n", xlab=xLab, ylab=yLab, ylim=yLim)
  lines(x, pdf, ...)
  abline(h=0, lty=2, col="gray50")
  if(!is.null(q)) {
    points(q, rep(max(pdf)/100,length(q)), cex=0.75, col="gray50", pch=6)
    text(q, par()$usr[3]/2, format(probs), cex=0.75, col="gray50")
    }
  if(!is.null(data)) {
    side <- if(is.null(probs)) 1 else 3  
    rug(data, side=side, ticksize = 0.02, col="gray50")
    }
  if (!is.null(main)) {
    font.main <- if("font.main" %in% nmdots) dots$font.main else par("font.main") 
    cex.main <- if("cex.main" %in% nmdots) dots$cex.main else par("cex.main") 
    title(main, line=2, cex.main=cex.main, font.main=font.main)
    }    
  invisible(list(object=obj, x=x, density=pdf))
}
      

plot.SECdistrMv <- function(x, range, probs, npt, landmarks="auto",
       main, comp, compLabs, data = NULL, data.par=NULL, gap = 0.5, ...) 
{# plot density of object of class "SECdistrMv"  
  obj <- x
  if(slot(obj, "class") != "SECdistrMv") stop("object of wrong class")
  compNames <- slot(obj, "compNames")
  d <- length(compNames)
  if(missing(probs)) probs <- c(0.25, 0.50, 0.75, 0.95)    
  if(any(probs <= 0 || probs >= 1)) stop("probs must be within (0,1)") 
  if(sum(probs > 0 && probs < 1) == 0) stop("invalid probs") 
  if(missing(npt)) npt <- rep(101, d)
  if(missing(main))  { main <- if(d==2) 
      paste("Density function of", slot(obj, "name")) else
      paste("Bivariate densities of", slot(obj, "name")) }
  if(missing(comp)) comp <- seq(1,d)      
  if(missing(compLabs)) compLabs <- compNames
  if(length(compLabs) != d) stop("wrong length of 'compLabs' or 'comp' vector")
  family <- toupper(obj@family)
  lc.family <- tolower(family)
  if(lc.family == "esn") lc.family <- "sn"
  dp <- slot(obj, "dp")
  if(missing(range)) {
    range <- matrix(NA,2,d)
    q.fn <- get(paste("q", lc.family, sep=""), inherits=TRUE)
    for(j in 1:d) {
      marg <- marginalSECdistr(obj, comp=j, drop=TRUE)
      q <- q.fn(c(0.05, 0.25, 0.75, 0.95), dp=marg@dp)
      dq <- diff(q)
      range[,j] <- c(q[1] - 1.5*dq[1], q[length(q)] + 1.5*dq[length(dq)])
      if(!is.null(data)) {
        range[1,j] <- min(range[1,j], min(data[,j]))
        range[2,j] <- max(range[2,j], max(data[,j]))
       }}
    }
  dots <- list(...)
  nmdots <- names(dots)  
  if(d == 1) {  
    message("Since dimension=1, plot as a univariate distribution")
    objUv <- marginalSECdistr(obj, comp=1, drop=TRUE)
    out <- plot(objUv, data=data, ...)
    }
  if(d == 2) 
      out <- list(object=obj,
                  plot=plot.SECdistrBv(x, range, probs, npt, compNames,  
                           compLabs, landmarks, data, data.par, main, ...))
  if(d > 2) {
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) 
      text(x, y, txt, cex = cex, font = font)
    localAxis <- function(side, x, y, xpd, bg, main,  oma, ...) {
      if (side%%2 == 1) Axis(x, side = side, xpd = NA,  ...) else 
      Axis(y, side = side, xpd = NA, ...)
      }
    localPlot <- function(..., oma, font.main, cex.main) plot.SECdistrBv(...)
    text.diag.panel <- compLabs
    oma <- if ("oma" %in% nmdots) dots$oma else NULL
    if (is.null(oma)) {
      oma <- c(4, 4, 4, 4)
      if (!is.null(main))  oma[3L] <- 6
      }    
    opar <- par(mfrow = c(length(comp), length(comp)), 
                mar = rep(c(gap,gap/2), each=2), oma=oma)
    on.exit(par(opar))
    out <- list(object=obj)
    count <- 1
    for (i in comp) 
      for (j in comp) {
        count <- count + 1
        if(i == j) {
          plot(1, type="n", xlab="", ylab="", axes=FALSE)
          text(1, 1, text.diag.panel[i], cex=2)
          box()
          out[[count]] <- list()
          names(out)[count] <- paste("diagonal component", compNames[i])
        } else {
        ji <- c(j,i) 
        marg <- marginalSECdistr(obj, comp=ji) 
        out[[count]] <- localPlot(x=marg, range=range[,ji], probs=probs,
             npt=npt[ji], compNames= compNames[ji], compLabs=compLabs[ji], 
             landmarks=landmarks, data=data[,ji],  data.par=data.par, 
             main="", yaxt="n", xaxt="n", ...)   
        names(out)[count] <- paste("plot of components (", j, ",", i, ")")
        # if(i==comp[1]) {axis(3); if(j==length(comp)) axis(4)}
        # if(j==comp[1]) {axis(2); if(i==length(comp)) axis(1)}
        if(i==comp[1]) axis(3) ; if(j==length(comp)) axis(4)
        if(j==comp[1]) axis(2) ; if(i==length(comp)) axis(1)    
        box() }
      }
    par(new = FALSE)
    if (!is.null(main)) {
      font.main <- if ("font.main" %in% nmdots) 
         dots$font.main else par("font.main") 
      cex.main <- if ("cex.main" %in% nmdots) 
         dots$cex.main  else par("cex.main") 
      mtext(main, side=3, TRUE, line=5, outer = TRUE, at=NA, cex=cex.main, 
            font=font.main, adj=0.5)
      }}
  invisible(out)
}

plot.SECdistrBv <- function(x, range, probs, npt=rep(101,2), compNames,   
            compLabs, landmarks, data=NULL, data.par, main, ...)
{# plot BiVariate SEC distribution
  obj <- x
  dp <- slot(obj, "dp")
  family <- slot(obj, "family")
  lc.family <- tolower(family)
  if(lc.family == "esn") lc.family <- "sn"
  d.fn <- get(paste("dm", lc.family, sep=""), inherits=TRUE) # density funct
  n1 <- npt[1]
  n2 <- npt[2]
  x1 <- seq(min(range[,1]), max(range[,1]), length=n1)
  x2 <- seq(min(range[,2]), max(range[,2]), length=n2)
  x1.x2 <- cbind(rep(x1, n2), as.vector(matrix(x2, n1, n2, byrow=TRUE)))
  X <- matrix(x1.x2, n1 * n2, 2, byrow = FALSE)
  pdf <- matrix(d.fn(X, dp=dp), n1, n2)
  Omega <- dp[[2]]
  Omega.bar <- cov2cor(Omega)
  alpha <- dp[[3]]
  alpha.star <- sqrt(sum(alpha * as.vector(Omega.bar %*% alpha)))
  omega <- sqrt(diag(Omega))
  if(lc.family == "sn") {
    k.tau <- if (length(dp) == 4) (zeta(2,dp[[4]])*pi)^2/4 else 1
    log.levels <- (log(1-probs) - log(2*pi)- 0.5*log(1-Omega.bar[1,2]^2)
                   + k.tau * log(1+exp(-1.544/alpha.star))) - sum(log(omega))
    }
  if(lc.family == "st" | lc.family == "sc") {
    nu <- if(lc.family == "st") obj@dp[[4]] else 1
    l.nu <- (-1.3/nu - 4.93)
    if(alpha.star > 0) {
      h <- 100 * log(exp(((1.005*alpha.star-0.045)* l.nu -1.5)/alpha.star)+1) 
      K <-  h *(1.005*alpha.star-0.1)*(1+nu)/(alpha.star * nu) }  else K <- 0
    qF <- qf(probs, 2, nu)
    log.levels <- (lgamma(nu/2+1) -lgamma(nu/2) - log(pi*nu) 
           -0.5*log(1-Omega.bar[1,2]^2) - (nu/2+1)*log(2*qF/nu + 1)  + K
           -sum(log(omega)))
    } 
  oo <- options()
  options(warn=-1)
  # contour(x1, x2, pdf, levels=exp(log.levels), 
  #  labels=paste("p=", as.character(probs), sep=""),
  #   main=main, xlab=compLabs[1], ylab=compLabs[2], ...)
  plot(x1, x2, type="n", main=main, xlab=compLabs[1], ylab=compLabs[2], ...)
  if(!is.null(data)) {
    col <- if(!is.null(data.par$col)) data.par$col else par()$col
    pch <- if(!is.null(data.par$pch)) data.par$pch else par()$pch
    cex <- if(!is.null(data.par$cex)) data.par$cex else par()$cex
    points(data, col=col, pch=pch, cex=cex)
    if(!is.null(id.i <- data.par$id.i)) 
      text(data[id.i,1], data[id.i,2], id.i, cex=cex/1.5, pos=1)
    }
  d.levels <- exp(log.levels)  
  names(d.levels) <- as.character(probs)
  contour(x1, x2, pdf, levels=d.levels, 
    labels=paste("p=", as.character(probs), sep=""), add=TRUE, ...)
  if(landmarks != "") {
    if(landmarks == "auto") { 
      mean.type <-  "proper"  
      if(lc.family == "sc") mean.type <- "pseudo"      
      if(lc.family == "st") { if(dp[[4]] <= 1)  mean.type <- "pseudo"}
     }
    else 
    mean.type <- landmarks
    landmarks.label <- 
       c("origin", "mode", if(mean.type == "proper")  "mean" else "mean~")
    cp <- dp2cpMv(dp, family, cp.type=mean.type, upto=1)
    mode <- modeSECdistrMv(dp, family)
    x.pts <- c(dp$xi[1], mode[1], cp[[1]][1])
    y.pts <- c(dp$xi[2], mode[2], cp[[1]][2])
    points(x.pts, y.pts, ...)
    text(x.pts, y.pts, landmarks.label, pos=2, offset=0.3, ...)
    lines(x.pts, y.pts, lty=2)
    }  
  options(oo) 
  cL <- contourLines(x1, x2, pdf, levels=d.levels)
  for(j in 1:length(probs)) cL[[j]]$prob <- probs[j]
  return(list(x=x1, y=x2, names=compNames, density=pdf, contourLines=cL))
}    

plot.selm <- function(x, param.type="CP", which = c(1:4), caption, 
    panel = if (add.smooth) panel.smooth else points, main = "", 
    # sub.caption = NULL, 
    ask = prod(par("mfcol")) < length(which) && dev.interactive(), ..., 
    id.n = 3, labels.id = names(x@residuals.dp), 
    cex.id = 0.75, identline = TRUE, add.smooth = getOption("add.smooth"), 
    label.pos = c(4, 2), cex.caption = 1) 
{
    if(class(x) != "selm") stop("object not of class 'selm'")
    show <- rep(FALSE, 4)
    show[which] <- TRUE
    dots <- list(...)
    nmdots <- names(dots)  
    p <- slot(x, "size")["p"]
    if(missing(caption))  { caption <-  if(p> 1) 
      c("Residuals vs Fitted Values", 
       "Residual values and fitted error distribution", 
       "Q-Q plot of (scaled DP residuals)^2",
       "P-P plot of (scaled DP residuals)^2") else
      c("Boxplot of observed values", 
       "Empirical values and fitted distribution", 
       "Q-Q plot of (scaled DP residuals)^2",
       "P-P plot of (scaled DP residuals)^2")}
    all.par <- slot(x, "param")
    param.type <- tolower(param.type)  
    param <- all.par[[param.type]]
    if(is.null(param)) { message(paste(
        "Requested param.type='", param.type, "' evaluates to NULL.", sep=""))
      if(param.type == "pseudo-cp" & x@family== "SN") 
        message("Pseudo-CP makes no sense for SN family")
      if(param.type == "cp" & x@family== "SC")
        message("CP makes no sense for SC family")
      if(param.type == "cp" & x@family== "ST")
        message("CP of ST family requires nu>4")  
      stop("Consider another choice of param.type (DP or pseudo-CP)")
      }
    r <- residuals(x, param.type)
    r.lab <- paste(toupper(param.type), "residuals")
    dp <- if(length(all.par$fixed) > 0) all.par$dp.complete else all.par$dp
    nu. <- switch(x@family, ST = dp[p+3], SN = Inf, SC=1)  
    rs <- slot(x,"residuals.dp")/dp[p+1]
    rs2 <- rs^2
    n <- slot(x, "size")["n.obs"]
    yh <- fitted(x, param.type)    
    w <- weights(x)
    if (!is.null(w)) {
        wind <- (w != 0)
        r <- r[wind]
        yh <- yh[wind]
        w <- w[wind]
        labels.id <- labels.id[wind]
    }
    else w <- rep(1,n)
    rw <- n*w/slot(x,"size")["nw.obs"]
    cex.pts <- rw * if("cex" %in% nmdots) dots$cex else par("cex")
    if (is.null(id.n)) 
        id.n <- 0
    else {
        id.n <- as.integer(id.n)
        if (id.n < 0 || id.n > n) 
            stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)
    }
    if (id.n > 0) {
        if (is.null(labels.id)) 
            labels.id <- paste(1:n)
        iid <- 1:id.n
        # show.r <- sort.list(abs(r), decreasing = TRUE)[iid]        
        show.rs <- sort.list(rs2, decreasing = TRUE)[iid]
        # rs2.lab <- paste("(scaled DP residuals)^2")
        text.id <- function(x, y, ind, adj.x = TRUE) {
            labpos <- if (adj.x) 
                label.pos[1 + as.numeric(x > mean(range(x)))]
            else 3
            text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE, 
                pos = labpos, offset = 0.25)
        }
    }
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if (show[1]) {
      if(all(is.na(r)) & p>1)  message(paste("CP residuals not available;",
         "consider param.type='DP' or 'pseudo-CP'"))
      else { 
        if(p == 1){ 
          y <-  (x@residuals.dp + x@fitted.values.dp) 
          boxplot(y, plot=TRUE, col="gray85", border="gray60")
          }
        else { # p>1
        # if (id.n > 0) 
        #    ylim <- extendrange(r = ylim, f = 0.08)        
        ylim <- range(r, na.rm = TRUE)
        plot(yh, r, xlab = "Fitted values", ylab = r.lab, main = main, 
            ylim = ylim, type = "n")
        panel(yh, r, ...)  # previously it included 'cex=cex.pts'
        # if (one.fig) title(sub = sub.caption, ...)
        if (id.n > 0) {
          y.id <- r[show.rs]
          y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
          text.id(yh[show.rs], y.id, show.rs)
          }
        abline(h = 0, lty = 2, col = "gray")
        } }
    mtext(caption[1], 3, 0.5, cex = cex.caption) }
    if (show[2]) {
      if(all(is.na(r)) & p>1) message(
        "CP residuals not available; consider param.type='DP' or 'pseudo-CP'")
      else {
        if (p == 1){
          y <-  (x@residuals.dp + x@fitted.values.dp) 
          dp0 <- dp
          xlab="observed variable"}
        else {
          y <- r
          dp0 <- as.numeric(c(dp[1]-param[1], dp[-(1:p)]))
          xlab=r.lab
        }
        h <- hist(rep(y, w), plot=FALSE)
        extr <- extendrange(x=h$breaks)
        x.pts <- seq(max(extr), min(extr), length=501) 
        d.fn <- get(paste("d", tolower(x@family), sep=""), inherits = TRUE)
        pdf <- d.fn(x.pts, dp=dp0)
        plot(c(h$mids, x.pts), c(h$density, pdf), type="n", main=main, 
          xlab=xlab,  ylab="probability density")
        hist(rep(y, w), col="gray95", border="gray60", probability=TRUE, 
          freq=FALSE, add=TRUE)
        lines(x.pts, pdf, ...)
        rug(y, ticksize=0.02, ...)
        # if (id.n > 0) {     rug(y, ticksize=0.015, ...)
        #   text(y[show.rs], 0, labels.id[show.rs], srt=90, cex=0.5, pos=1, 
        #   offset=0.2) } 
        mtext(caption[2], 3, 0.25, cex = cex.caption)
      }}
    if (show[3]) {
      ylim <- c(0, max(pretty(rs2)))
      q <- qf((1:n)/(n+1), 1, nu.)
      plot(q, sort(rs2), xlab="Theoretical values", ylab="Empirical values", 
        ylim=ylim, type="p", main=main, ...)   # cex=cex.pts
      if(identline) abline(0, 1, lty = 2, col = "gray50")
      # if (one.fig) title(sub = sub.caption, ...)
      mtext(caption[3], 3, 0.25, cex = cex.caption)
      if (id.n > 0) text.id(q[n+1-iid], rs2[show.rs], show.rs) 
    }
    if (show[4]) {
      p <- (1:n)/(n+1)
      pr <- pf(sort(rs2), 1, nu.)
      plot(p, pr, xlab="Theoretical values", ylab="Empirical values",
         xlim=c(0,1), ylim=c(0,1), main=main, ...) # cex=cex.pts,
      if(identline) abline(0, 1, lty = 2, col = "gray50")
      # if (one.fig)  title(sub = sub.caption, ...)
      mtext(caption[4], 3, 0.25, cex = cex.caption)
      if(identline) abline(0, 1, lty = 2, col = "gray50")
      if (id.n > 0)  text.id(p[n+1-iid], pr[n+1-iid], show.rs)
    } 
    # if (!one.fig && par("oma")[3] >= 1) 
    #     mtext(sub.caption, outer = TRUE, cex = 1.25)
    invisible()
  }


print.summary.selm <- function(object)
{
    obj <- object
    digits = max(3, getOption("digits") - 3)
    cat("Call: ")
    print(slot(obj, "call"))
    n <- obj@size["n.obs"]
    cat("Number of observations:", n, "\n")
    if(!is.null(slot(obj,"aux")$weights))
      cat("Weighted number of observations:", obj@size["nw.obs"], "\n")
    show.family <- slot(obj,"family")
    cat("Family:", show.family,"\n")
    fixed <- slot(obj, "param.fixed") 
    if(length(fixed) > 0) { fixed.char <-
         paste(names(fixed), format(fixed), sep=" = ", collapse=", ")
         cat("Fixed parameters:", fixed.char, "\n") }
    method <- slot(obj, "method")
    u <- if(length(method)==1) NULL else paste(", penalty function:", method[2])
    cat("Estimation method: ", method[1], u, "\n", sep="")
    logL.name <- paste(if(method[1] == "MLE") "Log" else "Penalized log",
         "likelihood:", sep="-")
    cat(logL.name, format(slot(obj,"logL"), nsmall=2), "\n")
    param.type <- slot(obj, "param.type")
    cat("Parameter type:", param.type,"\n") 
    if((note <- slot(object,"note")) != "") cat(paste("Note:", note, "\n"))
    if(obj@boundary) 
      cat("Estimates on/near the boundary of the parameter space\n")
    resid <- slot(obj, "resid")
    if(n > 5) {
      nam <- c("Min", "1Q", "Median", "3Q", "Max")
      rq <- if (length(dim(resid)) == 2) 
            structure(apply(t(resid), 1, quantile), dimnames = list(nam, 
                dimnames(resid)[[2]]))
      else structure(quantile(resid), names = nam)
      cat("\n", param.type, " residuals:\n", sep="")
      print(rq, digits = digits)
    }
    param <- slot(obj,"param.table")
    p <- obj@size["p"]
    cat("\nRegression coefficients\n")
    printCoefmat(param[1:p, ,drop=FALSE], digits = digits,  
      signif.stars = getOption("show.signif.stars"), na.print = "NA")
    cat("\nParameters of the SEC random component\n")
    printCoefmat(param[(p+1):nrow(param), 1:2, drop=FALSE], digits = digits,  
      signif.stars = FALSE, na.print = "NA")  
    if(!is.null(obj@aux$param.cor)) {
      cat("\nCorrelations of parameter estimates:\n")
      print(obj@aux$param.cor)
      }
    if(!is.null(obj@aux$param.cov)) {
      cat("\nCovariances of parameter estimates:\n")
      print(obj@aux$param.cov)
      } 
  invisible(object)
  }


plot.mselm <- function (x, param.type="CP", which, caption, 
    panel = if (add.smooth) panel.smooth else points, main = "", 
    # sub.caption = NULL, 
    ask = prod(par("mfcol")) < length(which) && dev.interactive(), ..., 
    id.n = 3, labels.id = names(x@residuals.dp), 
    cex.id = 0.75, identline = TRUE, add.smooth = getOption("add.smooth"), 
    label.pos = c(4, 2), cex.caption = 1) 
  { 
    p <- slot(x,"size")["p"]
    if(missing(which)) which <- if(p == 1) c(1,3,4) else 2:4
    show <- rep(FALSE, 4)
    show[which] <- TRUE
    if(missing(caption)) caption <- 
       c("Observed values and fitted distribution", 
       paste("Distribution of", param.type, "residual values"),
       "Q-Q plot of Mahalanobis distances",
       "P-P plot of Mahalanobis distances")
    if(!show[2]) param.type <- "DP"   # CP-residuals only used for show[2]
    param.type <- tolower(param.type)  
    all.par <- slot(x, "param")
    param <- all.par[[param.type]]
    dots <- list(...)
    if(is.null(param)) { message(paste(
        "Requested param.type='", param.type, "' evaluates to NULL.", sep=""))
      if(param.type == "pseudo-cp" & x@family== "SN") 
        message("Pseudo-CP makes no sense for SN family")
      if(param.type == "cp" & x@family== "SC")
        message("CP makes no sense for SC family")
      if(param.type == "cp" & x@family== "ST")
        message("CP of ST family requires nu>4")  
      stop("Consider another choice of param.type")
      }
    r <- residuals(x, param.type)
    r.lab <- paste(toupper(param.type), "residuals")
    # family <- x@family 
    dp <- if(length(all.par$fixed) > 0) all.par$dp.complete else all.par$dp
    cp <- dp2cpMv(dp, family=x@family, cp.type="auto") 
    nu. <- switch(x@family, ST = dp$nu, SN = Inf, SC=1)     
    n <- slot(x,"size")["n.obs"]
    d  <- x@size["d"]
    yh <- fitted(x, param.type)    
    w <- weights(x)
    if (!is.null(w)) {
        wind <- w != 0
        r <- r[wind]
        yh <- yh[wind]
        w <- w[wind]
        labels.id <- labels.id[wind]
    }
    else w <- rep(1,n)
    rw <- n*w/slot(x,"size")["nw.obs"]
    if (is.null(id.n)) 
        id.n <- 0
    else {
        id.n <- as.integer(id.n)
        if (id.n < 0 || id.n > n) 
            stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)
    }
    Omega.inv <- pd.solve(dp$Omega, silent=TRUE)
    r.dp <- t(slot(x, "residuals.dp"))
    rs2 <- colSums((Omega.inv %*% r.dp) * r.dp)
    if (id.n > 0) {
		if (is.null(labels.id)) labels.id <- paste(1:n)
		iid <- 1:id.n
		show.r <- sort.list(abs(r), decreasing = TRUE)[iid]		
		show.rs <- sort.list(rs2, decreasing = TRUE)[iid]
		text.id <- function(x, y, ind, adj.x = TRUE) {
			labpos <- if (adj.x) 
				label.pos[1 + as.numeric(x > mean(range(x)))]
			else 3
			text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE, 
				pos = labpos, offset = 0.25)
		} 
      } else show.rs <- NULL
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if (show[1]) { # data scatter matrix and fitted curves (only if p=1)
      if(p == 1) {
        y <- (x@residuals.dp + x@fitted.values.dp)  
        fitted.distr <- makeSECdistr(dp, family=x@family, 
          name="fitted distribution", compNames=colnames(x@param$dp[[1]]))
        data.par <- list(col=dots$col, pch=dots$pch, cex=dots$cex,
          id.i=show.rs)
        plot(fitted.distr, landmarks="", data=y, main=main, data.par=data.par,
             ...)   # previously it included cex=sqrt(rw)
        # text.id(..) se d=1, ma se d>1 si deve fare per ogni pannello (?!)
        mtext(caption[1], 3, 1.5, cex = cex.caption)
        } else  
      message(paste("plot of (observed data, fitted distribution)",
              "makes no sense if covariates 'x' exist", 
              "and fitted distribution varies with 'x'"))
      }
    if (show[2]) { # scatter matrix of residuals and fitted curves
      dp0 <- dp
      dp0[[1]] <-  as.numeric((dp[[1]]-param[[1]])[1,])
      data.par <- list(col=dots$col, pch=dots$pch, cex=dots$cex,
            id.i=show.rs)
      resid.distr <- makeSECdistr(dp0, family=x@family, 
         name="Residual distribution", compNames=colnames(x@residuals.dp))
      plot(resid.distr, landmarks="", data=residuals(x, param.type), 
         main=main, data.par=data.par)
      # mtext(caption[2], 3, 0.25, cex = cex.caption)
      mtext(caption[2], 3, 1.5, cex = cex.caption)
      }
    if (show[3]) { # QQ-plot
      # ylim <- c(0, max(pretty(rs2)))
      q <- qf((1:n)/(n+1), d, nu.) * d
      plot(q, sort(rs2), xlab="theoretical values", ylab="empirical values",
           main=main, ...)  # cex=sqrt(rw) now dropped
      if(identline) abline(0, 1, lty = 2, col = "gray70")
      # if (one.fig) title(sub = sub.caption, ...)
      mtext(caption[3], 3, 0.25, cex = cex.caption)
      if (id.n > 0)  text.id(q[n+1-iid], rs2[show.rs], show.rs)
      }
    if (show[4]) { # PP-plot
      p <- pf(rs2/d, d, nu.)
      p0 <- (1:n)/(n+1) 
      plot(p0, sort(p),  xlab="theoretical values", ylab="empirical values",
         xlim=c(0,1), ylim=c(0,1), main=main, ...)  # cex=sqrt(rw) now dropped
      if(identline) abline(0, 1, lty = 2, col = "gray70")
      # if (one.fig) title(sub = sub.caption, ...)
      mtext(caption[4], 3, 0.25, cex = cex.caption)
      if (id.n > 0) text.id(p[show.rs], p0[n+1-iid], show.rs)
      } 
    # if (!one.fig && par("oma")[3] >= 1) 
    #    mtext(sub.caption, outer = TRUE, cex = 1.25)
    invisible()
  }


print.summary.mselm <-  function(object)
{
    obj <- object
    digits = max(3, getOption("digits") - 3)
    # cat("Obj: ", deparse(substitute(obj)),"\n")
    cat("Call: ")
    print(slot(obj,"call"))
    n <- obj@size["n.obs"]
    d <- obj@size["d"]
    # p <- obj@size["p"]
    cat("Number of observations:", n, "\n")
    nw <- obj@size["nw.obs"]
    if(n != nw)  cat("Weighted number of observations:", nw, "\n")
    family <- slot(obj, "family")
    cat("Family:", family, "\n")
    method <- slot(object, "method") 
    u <- if(length(method)==1) NULL else 
         paste(", penalty function:", method[2])
    cat("Estimation method: ", method[1], u, "\n", sep="")
    fixed <- slot(obj, "param.fixed")
    if(length(fixed) > 0) {fixed.char <- 
        paste(names(fixed), format(fixed), sep=" = ", collapse=", ")
      cat("Fixed parameters:", fixed.char, "\n") }
    cat("Log-likelihood:", format(slot(obj,"logL"), nsmall=2), "\n")
    cat("Parameter type:", obj@param.type,"\n") 
    if(obj@boundary) 
      cat("Estimates on/near the boundary of the parameter space\n")
    names <- dimnames(obj@scatter$matrix)[[1]]
    for(j in 1:d) {
      param <- obj@coef.tables[[j]]
      cat("\n--- Response variable No.", j, ": ", names[j],"\n",sep="")
      resid <- obj@resid[,j]
      if(n>5) {
        nam <- c("Min", "1Q", "Median", "3Q", "Max")
        rq <- if (length(dim(resid)) == 2) 
              structure(apply(t(resid), 1, quantile), dimnames = list(nam, 
                  dimnames(resid)[[2]]))
        else structure(quantile(resid), names = nam)
        cat(obj@param.type, "residuals\n")
        print(rq, digits = digits)
      }
      cat("\nRegression coefficients\n")
      printCoefmat(param[, ,drop=FALSE], digits = digits,  
        signif.stars = getOption("show.signif.stars"), na.print = "NA")
    }
    cat("\n--- Parameters of the SEC random component\n")
    cat("Scatter matrix: ", obj@scatter$name,"\n", sep="")
    print(obj@scatter$matrix) 
    if(length(obj@slant) > 0) {
      cat("\nSlant parameter: ", obj@slant$name, "\n", sep="")
      print(cbind(estimate=obj@slant$param, std.err=obj@slant$se))
      }
    if(length(obj@tail) > 0) {
       cat("\nTail-weight parameter: ", obj@tail$name, "\n", sep="")
       print(c(estimate=obj@tail$param, std.err=obj@tail$se))
      }
    if(!is.null(obj@aux$param.cor)) {
      cat("\nCorrelations of parameter estimates:\n")
      print(obj@aux$param.cor)
      }
    if(!is.null(obj@aux$param.cov)) {
      cat("\nVar-covariance matrix of parameter estimates:\n")
      print(obj@aux$param.cov)
      }
}

dp2op <- function(dp, family)
{
  nt <- switch(tolower(family), "sn" = 3, "esn" = 4, "st" = 4, "sc" = 3, NULL)
  if(is.null(nt)) stop("unknown family")
  op <- dp
  if (is.list(dp))  { # multivariate case 
    if(length(dp) != nt) stop("wrong length of 'dp'")
    Omega <- dp[[2]] 
    alpha <- dp[[3]]
    d <- length(alpha)
    tmp <- delta.etc(alpha, Omega)
    delta <- tmp$delta
    Omega.cor <- tmp$Omega.cor
    D.delta <- sqrt(1 - delta^2)    # (5.18) of SN book, but as vector
    lambda <- delta/D.delta         # (5.20)
    omega <- sqrt(diag(as.matrix(Omega)))
    Psi <- Omega - outer(omega*delta, omega*delta)  # four lines before (5.30)
    op[[2]] <- Psi
    op[[3]] <- lambda
    names(op)[2:3] <- c("Psi", "lambda")
    }
  else {  # univariate case
    p <- length(dp) - nt + 1
    if(p < 1) stop("wrong length of 'dp'")
    delta <- delta.etc(dp[p+2])
    op[p+1] <- dp[p+1] * sqrt(1 - delta^2)
    names(op)[(p+1):(p+2)] <- c("psi", "lambda")
    } 
  op
}

op2dp <- function(op, family) 
{
  nt <- switch(tolower(family), "sn" = 3, "esn" = 4, "st" = 4, "sc" = 3, NULL)
  if(is.null(nt)) stop("unknown family")
  dp <- op
  if(is.list(op)) { # multivariate case 
    if(length(op) != nt) stop("wrong length of 'op'")
    Psi <- op[[2]] 
    psi <- sqrt(diag(Psi))
    lambda <- op[[3]]
    delta <- lambda/sqrt(1 + lambda^2)
    D.delta <- sqrt(1 - delta^2)
    Psi.bar <- cov2cor(Psi)
    omega <- psi/D.delta
    tmp <- as.vector(pd.solve(Psi.bar) %*% lambda)
    dp[[2]] <- Psi + outer(psi*lambda, psi*lambda)  # four lines before (5.30)
    dp[[3]] <- (tmp/D.delta)/sqrt(1 + sum(lambda*tmp))  # (5.22)
    names(dp)[2:3] <- c("Omega", "alpha")
    } 
  else { # univariate case
    p <- length(op) - nt + 1
    if(p < 1) stop("wrong length of 'dp'")
    delta <- delta.etc(dp[p+2])
    dp[p+1] <- op[p+1]/sqrt(1 - delta^2)
    names(dp)[(p+1):(p+2)] <- c("omega", "alpha")
    }  
  dp
}
 
coef.selm <- function(object, param.type="CP", ...) {
    param <- slot(object,"param")[[tolower(param.type)]]
    if(is.null(param) & tolower(param.type)=="cp") {
        message("CP not defined, consider param.type='DP' or 'pseudo-CP'")
        return(NULL)}
    param} 
 
coef.mselm <- function(object, param.type="CP", vector=TRUE, ...) 
{
    list <- slot(object,"param")[[tolower(param.type)]]
    if(is.null(list) & tolower(param.type)=="cp") {
        message("CP not defined, consider param.type='DP' or 'pseudo-CP'")
        return(NULL)}
    if(!vector) return(list)
    as.vector(c(list[[1]], vech(list[[2]]), unlist(list[3:length(list)])))
}

extractSECdistr <- function(object, name, compNames) 
{
  obj.class <- class(object)
  if(!(obj.class %in% c("selm", "mselm")))
    stop(gettextf("wrong object class: '%s'", obj.class), domain = NA)
  param <- slot(object, "param")   
  dp <- if(length(param$dp.complete) > 0) param$dp.complete else param$dp
  p <- slot(object, "size")[2]
  if(obj.class == "selm")  {
    lead <- if(p > 1) 0 else dp[1]
    dp0 <- c(lead, dp[-(1:p)])  
    names(dp0)[1] <- "xi"
    }
  else { # class = "mselm"
    dp0 <- dp
    names(dp0)[1] <- "xi"
    dp0[[1]] <- if(p == 1) as.vector(dp0[[1]]) else 
                rep(0, slot(object, "size")[1])
    }
  if((obj.class == "mselm") & missing(compNames)) compNames <- names(dp$alpha)  
  if(missing(name)) {
    name <- paste("SEC distribution of", deparse(substitute(object))) 
    name <- if(p > 1) paste("Residual", name) else paste("Fitted", name)
    }   
  if(obj.class == "selm")  
    new("SECdistrUv", dp=dp0, family=slot(object, "family"), name=name) else 
    new("SECdistrMv", dp=dp0, family=slot(object, "family"), name=name,
        compNames=compNames)   
}


# introduce sd generic function, in the same fashion of package circular
#
sd <- function(x, ...) UseMethod("sd")
sd.default <- function(x, na.rm = FALSE, ...) stats::sd(x=x, na.rm=na.rm)

mean.SECdistrUv <- function(x) dp2cp(object=x, upto=1)
mean.SECdistrMv <- function(x) dp2cp(object=x, upto=1)[[1]]
sd.SECdistrUv <- function(x) dp2cp(object=x, upto=2)[2]
vcov.SECdistrMv <- function(object) dp2cp(object=object, upto=2)[[2]]

#---
#
profile.selm <- function(fitted, param.type, param.name, param.values, npt,
  opt.control=list(), plot.it=TRUE, log=TRUE, level, trace=FALSE, ...)
{ obj <- fitted
  obj.class <- class(obj)
  if(obj.class != "selm" | attr(obj.class, "package") != "sn") 
    stop(gettextf("wrong object class: '%s'", obj.class), domain = NA)
  param.type <- match.arg(toupper(param.type), c("DP", "CP"))
  family <- slot(obj, "family")
  obj.par <- slot(obj, "param")
  dp.full <- if(length(obj.par$fixed)==0) obj.par$dp else obj.par$dp.complete
  if(param.type == "CP") { 
    cp.full <- mle.full <- dp2cpUv(dp.full, family)
    profile.comp <- match(param.name, names(cp.full))
    }
  else {
    mle.full <- dp.full
    profile.comp <- match(param.name, names(dp.full))
    }
  fixed.names <- setdiff(names(obj.par$dp.complete), names(obj.par$dp))
  if(length(fixed.names) > 0) {
    fixed.comp <- match(fixed.names, names(dp.full))
    fixed.values <- mle.full[fixed.comp] 
    }
    else fixed.comp <- fixed.values <- NULL
  clash <- intersect(fixed.comp, profile.comp)
  if(length(clash) > 0)  stop(paste("parameter component No.", clash,
       "is fixed in the model, it cannot be profiled"))
  p <- slot(obj, "size")["p"]
  method <- slot(obj, "method")
  penalty <- if(method[1] == "MPLE") method[2] else NULL
  constr.comp <- c(profile.comp, fixed.comp)
  free.comp <- setdiff(1:length(dp.full), constr.comp)
  if(anyNA(profile.comp)) stop("some wrong item in param.name")
  npc <- length(profile.comp) # number of terms in profile.comp (either 1 or 2)
  if((npc != 1) && (npc != 2)) stop("wrong length(param.name)")
  if(missing(npt)) npt <- rep((50+npc) %/% npc, npc) else
     if(length(npt) != npc) npt <- rep(npt[1], npc)
  log.comp <- if(!log)  rep(NA, npc) else { 
    if(param.type == "DP") match(c("omega", "nu"), param.name, NULL)   
    else match(c("s.d.", "gamma2"), param.name, NULL) }  
  logScale <- (1:2) %in% which(!is.na(log.comp))   
  m <- slot(obj,"input")$model
  x <- model.matrix(attr(m, "terms"), data=m)
  w <- slot(obj, "input")$model$"(weights)"
  weights <- if(is.null(w)) rep(1, nrow(x)) else w
  opt.control$fnscale <- (-1)
  par.val <- param.values
  if(npc == 1) { # one-parameter profile logLik
    par.val <- as.vector(par.val)
    if(any(diff(par.val) <= 0)) stop("param.values not an increasing sequence")
    if(prod(range(par.val) - mle.full[profile.comp]) > 0)
      stop(gettextf("param range does not bracket the MLE/MPLE point: '%s'"),
        format(mle.full[profile.comp]), domain=NA)
    logScale <- logScale[1]
    if(length(par.val) > 2) npt <- length(par.val) else 
      par.val <- seqLog(par.val[1], par.val[2], length=npt, logScale)
    logL <- numeric(npt)
    for(k in 1:npt) {
      constr.values <- c(par.val[k], fixed.values)
      free.values <- mle.full[-constr.comp] 
      opt <- optim(free.values, constrained.logLik,  method="BFGS",
        control=opt.control, param.type=param.type, x=x, y=m[[1]],
        weights=weights, family=family, constr.comp=constr.comp, 
        constr.values=constr.values, penalty=penalty, trace=trace)
      logL[k] <- opt$value  
      }
    out <- list(call=match.call(), param=par.val, logLik=logL)
    names(out)[2] <- param.name  
    deviance <- 2*(logLik(obj) - logL)
    if(any(deviance + sqrt(.Machine$double.eps) < 0)) warning(paste(
      "A relative maximum of the (penalized) likelihood seems to have been",
      "taken as\n the MLE (or MPLE).",
      "Re-fit the model with starting values suggested by the plot."))
    s <- diff((sign(diff(deviance))))
    if(length(which(s != 0)) > 1) {
       message(paste("The log-likelihood function appears to have multiple",
        "maxima.\n", "Confidence intervals may be handled improperly.\n"))
       # readline("Press <Enter> to continue<cr>")
       # browser()
       }
    if(missing(level)) level <- 0.95
    level <- level[1]
    if(is.na(level) | level <= 0 | level >= 1) {
      message("illegal level value is reset to default value")
      level <- 0.95 }
    if(obj.par$boundary) {
      message("parameter estimates at the boundary, no confidence interval")
      level <- NULL
      } 
    if(!is.null(level)) {
      q <- qchisq(level[1], 1)
      if(deviance[1] < q | deviance[npt] < q) warning(
        "parameter range seems short; confidence interval may be inaccurate")
      dev.fn <- splinefun(par.val, deviance - q, method="monoH.FC")
      rootL <- try(uniroot(dev.fn, lower=min(par.val),  check.conv=TRUE, 
                   upper=mle.full[profile.comp],  extendInt="downX"))
      rootH <- try(uniroot(dev.fn, lower=mle.full[profile.comp], 
                   upper=max(par.val), check.conv=TRUE, extendInt="upX")) 
      fail.confint <- (class(rootL)=="try-error" | class(rootH)=="try-error")                   
      out$confint <- if(fail.confint) rep(NULL,2) else c(rootL$root, rootH$root)   
      out$level <- level                         
      }
    if(plot.it) {  
      if(logScale) { 
        par.val <-  log(par.val)
        param.name <- paste("log(", param.name, ")", sep="")
        }
      plot(par.val, deviance, type="l", xlab=param.name,
          ylab="2*{max(logLik) - logLik}", ...)
      if(logScale) {
          rug(log(mle.full[profile.comp]), ticksize = 0.02)
          if(is.null(level) | fail.confint) low <- hi <- NULL else { 
            low <- log(rootL$root)
            hi <- log(rootH$root) }}
        else {
          rug(mle.full[profile.comp], ticksize = 0.02)
          if(is.null(level)| fail.confint) low <- hi <- NULL else { 
            low <- rootL$root
            hi <- rootH$root
          }}
      if(!is.null(level) & !fail.confint) { 
        abline(h=q, lty=3, ...)
        lines(rep(low, 2), c(par()$usr[3], q), lty=3, ...)
        lines(rep(hi, 2), c(par()$usr[3], q), lty=3, ...)
        }
      }
    }
  else { # npc==2, two-parameter profile logLik
    if(length(par.val) != 2) stop("wrong dimension of param.values")
    u <- unlist(lapply(par.val, length))
    param1 <- par.val[[1]]
    param2 <- par.val[[2]]
    if(prod(range(param1) - mle.full[profile.comp][1]) > 0 |
          prod(range(param2) - mle.full[profile.comp][2]) > 0) stop(
       gettextf("parameter range does not bracket the MLE/MPLE point: '%s'",
           paste(format(mle.full[profile.comp]), collapse=",")), domain=NA)  
    if(u[1] > 2) npt[1] <- u[1] else 
      param1 <- seqLog(param1[1], param1[2], length=npt[1], logScale[1])
    if(u[2] > 2) npt[2] <- u[2] else   
      param2 <- seqLog(param2[1], param2[2], length=npt[2], logScale[2]) 
    logL <- matrix(NA, npt[1], npt[2])
    if(any(diff(param1) <= 0)) 
       stop("param.values[[1]] not an increasing sequence")
    if(any(diff(param2) <= 0)) 
       stop("param.values[[2]] not an increasing sequence")
    for(k1 in 1:npt[1]) for(k2 in 1:npt[2]){
      constr.values <- c(param1[k1], param2[k2], fixed.values)
      free.values <- mle.full[-constr.comp] 
      opt <- optim(free.values, constrained.logLik,  method="BFGS",
        control=opt.control, param.type=param.type, x=x, y=m[[1]],
        weights=weights, family=family, constr.comp=constr.comp, 
        constr.values=constr.values, penalty=penalty, trace=trace)  
      logL[k1,k2] <- opt$value  
      }
    out <- list(call=match.call(), param1=param1, param2=param2, logLik=logL)
    names(out)[2:3] <- param.name  
    if(missing(level)) level <- c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
    if(anyNA(level) | any(level<=0) | any(level>=1)) {
      message("illegal level values; vector 'level' reset to default value")
      level <- c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99) }
    if(obj.par$boundary) {
      message("parameter estimates at the boundary, no confidence regions")
      level <- NULL
      } 
    q <- if(is.null(level)) c(1, 2, 5, 8, 15) else  qchisq(level, 2) 
    deviance <- 2*(logLik(obj)-logL)
    if(any(deviance + sqrt(.Machine$double.eps) < 0)) warning(paste(
      "A relative maximum of the (penalized) likelihood seems to have taken",
      "as\n the MLE (or MPLE).",
      "Re-fit the model with starting values suggested by the plot."))
    cL <- contourLines(param1, param2, deviance, levels=q)
    out$deviance.contour <- cL
    if(!is.null(level)) for(j in 1:length(cL)) {
      k <- which(q == cL[[j]]$level)
      out$deviance.contour[[j]]$prob <- level[k]
      }
    if(plot.it) {
      if(logScale[1]) { 
	    param1 <-  log(param1)
    	param.name[1] <- paste("log(", param.name[1], ")", sep="")
    	}		  
      if(logScale[2]) {
	    param2 <-  log(param2)
	    param.name[2] <- paste("log(", param.name[2], ")", sep="")
	    }      
      contour(param1, param2, deviance, levels=q, labels=level,
          xlab=param.name[1], ylab=param.name[2], ...)
      mark <- mle.full[profile.comp]
      mark[logScale] <- log(mark[logScale])
      points(mark[1], mark[2], ...)  
      }
    }
  invisible(out)
}

constrained.logLik <- function(free.param, param.type, x, y, weights, family, 
    constr.comp=NA, constr.values=NA, penalty=NULL, trace=FALSE)
{
 if(trace) cat("constrained.logLik, free.param:", format(free.param))
  n <- sum(weights)
  p <- ncol(x)
  param <- numeric(length(free.param) + length(constr.values))
  param[constr.comp] <- constr.values
  param[-constr.comp] <- free.param
  par0 <- c(0, param[-(1:p)])
  # if(par0[2] <= 0) return(-Inf)
  # if(family=="ST" & par0[4] <= 0) return(-Inf) 
  # if(family=="ST" & par0[4] > 1e4) par0[4] <- Inf
  dp0 <- if(param.type =="DP") par0 else 
     cp2dpUv(par0, family, tol=1e-7, silent=TRUE)
  if(anyNA(dp0)) {
      if(is.null(dp0)) {message("null dp0"); browser()}
    excess <- attr(dp0, "excess")
    if(length(excess) == 0) {message("0-length excess"); browser() }
    if(is.null(excess) | is.na(excess) | abs(excess)==Inf ) 
        excess <- (.Machine$double.xmax)^(1/3)
        # {message("bad excess"); browser()}
    return(-1e9 * (1+ excess)^2)
    } 
  d.fn <- get(paste("d", tolower(family), sep=""), inherits = TRUE)
  logL <- try(d.fn((y - x %*% param[1:p]), dp=dp0, log=TRUE))
  if(class(logL) == "try-error") browser()
  Q <- if(is.null(penalty)) 0 else {
    penalty.fn <- get(penalty, inherits = TRUE)
    nu <- if(family=="ST") par0[4] else NULL
    penalty.fn(dp0[3], nu)
    } 
  out <- if(anyNA(logL)) -Inf else sum(logL * weights) - Q 
  if(trace) cat(", logL:", format(out), "\n")
  return(out)
}

seqLog <- function(from, to, length, logScale=FALSE) {
  if(logScale & any(c(from, to) <= 0)) 
    stop("logScale requires positive arguments 'from' and 'to'")
  if(logScale)  exp(seq(log(from), log(to), length.out=length)) else
    seq(from, to, length.out=length) 
  }
