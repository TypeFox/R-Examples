## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


### List of supported Archimedean copulas

## NOTA BENE:  Write psi(), tau(), ... functions such that they *vectorize*
## ---------   *and* do work even for (non-numeric) NULL argument
##		{now checked in "acopula" validityMethod -- see ./AllClass.R }

### Ali-Mikhail-Haq, see Nelsen (2007) p. 116, # 3 #############################

##' AMH object
copAMH <-
    (function() { ## to get an environment where  C.  itself is accessible
	C. <- new("acopula", name = "AMH",
		  ## generator
		  psi = function(t, theta) { (1-theta)/(exp(t+0)-theta) },
		  iPsi = function(u, theta, log=FALSE) {
                      res <- log((1-theta*(1-u))/u) # alternative: log1p((1-theta)*(1/u-1))
                      if(log) log(res) else res
		  },
		  ## parameter interval
		  paraInterval = interval("[0,1)"),
		  ## absolute value of generator derivatives
		  absdPsi = function(t, theta, degree = 1, n.MC = 0, log = FALSE,
				     is.log.t = FALSE,
				     method = "negI-s-Eulerian", Li.log.arg=TRUE)
	      {
		  lth <- log(theta)
		  if(n.MC > 0) {
		      if(is.log.t) t <- exp(t) # very cheap for now
		      absdPsiMC(t, family="AMH", theta=theta, degree=degree,
				n.MC=n.MC, log=log)
		  } else {
                      ## FIXME: deal with  is.log.t
		      if(is.log.t) t <- exp(t) # very cheap for now

		      ## Note: absdPsi(0, ...) is correct, namely (1-theta)/theta * polylog(theta, s=-degree)
		      if(theta == 0) return(if(log) -t else exp(-t)) # independence
		      Li.arg <- if(Li.log.arg) lth - t else theta*exp(-t)
		      Li. <- polylog(Li.arg, s = -degree, method=method, is.log.z = Li.log.arg, log=log)
		      if(log)
			  Li. + log1p(-theta)-lth
		      else
			  Li. * (1-theta)/theta
		  }
	      },
		  ## derivatives of the generator inverse
		  absdiPsi = function(u, theta, degree=1, log=FALSE) {
		      switch(degree,
			     ## 1 :
			     if(log) log1p(-theta)-log(u)-log1p(theta*(u-1))
			     else (1-theta)/(u*(1-theta*(1-u))),
			     ## 2 :
			     if(log) log1p(-theta)+ log1p(1 + theta * (2*u - 1)) -2*(log(u)+log1p(theta*(u-1)))
			     else (1-theta) * (1 + theta * (2*u - 1)) / (u*(1 + theta * (u-1)))^2,
			     ## >= 3:
			     stop("not yet implemented for degree > 2"))
		  },
		  ## density of the diagonal
		  dDiag = function(u, theta, d, log=FALSE) {
                      x <- (1-(1-u)*theta)/u
		      if(any(iI <- is.infinite(x))) { ## for u == 0 (seems unneeded for small u << 1 ?)
			  u[iI] <- if(log) -Inf else 0
			  ok <- !iI
			  u[ok] <- C.@dDiag(u[ok], theta, d=d, log=log)
			  return(u)
		      }
		      if(log) log(d)+ (d-1)*log(x) + 2*(log(x-theta)-log(x^d-theta))
		      else d* x^(d-1) * ((x-theta)/(x^d-theta))^2
                  },
		  ## density  AMH
		  dacopula = function(u, theta, n.MC=0, log=FALSE, checkPar=TRUE,
				      method = "negI-s-Eulerian", Li.log.arg=TRUE)
              {
		  ## if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
                  if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                  n <- nrow(u)
		  if(d > 2 && !C.@paraConstr(theta)) {
		      if(checkPar) stop("theta is outside its defined interval")
		      return(rep.int(if(log) -Inf else 0, n))
		  }
                  ## f() := NaN outside and on the boundary of the unit hypercube
                  res <- rep.int(NaN, n)
		  n01 <- u.in.01(u)## indices for which density has to be evaluated
                  ## if(!any(n01)) return(res)
                  if(theta == 0) { res[n01] <- if(log) 0 else 1; return(res) } # independence
                  ## auxiliary results
                  u. <- u[n01,, drop=FALSE]
                  tIu <- -theta*(1-u.)
                  sum. <- rowSums(log(u.))
                  sumIu <- rowSums(log1p(tIu))
                  ## main part
                  if(n.MC > 0) { # Monte Carlo
                      V <- C.@V0(n.MC, theta)
                      res[n01] <- colMeans(exp(d*(log1p(-theta) + log(V)) +
                                               (V-1) %*% t(sum.) - (V+1) %*% t(sumIu)))
                      if(log) log(res) else res
                  } else { # explicit
                      Li.arg <-
                          if(Li.log.arg) log(theta) + sum. - sumIu else theta* apply(u./(1+tIu), 1, prod)
                      Li. <- polylog(Li.arg, s = -d, method=method, is.log.z = Li.log.arg, log=TRUE)
                      res[n01] <- (d+1)*log1p(-theta)-log(theta)+ Li. -sum. -sumIu
                      if(log) res else exp(res)
                  }
              },
		  ## score function
		  score = function(u, theta) {
		      if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
		      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
		      if(d > 2) stopifnot(C.@paraConstr(theta))
		      omu <- 1-u
		      b <- rowSums(omu/(1-theta*omu))
		      h <- theta*apply(u/(1-theta*omu), 1, prod)
		      -(d+1)/(1-theta) - 1/theta + b + (b+1/theta) *
			  polylog(h, s=-(d+1), method="negI-s-Stirling") /
			      polylog(h, s=-d, method="negI-s-Stirling")
		  },
                  ## uscore function
                  uscore = function(u, theta, d, method = "negI-s-Eulerian") {
                      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      omu <- 1-u
                      h <- theta*apply(u/(1-theta*omu), 1, prod)
                      Li.md1 <- polylog(h, s=-(d+1), method=method, log=TRUE)
                      Li.md <- polylog(h, s=-d, method=method, log=TRUE)
                      (1-theta) / (u*(1-theta*omu)) * (theta*exp(Li.md1-Li.md) - 1)
                  },
		  ## nesting constraint
		  nestConstr = function(theta0,theta1) {
		      C.@paraConstr(theta0) &&
		      C.@paraConstr(theta1) && theta1 >= theta0
		  },
		  ## V0 with density dV0 and V01 with density dV01 corresponding to
		  ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
		  V0 = function(n,theta, log=FALSE) {
                      if(log) stop("'log=TRUE' not yet implemented")
                      else
                          rgeom(n, 1-theta) + 1
                  },
		  dV0 = function(x,theta,log = FALSE) dgeom(x-1, 1-theta, log),
                  V01 = function(V0,theta0,theta1) {
                      rnbinom(length(V0),V0,(1-theta1)/(1-theta0))+V0
                  },
                  dV01 = function(x,V0,theta0,theta1,log = FALSE) {
                      stopifnot(length(V0) == 1 || length(x) == length(V0))
                      dnbinom(x-V0,V0,(1-theta1)/(1-theta0),log = log)
                  },
                  ## Kendall's tau
                  tau = tauAMH, ##-> ./aux-acopula.R
                  ## function(th)  1 - 2*((1-th)*(1-th)*log(1-th)+th)/(3*th*th)
                  ## but numerically stable, including, theta -> 0
                  iTau = function(tau, tol=.Machine$double.eps^0.25, check=TRUE, warn=TRUE, ...) {
		      if((L <- any(tau > 1/3)) || any(tau < 0)) {
			  ct <- if(length(ct <- sort(tau, decreasing = L)) > 3)
			      paste0(paste(format(ct[1:3]),collapse=", "),", ...") else format(ct)
			  msg <- "For AMH, Kendall's tau must be in [0, 1/3], but ("
			  if(check)
			      stop(msg, if(L) "largest" else "smallest", " sorted) tau = ", ct)
			  else {
                              if(warn)
                                  warning(msg, if(L) "largest" else "smallest", " sorted) tau = ", ct)
			      r <- tau
			      r[sml <- tau <=  0 ] <- 0
			      r[lrg <- tau >= 1/3] <- 1
			      ok <- !sml & !lrg
			      r[ok] <- C.@iTau(tau[ok], tol=tol, ...)
			      return(r)
			  }
		      }
		      vapply(tau, function(tau) {
			  if(abs(tau - 1/3) < 1e-10) return(1.) ## else:
			  r <- safeUroot(function(th) tauAMH(th) - tau,
					 interval = c(0, 1-1e-12),
					 Sig = +1, tol = tol, check.conv=TRUE, ...)
			  r$root
		      }, 0.)
		  },
                  ## lower tail dependence coefficient lambda_l
		  lambdaL = function(theta) { 0*theta },
		  lambdaLInv = function(lambda) {
		      if(any(lambda != 0))
			  stop("Any parameter for an Ali-Mikhail-Haq copula gives lambdaL = 0")
		      NA * lambda
		  },
		  ## upper tail dependence coefficient lambda_u
		  lambdaU = function(theta) { 0*theta },
		  lambdaUInv = function(lambda) {
		      if(any(lambda != 0))
			  stop("Any parameter for an Ali-Mikhail-Haq copula gives lambdaU = 0")
		      NA * lambda
		  }
		  )
	C.
    })()# {copAMH}


### Clayton, see Nelsen (2007) p. 116, #1 (slightly simpler form) ##############

##' Clayton object
copClayton <-
    (function() { ## to get an environment where  C.  itself is accessible
	C. <- new("acopula", name = "Clayton",
		  ## generator
		  psi = function(t, theta) { (1+t)^(-1/theta) },
		  iPsi = function(u, theta, log=FALSE) {
                      res <- u^(-theta) - 1
                      if(log) log(res) else res
		  },
		  ## parameter interval
		  paraInterval = interval("(0,Inf)"),
		  ## absolute value of generator derivatives
		  absdPsi = function(t, theta, degree=1, n.MC=0, log=FALSE) {
                      if(n.MC > 0) {
                          absdPsiMC(t, family="Clayton", theta=theta, degree=degree,
                                    n.MC=n.MC, log=log)
                      } else {
                          ## Note: absdPsi(0, ...) is correct, namely gamma(d+1/theta)/gamma(1/theta)
                          alpha <- 1/theta
                          res <- lgamma(degree+alpha)-lgamma(alpha)-(degree+alpha)*log1p(t)
                          if(log) res else exp(res)
                      }
                  },
                  ## derivatives of the generator inverse
		  absdiPsi = function(u, theta, degree=1, log=FALSE) {
		      switch(degree,
			     ## 1 :
			     if(log) log(theta)-(1+theta)*log(u) else theta*u^(-(1+theta)),
			     ## 2 :
			     if(log) log(theta)+log1p(theta)-(theta+2)*log(u)
			     else theta * (1 + theta) * u^-(theta + 2),
			     ## >= 3:
			     stop("not yet implemented for degree > 2"))
		  },
		  ## density of the diagonal
		  dDiag = function(u, theta, d, log=FALSE){
                      if(log) log(d)-(1+1/theta)*log(1+(d-1)*(1-u^theta)) else
                      d*(1+(d-1)*(1-u^theta))^(-(1+1/theta))
                  },
		  ## density  Clayton
		  dacopula = function(u, theta, n.MC=0, log=FALSE, checkPar=TRUE) {
		      ## if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
		      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
		      n <- nrow(u)
		      if(d > 2 && !C.@paraConstr(theta)) {
			  if(checkPar) stop("theta is outside its defined interval")
			  return(rep.int(if(log) -Inf else 0, n))
		      }
		      ## f() := NaN outside and on the boundary of the unit hypercube
		      res <- rep.int(NaN, n)
		      n01 <- u.in.01(u)## indices for which density has to be evaluated
		      ## if(!any(n01)) return(res)
		      ## auxiliary results
		      u. <- u[n01,, drop=FALSE]
		      lu <- rowSums(log(u.))
		      t <- rowSums(C.@iPsi(u., theta))
		      ## main part
		      if(n.MC > 0) { # Monte Carlo
			  lu.mat <- matrix(rep(lu, n.MC), nrow=n.MC, byrow=TRUE)
			  V <- C.@V0(n.MC, theta)
			  ## stably compute log(colMeans(exp(lx)))
			  lx <- d*(log(theta) + log(V)) - log(n.MC) - (1+theta)*lu.mat - V %*% t(t) # matrix of exponents; dimension n.MC x n ["V x u"]
			  ## note: smle goes wrong if:
			  ##       (1) d*log(theta*V) [old code]
			  ##       (2) U is small (simultaneously)
			  ##       (3) theta is large
			  res[n01] <- lsum(lx)
		      } else { # explicit
			  res[n01] <- sum(log1p(theta*(0:(d-1)))) - (1+theta)*lu -
                              (d+1/theta)*log1p(t)
		      }
		      if(log) res else exp(res)
		  },
		  ## score function
		  score = function(u, theta) {
		      if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
		      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
		      if(d > 2) stopifnot(C.@paraConstr(theta))
		      lu <- log(u)
		      t <- rowSums(C.@iPsi(u, theta=theta))
		      ltp1 <- log(1+t)
		      lx <- t(log(-lu)-theta*lu) # caution: lsum needs an (d,n)-matrix
		      ldt <- lsum(lx) # log of the derivative of t w.r.t. theta
		      alpha <- 1/theta
		      k <- 0:(d-1)
		      sum(k/(theta*k+1))-rowSums(lu)+alpha^2*ltp1-(d+alpha)*
                          exp(ldt-ltp1)
		  },
                  ## uscore function
                  uscore = function(u, theta, d) {
                      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      t <- rowSums(iPsi(u, theta=theta))
                      tht <- (theta*d+1)/(1+t)
                      t1 <- theta+1
                      tht*u^(-t1)-t1/u
                  },
		  ## nesting constraint
		  nestConstr = function(theta0,theta1) {
		      C.@paraConstr(theta0) &&
		      C.@paraConstr(theta1) && theta1 >= theta0
		  },
		  ## V0 with density dV0 and V01 with density dV01 corresponding to
		  ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
		  V0 = function(n,theta, log=FALSE) {
                      if(log) stop("'log=TRUE' not yet implemented")
                      else rgamma(n, shape = 1/theta) },
		  dV0 = function(x,theta,log = FALSE) dgamma(x, shape = 1/theta, log),
		  V01 = function(V0,theta0,theta1) { retstable(alpha=theta0/theta1, V0) },
		  dV01 = function(x,V0,theta0,theta1,log = FALSE) {
		      stopifnot(length(V0) == 1 || length(x) == length(V0))
		      alpha <- theta0/theta1
		      gamma <- (cospi2(alpha)*V0)^(1/alpha)
		      delta <- V0*(alpha == 1)
		      ## NB: new dstable() is vectorized in (x, gamma, delta) [but not the others]
		      dst <- dstable(x, alpha=alpha, beta = 1, gamma=gamma, delta=delta,
				     pm = 1, log=log, tol = 128* .Machine$double.eps)
		      if(log) V0-x + dst else exp(V0-x) * dst
		  },
		  ## Kendall's tau
		  tau = function(theta) { theta/(theta+2) },
		  iTau = function(tau) { 2*tau/(1-tau) },
		  ## lower tail dependence coefficient lambda_l
		  lambdaL = function(theta) { 2^(-1/theta) },
		  lambdaLInv = function(lambda) { -1/log2(lambda) },
		  ## upper tail dependence coefficient lambda_u
		  lambdaU = function(theta) { 0*theta },
		  lambdaUInv = function(lambda) {
		      if(any(lambda != 0))
			  stop("Any parameter for a Clayton copula gives lambdaU = 0")
		      NA * lambda
		  }
		  )
	C.
    })()# {copClayton}


### Frank, see Nelsen (2007) p. 116, # 5 #######################################

##' Frank object
copFrank <-
    (function() { ## to get an environment where  C.  itself is accessible
	C. <- new("acopula", name = "Frank",
		  ## generator
		  psi = function(t, theta) {
                      ## = -log(1-(1-exp(-theta))*exp(-t))/theta
		      ## -log1p(expm1(-theta)*exp(0-t))/theta # fails really small t, theta > 38
		      -log1mexp(t-log1mexp(theta))/theta
		  },
		  iPsi = function(u, theta, log=FALSE) {
		      ## == -log( (exp(-theta*u)-1) / (exp(-theta)-1) )
		      thu <- u*theta # (-> recycling args)
		      if(!length(thu)) return(thu) # {just for numeric(0) ..hmm}
		      et1 <- expm1(-theta) # e^{-th} - 1 < 0
		      ## FIXME ifelse() is not quite efficient
		      ## FIXME(2): the "> c* th" is pi*Handgelenk
### FIXME: use  delta = exp(-thu)*(1 - exp(thu-th)/ (-et1) =
###                   = exp(-thu)* expm1(thu-theta)/et1   (4)

###-- do small Rmpfr-study to find the best form -- (4) above and the three forms below

		      r <- ifelse(abs(thu) > .01*abs(theta), # thu = u*th > .01*th <==> u > 0.01
                              {   e.t <- exp(-theta)
                                  ifelse(e.t > 0 & abs(theta - thu) < 1/2,# th -th*u = th(1-u) < 1/2
                                         -log1p(e.t * expm1(theta - thu)/et1),
                                         -log1p((exp(-thu)- e.t) / et1))
                              },
                                  ## for small u (u < 0.01) :
                                  -log(expm1(-thu)/et1))
                      if(log) log(r) else r
		  },
		  ## parameter interval
		  paraInterval = interval("(0,Inf)"),
		  ## absolute value of generator derivatives
		  absdPsi = function(t, theta, degree = 1, n.MC = 0, log = FALSE, is.log.t = FALSE,
				     method = "negI-s-Eulerian", Li.log.arg = TRUE)
              {
                  if(n.MC > 0) {
                      absdPsiMC(t, family="Frank", theta=theta, degree=degree,
                                n.MC=n.MC, log=log)
                  } else {
                      ## Note: absdPsi(0, ...) is correct, namely (1/theta)*polylog(1-exp(-theta), s=-(degree-1))
		      Li.arg <- if(Li.log.arg) log1mexp(theta) - t else -expm1(-theta)*exp(-t)
		      Li. <- polylog(Li.arg, s = -(degree-1), log=log,
				     method=method, is.log.z = Li.log.arg)
		      if(log) Li. - log(theta) else Li. / theta
		  }
	      },
		  ## derivatives of the generator inverse
		  absdiPsi = function(u, theta, degree=1, log=FALSE) {
		      ut <- u*theta
		      switch(degree,
			     ## 1 :
			     if(log) log(theta) - (ut + log1mexp(ut))
			     else theta/expm1(ut),
			     ## 2 :
			     if(log) 2*log(theta) + ut - 2*log1mexp(-ut)
			     else (theta^2 * exp(ut))/expm1(ut)^2,
			     ## >= 3:
			     stop("not yet implemented for degree > 2"))
		  },
		  ## density of the diagonal
		  dDiag = function(u, theta, d, log=FALSE) {
		      r <- ut <- u*theta
		      ## h <- -expm1(-ut) # 1 - exp(-u * theta)
		      ## L <- ut > log(2)*.Machine$double.digits # <==> exp(-ut) < Mach..eps <==> h=1
                      ## empirically, this is "better:"
		      ## (more in ../demo/dDiag-plots.R -- should also depend on d !
		      L <- ut > 25   # [L]arge
		      S <- ut < 0.10 # [S]mall
		      M <- !L & !S   # [M]edium
		      ## recycle (u, theta):
		      u	 <- rep(u,     length.out=length(ut))
		      th <- rep(theta, length.out=length(ut))
		      r[L] <- dDiagFrank(u[L], th[L], d, log=log, method = "poly2")
		      r[M] <- dDiagFrank(u[M], th[M], d, log=log, method = "poly4")
		      r[S] <- dDiagFrank(u[S], th[S], d, log=log, method = "m1")
		      r
		  },
		  ## density  Frank
		  dacopula = function(u, theta, n.MC=0, log=FALSE, checkPar=TRUE,
				      method = "negI-s-Eulerian", Li.log.arg=TRUE)
	      {
		  ## if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
		  if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
		  n <- nrow(u)
		  if(d > 2 && !C.@paraConstr(theta)) {
		      if(checkPar) stop("theta is outside its defined interval")
		      return(rep.int(if(log) -Inf else 0, n))
		  }
                  ## f() := NaN outside and on the boundary of the unit hypercube
                  res <- rep.int(NaN, n)
		  n01 <- u.in.01(u)## indices for which density has to be evaluated
		  ## if(!any(n01)) return(res)
		  ## auxiliary results
		  u. <- u[n01,, drop=FALSE]
		  u.sum <- rowSums(u.)
		  lp  <- log1mexp(theta)    # log(1 - exp(-theta))
		  lpu <- log1mexp(theta*u.) # log(1 - exp(-theta * u))
		  lu <- rowSums(lpu)
		  ## main part
		  res[n01] <-
		      if(n.MC > 0) { # Monte Carlo
			  V <- C.@V0(n.MC, theta)
			  lx <- rep(-theta*u.sum, each=n.MC) + d*(log(theta) +
                                                  log(V) - V*lp) - log(n.MC) +
                                                      (V-1) %*% t(lu)
			  lsum(lx)
		      } else { # explicit
			  Li.arg <-
			      if(Li.log.arg) lp + rowSums(lpu-lp)
			      else -expm1(-theta)*exp(rowSums(lpu-lp))
			  Li. <- polylog(Li.arg, s = -(d-1), log=TRUE,
					 method=method, is.log.z = Li.log.arg)
			  (d-1)*log(theta) + Li. - theta*u.sum - lu
		      }
		  if(log) res else exp(res)
	      },
		  ## score function
		  score = function(u, theta) {
		      if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
		      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
		      if(d > 2) stopifnot(C.@paraConstr(theta))
		      e <- exp(-theta)
		      Ie <- -expm1(-theta) # == 1 - e == 1-e^{-theta}
		      etu <- exp(mtu <- -theta*u) # exp(-theta*u)
		      Ietu <- -expm1(mtu) ## = 1 - etu == 1 - exp(-theta*u)
		      ## FIXME: allow 'Li.log.arg' -> psilog(*, is.log.z) as for dacopula() above
		      h <- Ie*apply(Ie/etu, 1, prod)
		      factor <- rowSums(u*etu/Ietu) - (d-1)*e/Ie
		      (d-1)/theta - rowSums(u/Ietu) + factor *
			  polylog(h, s=-d, method="negI-s-Stirling") /
                              polylog(h, s=-(d-1), method="negI-s-Stirling")
		  },
                  ## uscore function
                  uscore = function(u, theta, d, method = "negI-s-Eulerian") {
                      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      Ie <- -expm1(-theta) # == 1 - e == 1-e^{-theta}
                      etu <- exp(mtu <- -theta*u) # exp(-theta*u)
                      h <- Ie*apply(Ie/etu, 1, prod)
                      factor <- theta/(-expm1(mtu))
                      Li.md <- polylog(h, s=-d, method=method, log=TRUE)
                      Li.mdm1 <- polylog(h, s=-(d-1), method=method, log=TRUE)
                      factor * (exp(Li.md+log(h)-theta*u - Li.mdm1) - 1)
                  },
		  ## nesting constraint
		  nestConstr = function(theta0,theta1) {
		      C.@paraConstr(theta0) &&
		      C.@paraConstr(theta1) && theta1 >= theta0
		  },
		  ## V0 (algorithm of Kemp (1981)) with density dV0 and V01 with density
		  ## dV01 corresponding to LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
		  V0 = function(n,theta, log=FALSE) {
                      if(log) stop("'log=TRUE' not yet implemented")
                      else {
                          if(-theta < .Machine$double.min.exp*log(2))
                              stop("'theta' is too large in Frank's V0();  maybe use log=TRUE")
                          rlog(n, -expm1(-theta), exp(-theta))
                      }
                  },
		  dV0 = function(x,theta, log = FALSE) {
		      if(any(x != (x. <- round(x)))) {
			  x <- x.; warning("x has been rounded to integer")
		      }
		      if(log) {
			  x*log1mexp(theta) - log(x*theta)
		      } else {
			  p <- -expm1(-theta) # == 1-exp(-theta)
			  p^x/(x*theta)
		      }
		  },
		  V01 = function(V0, theta0, theta1, rej = 1, approx = 10000) {
		      ## rej method switch: if V0*theta_0*p0^(V0-1) > rej a rejection
		      ## from F01 of Joe is applied (otherwise, the sum is
		      ## sampled via a logarithmic envelope for the summands)
		      ## approx is the largest number of summands before asymptotics is used (see copJoe@V01)
		      ## FIXME: optimal value of rej (for approx = 10000) is not clear yet; rej = 1 is not bad, however;
		      ##	lgammacor gives underflow warnings for rej < 1
		      rF01Frank(V0, theta0, theta1, rej, approx)
		  },
		  dV01 = function(x,V0,theta0,theta1,log = FALSE) {
		      stopifnot(length(V0) == 1 || length(x) == length(V0), all(x >= V0))
		      lfactor <- x*log1mexp(theta1) - V0*log1mexp(theta0)
		      res <- lfactor + dsumSibuya(x, V0, theta0/theta1, log=TRUE)
		      if(log) res else exp(res)
		  },
		  ## Kendall's tau; debye1() calls debye_1() from package 'gsl' :
		  tau = function(theta) {
		      if((l <- length(theta)) == 0) return(numeric(0)) # to work with NULL
		      res <- numeric(l)
		      res[isN <- theta == 0] <- 0 # limiting case
		      res[na <- is.na(theta)] <- NA
		      if(any(i <- !(na | isN)))
			  res[i] <- 1 + 4*(debye1(theta[i]) - 1)/theta[i]
		      res
		  },
		  iTau = function(tau, tol = .Machine$double.eps^0.25, ...) {
		      res <- tau
		      isN <- tau == 0 ## res[isN] <- 0 # limiting case
		      if(length(nn <- which(!isN)))
			  res[nn] <- vapply(tau[nn], function(tau) {
			      r <- safeUroot(function(th) C.@tau(th) - tau,
					     ## interval: from experiments on x86_64 Linux (2012-08)
					     interval = if(tau > 0) c(0, 7.21e+16) else c(-1.81e+16, 0),
					     Sig = +1, tol=tol,
					     check.conv=TRUE, ...)
			      r$root
			  }, 0.)
		      res
		  },
		  ## lower tail dependence coefficient lambda_l
		  lambdaL = function(theta) { 0*theta },
		  lambdaLInv = function(lambda) {
		      if(any(lambda != 0))
			  stop("Any parameter for a Frank copula gives lambdaL = 0")
		      NA * lambda
		  },
		  ## upper tail dependence coefficient lambda_u
		  lambdaU = function(theta) { 0*theta },
		  lambdaUInv = function(lambda) {
		      if(any(lambda != 0))
			  stop("Any parameter for a Frank copula gives lambdaU = 0")
		      NA * lambda
		  }
		  )
	C.
    })()# {copFrank}


### Gumbel, see Nelsen (2007) p. 116, # 4 ######################################

##' Gumbel object
copGumbel <-
    (function() { ## to get an environment where  C.  itself is accessible
	C. <- new("acopula", name = "Gumbel",
		  ## generator
		  psi = function(t, theta) { exp(-t^(1/theta)) },
		  iPsi = function(u, theta, log=FALSE) {
                      if(log) theta*log(-log(u)) else (-log(u+0))^theta
		  },
		  ## parameter interval
		  paraInterval = interval("[1,Inf)"),
		  ## absolute value of generator derivatives
		  absdPsi = function(t, theta, degree=1, n.MC=0,
				     method = eval(formals(polyG)$method), log = FALSE) {
	              is0 <- t == 0
                      isInf <- is.infinite(t)
	              res <- numeric(n <- length(t))
	              res[is0] <- if(theta==1) 0 else Inf # Note: absdPsi(0, ...) is correct (even for n.MC > 0)
                      res[isInf] <- -Inf
                      n0Inf <- !(is0 | isInf)
                      if(all(!n0Inf)) return(if(log) res else exp(res))
                      t. <- t[n0Inf]
                      if(n.MC > 0) {
                          res[n0Inf] <- absdPsiMC(t, family="Gumbel", theta=theta,
                                                  degree=degree, n.MC=n.MC, log=TRUE)
                      } else {
                          if(theta == 1) {
                              res[n0Inf] <- -t. # independence
                          } else {
                              alpha <- 1/theta
                              lt <- log(t.)
                              res[n0Inf] <- -degree*lt -t.^alpha +
                                  polyG(alpha*lt, alpha=alpha, d = degree,
                                        method=method, log=TRUE)
                          }
                      }
                      if(log) res else exp(res)
                  },
		  ## derivatives of the generator inverse
		  absdiPsi = function(u, theta, degree=1, log=FALSE) {
		      lu <- log(u)
		      switch(degree,
			     ## 1 :
			     if(log) log(theta)+(theta-1)*log(-lu)-lu
			     else theta*(-lu)^(theta-1)/u,
			     ## 2 :
			     if(log) log(theta)+ log(theta-1 - lu) + (theta-2)*log(-lu) - 2*lu
			     else theta * (theta-1 - lu) * (-lu)^(theta-2) / u^2,
			     ## >= 3:
			     stop("not yet implemented for degree > 2"))
		  },
		  ## density of the diagonal
		  dDiag = function(u, theta, d, log=FALSE) {
                      alpha <- 1/theta
                      da <- d^alpha
                      if(log) (da-1)*log(u) + alpha*log(d) else da*u^(da-1)
                  },
		  ## density  Gumbel
		  dacopula = function(u, theta, n.MC=0, checkPar=TRUE,
				      method = eval(formals(polyG)$method), log=FALSE)
              {
                  ## if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
                  if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
		  n <- nrow(u)
		  if(!C.@paraConstr(theta)) {
		      if(checkPar) stop("theta is outside its defined interval")
		      return(rep.int(if(log) -Inf else 0, n))
		  }
		  ## f() := NaN outside and on the boundary of the unit hypercube
		  res <- rep.int(NaN, n)
		  n01 <- u.in.01(u)## indices for which density has to be evaluated
                  ## if(!any(n01)) return(res)
                  if(theta == 1) { res[n01] <- if(log) 0 else 1; return(res) } # independence
                  ## auxiliary results
                  u. <- u[n01,, drop=FALSE]
                  mlu <- -log(u.) # -log(u)
                  lmlu <- log(mlu) # log(-log(u))
                  ## may overflow to Inf: t <- rowSums(C.@iPsi(u., theta))
                  l.iP <- C.@iPsi(u., theta, log=TRUE)
                  ln.t <- copula:::lsum(t(l.iP))
                  ## main part
                  if(n.MC > 0) { # Monte Carlo
                      V <- C.@V0(n.MC, theta)
                      sum. <- rowSums((theta-1)*lmlu + mlu)
                      sum.mat <- matrix(rep(sum., n.MC), nrow=n.MC, byrow=TRUE)
                      ## stably compute log(colMeans(exp(lx)))
                      t <- exp(ln.t)## quickly overflows to +Inf ( <==> lx[.] == -Inf )
                      lx <- d*(log(theta)+log(V)) - log(n.MC) - V %*% t(t) + sum.mat  # matrix of exponents; dimension n.MC x n ["V x u"]
                      res[n01] <- lsum(lx)
                      if(log) res else exp(res)
                  } else { # explicit
                      alpha <- 1/theta
                      ## compute lx = alpha*log(sum(iPsi(u., theta)))
		      lx <- alpha*ln.t ## == log(t^alpha) == log( t^(1/theta) )
                      ## ==== former version [start] (numerically slightly more stable but slower) ====
                      ## im <- apply(u., 1, which.max)
                      ## mat.ind <- cbind(seq_len(n), im) # indices that pick out maxima from u.
                      ## mlum <- mlu[mat.ind] # -log(u_max)
                      ## mlum.mat <- matrix(rep(mlum, d), ncol = d)
                      ## lx <- lmlu[mat.ind] + alpha*log(rowSums((mlu/mlum.mat)^theta)) # alpha*log(sum(iPsi(u, theta)))
                      ## ==== former version [end] ====
                      ## compute sum
                      ls. <- polyG(lx, alpha, d, method=method, log=TRUE)-d*lx/alpha
		      ## the rest
                      ## C() = psi(t(.)) = exp(-t ^ alpha)
		      lnC <- -exp(lx) ## = - t^alpha
		      res[n01] <- lnC + d*log(theta) + rowSums((theta-1)*lmlu + mlu) + ls.
		      if(!log) res[n01] <- exp(res[n01])
                      res
                  }
              },
		  ## score function
		  score = function(u, theta) {
		      if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
		      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
		      stopifnot(C.@paraConstr(theta))
		      stop("The score function is currently not implemented for Gumbel copulas")
		  },
                  ## uscore function
                  uscore = function(u, theta, d, Pmethod=eval(formals(polyG)$method),
                                    P.method=eval(formals(coeffG)$method)) {
                      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      t <- rowSums(C.@iPsi(u, theta))
                      alpha <- 1/theta
                      lx <- alpha*log(t)
                      lP <- polyG(lx, alpha, d, method=Pmethod, log=TRUE)
                      ## compute P' at lx; for now, only implemented all coeffG() methods (see polyG())
                      ## as they already cover a wide range of "stable" values
                      ## Note: this code is copied and adjusted from aux-acopula.R
                      ##     (2), buggy: still using 'lx' (which is a *vector* !)
                      if(missing(P.method)) P.method <- function(d, alpha) { # adjusted meth2012 in polyG()
                          if (d <= 30) "direct"
                          else if (d <= 50) {
                              if (alpha <= 0.8) "direct" else "dsSib.log"
                          }
                          else if (d <= 70) {
                              if (alpha <= 0.7) "direct" else "dsSib.log"
                          }
                          else if (d <= 90) {
                              if (alpha <= 0.5) "direct"
                              else if (alpha >= 0.8) "dsSib.log"
                              ## else if (lx <= 4.08) "pois"
                              else if (lx >= 5.4) "direct"
                              else "dsSib.Rmpfr"
                          }
                          else if (d <= 120) {
                              if (alpha < 0.003) "sort"
                              else if (alpha <= 0.4) "direct"
                              else if (alpha >= 0.8) "dsSib.log"
                              ## else if (lx <= 3.55) "pois"
                              else if (lx >= 5.92) "direct"
                              else "dsSib.Rmpfr"
                          }
                          else if (d <= 170) {
                              if (alpha < 0.01) "sort"
                              else if (alpha <= 0.3) "direct"
                              else if (alpha >= 0.9) "dsSib.log"
                              ## else if (lx <= 3.55) "pois"
                              else "dsSib.Rmpfr"
                          }
                          else if (d <= 200) {
                              if (lx <= 2.56) "pois"
                              else if (alpha >= 0.9) "dsSib.log"
                              else "dsSib.Rmpfr"
                          }
                          else "dsSib.Rmpfr"
                      }
                      ## compute the coefficients (adjusted for P')
                      k <- 1:d
                      l.a.dk. <- log(k) + coeffG(d, alpha, method=P.method, log = TRUE) # new: log(k)
                      logx. <- l.a.dk. + (k-1) %*% t(lx) # new: new l.a.dk. and k-1 instead of k
                      lP. <- lsum(logx.) # P'(lx)
                      ## put the pieces together
                      t.a <- t^alpha
                      mlu <- -log(u)
                      (mlu^(theta-1)*(t.a*(1-exp(lP.-lP))+d/alpha)/t - (theta-1)/mlu + 1) / u
                  },
		  ## nesting constraint
		  nestConstr = function(theta0,theta1) {
		      C.@paraConstr(theta0) &&
		      C.@paraConstr(theta1) && theta1 >= theta0
		  },
		  ## V0 with density dV0 and V01 with density dV01 corresponding to
		  ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
		  V0 = function(n,theta, log=FALSE) {
		      if(theta == 1) {
			  ## Sample from S(1,1,0,1;1)
			  ## with Laplace-Stieltjes transform exp(-t)
			  rep.int(if(log) 0. else 1., n)
		      } else {
                          if(log) stop("'log=TRUE' not yet implemented")
			  alpha <- 1/theta
			  ## Sample from S(alpha,1,(cos(alpha*pi/2))^(1/alpha),0;1)
			  ## with Laplace-Stieltjes transform exp(-t^alpha)
			  rstable1(n, alpha, beta=1,
				   gamma = cospi2(alpha)^(1/alpha))
			  ## Note: calling sequence:
			  ##	   rstable1 == rstable1C (in rstable1.R)
			  ##	   -> rstable_c (in retstable.c) -> rstable_vec (in retstable.c)
			  ##	   -> rstable0 (in retstable.c)
		      }
		  },
		  dV0 = function(x,theta,log = FALSE) C.@dV01(x,1,1,theta,log),
		  V01 = function(V0,theta0,theta1) {
		      alpha <- theta0/theta1
		      if(alpha == 1) {
			  ## Sample from S(1,1,0,V0;1)
			  ## with Laplace-Stieltjes transform exp(-V0*t)
			  V0
		      } else {
			  rstable1(length(V0), alpha, beta=1,
				   gamma = (cospi2(alpha)*V0)^(1/alpha))
			  ## Sample from S(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),0;1)
			  ## with Laplace-Stieltjes transform exp(-V0*t^alpha)
		      }
		  },
		  dV01 = function(x,V0,theta0,theta1,log = FALSE) {
		      stopifnot(length(V0) == 1 || length(x) == length(V0))
		      alpha <- theta0/theta1
		      gamma <- (cospi2(alpha)*V0)^(1/alpha)
		      delta <- V0*(alpha == 1)
		      ## NB: new dstable() is vectorized in (x, gamma, delta) [but not the others]
		      dstable(x, alpha=alpha, beta = 1, gamma=gamma, delta=delta, pm = 1, log=log,
			      tol = 128* .Machine$double.eps)
		  },
		  ## Kendall's tau
		  tau = function(theta) { (theta-1)/theta },
		  iTau = function(tau) { 1/(1-tau) },
		  ## lower tail dependence coefficient lambda_l
		  lambdaL = function(theta) { 0*theta },
		  lambdaLInv = function(lambda) {
		      if(any(lambda != 0))
			  stop("Any parameter for a Gumbel copula gives lambdaL = 0")
		      NA * lambda
		  },
		  ## upper tail dependence coefficient lambda_u
		  lambdaU = function(theta) { 2 - 2^(1/theta) },
		  lambdaUInv = function(lambda) { 1/log2(2-lambda) }
		  )
	C.
    })()# {copGumbel}


### Joe, see Nelsen (2007) p. 116, # 6 #########################################

##' Joe object
copJoe <-
    (function() { ## to get an environment where C. itself is accessible
	C. <- new("acopula", name = "Joe",
		  ## generator
		  psi = function(t, theta) {
		      1 - (-expm1(0-t))^(1/theta)
		      ## == 1 - (1-exp(-t))^(1/theta)
		  },
		  iPsi = function(u, theta, log=FALSE) {
                      res <- -log1p(-(1-u)^theta)
                      if(log) log(res) else res
		  },
		  ## parameter interval
		  paraInterval = interval("[1,Inf)"),## [0.24, Inf) for neg.tau
		  ## absolute value of generator derivatives
		  absdPsi = function(t, theta, degree = 1, n.MC = 0,
				     method= eval(formals(polyJ)$method), log = FALSE)
              {
                  is0 <- t == 0
                  isInf <- is.infinite(t)
                  res <- numeric(n <- length(t))
                  res[is0] <- if(theta==1) 0 else Inf # Note: absdPsi(0, ...) is correct (even for n.MC > 0)
                  res[isInf] <- -Inf
                  n0Inf <- !(is0 | isInf)
                  if(all(!n0Inf)) return(if(log) res else exp(res))
                  t. <- t[n0Inf]
                  if(n.MC > 0) {
                      res[n0Inf] <- absdPsiMC(t, family="Joe", theta=theta,
                                              degree=degree, n.MC=n.MC, log=TRUE)
                  } else {
                      if(theta == 1) {
                          res[n0Inf] <- -t. # independence
                      } else {
                          alpha <- 1/theta
                          mt <- -t.
                          l1mt <- log1mexp(t.) # log(1-exp(-t))
                          sum. <- polyJ(mt-l1mt, alpha, degree, method=method, log=TRUE)
                          res[n0Inf] <- -log(theta) + mt - (1-alpha)*l1mt + sum.
                      }
                  }
                  if(log) res else exp(res)
              },
		  ## derivatives of the generator inverse
		  absdiPsi = function(u, theta, degree=1, log=FALSE) {
		      Iu <- 1-u
		      Iuth <- Iu^theta
		      switch(degree,
			     ## 1 :
			     if(log) log(theta) + (theta-1)*log1p(-u) - log1p(-Iuth)
			     else theta / (Iu/Iuth - Iu),
			     ## 2 :
			     if(log) log(theta) + (theta-2)*log1p(-u) + log(theta-1+Iuth) - 2*log1p(-Iuth)
			     else theta * Iu^(theta-2) * (theta-1+Iuth) / (1-Iuth)^2,
			     ## >= 3:
			     stop("not yet implemented for degree > 2"))
		  },
		  ## density of the diagonal
		  dDiag = function(u, theta, d, log=FALSE) {
		      ##  d* x^(d-1) * ((1-x^d)/(1-x)) ^ (1/theta-1), where  x = 1 - (1-u)^theta
		      ##  we use  (1-x^d)/(1-x) === circRat((1-u)^theta, d)
		      Ix <- (1-u)^theta
		      x <- 1-Ix
		      if(log) log(d)+ (d-1)*log(x) + (1/theta-1)*log(circRat(Ix, d))
		      ## FIXME? for log-case: log(circRat(Ix, d)) = log1p(-x^d)-theta*log1p(-u))
		      else d* x^(d-1) * circRat(Ix, d)^(1/theta-1)
		  },
		  ## density  Joe
		  dacopula = function(u, theta, n.MC=0, checkPar=TRUE,
				      method = eval(formals(polyJ)$method), log = FALSE)
              {
                  ## if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
                  if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
		  n <- nrow(u)
		  if(d > 2 && !C.@paraConstr(theta)) {
		      if(checkPar) stop("theta is outside its defined interval")
		      return(rep.int(if(log) -Inf else 0, n))
		  }
		  ## f() := NaN outside and on the boundary of the unit hypercube
		  res <- rep.int(NaN, n)
		  n01 <- u.in.01(u)## indices for which density has to be evaluated
                  ## if(!any(n01)) return(res)
                  if(theta == 1) { res[n01] <- if(log) 0 else 1; return(res) } # independence
                  ## auxiliary results
                  u. <- u[n01,, drop=FALSE]
                  l1_u <- rowSums(log1p(-u.)) # sum_j log(1-u_j)
                  lh <- rowSums(log1p(-(1-u.)^theta)) # rowSums(log(1-(1-u)^theta)) = log(h)
                  ## main part
                  if(n.MC > 0) { # Monte Carlo
                      V <- C.@V0(n.MC, theta)
                      sum. <- (theta-1)*l1_u
                      sum.mat <- matrix(rep(sum., n.MC), nrow=n.MC, byrow=TRUE)
                      ## stably compute log(colMeans(exp(lx)))
                      ## matrix of exponents; dimension n.MC x n ["V x u"] :
                      lx <- d*(log(theta) + log(V)) + (V-1) %*% t(lh) +
                          sum.mat - log(n.MC)
                      res[n01] <- lsum(lx)
                  } else {
                      alpha <- 1/theta
                      l1_h <- log1mexp(-lh) # log(1-h)
                      lh_l1_h <- lh - l1_h # log(h/(1-h))
                      res[n01] <- (d-1)*log(theta) + (theta-1)*l1_u -
                          (1-alpha)*l1_h + polyJ(lh_l1_h, alpha, d, method=method,
                                                 log=TRUE)
                  }
                  if(log) res else exp(res)
              },
		  ## score function
		  score = function(u, theta, method=eval(formals(polyJ)$method)) {
		      if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
		      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
		      if(d > 2) stopifnot(C.@paraConstr(theta))
		      l1_u <- rowSums(log1p(-u)) # log(1-u)
		      u.th <- (1-u)^theta # (1-u)^theta
		      lh <- rowSums(log1p(-u.th)) # rowSums(log(1-(1-u)^theta)) = log(h)
		      l1_h <- log1mexp(-lh) # log(1-h)
		      lh_l1_h <- lh - l1_h # log(h/(1-h))
		      b <- rowSums(-l1_u*u.th/(1-u.th))
		      lP <- polyJ(lh_l1_h, alpha, d, method=method, log=TRUE)
		      k <- 1:d
		      alpha <- 1/theta
		      s <- alpha * unlist(lapply(k, function(k.) sum(1/(theta*(1:k.)-1))))
		      ls <- log(s. + (k-1)*exp(log(b) - l1_h))
		      l.a.k <- log(Stirling2.all(d)) + lgamma(k-alpha) - lgamma(1-alpha) + ls
		      lQ <- lsum(l.a.k + (k-1) %*% t(lh_l1_h))
		      (d-1)/theta + rowSums(l1_u) - l1_h/theta^2 + (1-1/theta)*
			  lh_l1_h*b + exp(lQ-lP)
		  },
                  ## uscore function
                  uscore = function(u, theta, d, method=eval(formals(polyJ)$method)) {
                      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      alpha <- 1/theta
                      ## compute the log of the argument of P and P'
                      omu <- 1-u # 1-u
                      u.th <- omu^theta # (1-u)^theta
		      lh <- rowSums (log1p(-u.th)) # rowSums(log(1-(1-u)^theta)) = log(h)
		      l1_h <- log1mexp(-lh) # log(1-h)
		      lh_l1_h <- lh - l1_h # log(h/(1-h))
                      ## compute log(P(log(h/(1-h))))
                      lP <- polyJ(lh_l1_h, alpha, d, method=method, log=TRUE)
                      ## compute log(P'(log(h/(1-h))))
                      ## Note: this is similar to polyJ() (see there for the comments!)
                      k <- 2:d # 2:d instead of 1:d
                      l.a.k <- log(Stirling2.all(d)) + lgamma(k-alpha) - lgamma(1-alpha) # log(a_{dk}(theta)*(k+1)), k = 1,..,d; note: these are not the a's of Hofert, Maechler, McNeil (2013); see polyJ()
                      l.a.k. <- log(k-1) + l.a.k
                      B <- l.a.k. + (k-2) %*% t(lh_l1_h) # new: new l.a.k. and k-2 instead of k-1
                      ## the following part is taken from polyJ() (but only the log cases)
                      lP. <- switch(method,
                                 "log.poly" = {
                                     lsum(B)
                                 },
                                 "log1p" = {
                                     im <- apply(B, 2, which.max) # indices (vector) of maxima
                                     n <- length(lh_l1_h) ; d1 <- d-1L
                                     max.B <- B[cbind(im, seq_len(n))] # get max(B[,i])_{i=1,..,n} == apply(B, 2, max)
                                     B.wo.max <- matrix(B[unlist(lapply(im, function(j) k[-j])) +
                                                          d*rep(0:(n-1), each = d1)], d1, n) # matrix B without maxima
                                     max.B + log1p(colSums(exp(B.wo.max - rep(max.B, each = d1))))
                                 },
                                 "poly" = {
                                     log(colSums(exp(B)))
                                 },
                                 stop(gettextf("unsupported method '%s' in polyJ", method)), domain=NA)
                      ## put the pieces together
                      theta*omu^(theta-1)/(1-u.th) * (1-alpha+exp(-l1_h)*exp(lP.-lP)) - (theta-1)/omu
                  },
		  ## nesting constraint
		  nestConstr = function(theta0,theta1) {
		      C.@paraConstr(theta0) &&
		      C.@paraConstr(theta1) && theta1 >= theta0
		  },
		  ## V0 with density dV0 and V01 with density dV01 corresponding to
		  ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
		  V0 = function(n,theta, log=FALSE) {
                      if(log) stop("'log=TRUE' not yet implemented")
                      else rSibuya(n, 1/theta) },
		  dV0 = function(x,theta,log = FALSE) {
		      if(log) lchoose(1/theta,x) else abs(choose(1/theta,x))
		  },
		  V01 = function(V0, theta0, theta1, approx = 10000) {
		      ## approx is the largest number of summands before asymptotics is used
		      alpha <- theta0/theta1
		      rF01Joe(V0, alpha, approx)
		  },
		  dV01 =
		  function(x, V0, theta0, theta1,
			   method= eval(formals(dsumSibuya)$method), log = FALSE) {
		      stopifnot(length(V0) == 1 || length(x) == length(V0))
		      ## also holds for theta0 == theta1
		      ## note: this is numerically challenging
		      dsumSibuya(x, V0, theta0/theta1, method=method, log=log)
		  },
		  ## Kendall's tau
		  ## noTerms: even for theta==0, the approximation error is < 10^(-5)
                  ## MM: "FIXME" , using  http://dlmf.nist.gov/2.10#E1  (or better?)
                  ## + maxima   integrate(1/(x*(t*x+2)*(t*x+2-t)), x)
		  tau = tauJoe,
		  iTau = function(tau, tol = .Machine$double.eps^0.25, ...) {
		      sapply(tau,function(tau) {
			  r <- safeUroot(function(th) C.@tau(th) - tau,
					 interval = c(1, 98),
					 Sig = +1, tol=tol, check.conv=TRUE, ...)
			  r$root
		      })
		  },
		  ## lower tail dependence coefficient lambda_l
		  lambdaL = function(theta) { 0*theta },
		  lambdaLInv = function(lambda) {
		      if(any(lambda != 0))
			  stop("Any parameter for a Joe copula gives lambdaL = 0")
		      NA * lambda
		  },
		  ## upper tail dependence coefficient lambda_u
		  lambdaU = function(theta) { 2-2^(1/theta) },
		  lambdaUInv = function(lambda) { log(2)/log(2-lambda) }
		  )
	C.
    })()# {copJoe}


### naming stuff ###############################################################

cNms <- c("copAMH", "copClayton", "copFrank", "copGumbel", "copJoe")
## == dput(ls("package:copula",pat="^cop"))
nmsC <- unlist(lapply(cNms, function(.)get(.)@name))
sNms <- abbreviate(nmsC, 1)
## keep these {hidden, for now}:
.ac.shortNames <- structure(sNms, names = nmsC)
.ac.longNames  <- structure(nmsC, names = sNms)
.ac.objNames   <- structure(cNms, names = nmsC)
.ac.classNames <- structure(paste0(tolower(nmsC), "Copula"), names = nmsC)
rm(cNms, nmsC, sNms)

