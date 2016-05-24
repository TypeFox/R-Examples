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


##' Outer power transformation of Archimedean copulas
##'
##' @title Outer power Archimedean copulas
##' @param copbase a "base" copula, i.e. of class "acopula";
##'    must be one of the predefined families
##' @param thetabase the (univariate) parameter 'theta' for the base copula
##' @return a new "acopula" object; the outer power copula [with generator psi(t^(1/theta)]
##' @author Marius Hofert
opower <- function(copbase, thetabase) {
    ## create object with name in here so it's in the environment and we can access it
    C. <- new("acopula", name = paste("opower", copbase@name, sep=":"),
              ## generator
              psi = function(t,theta) { copbase@psi(t^(1/theta), thetabase) },
	      iPsi = function(t,theta, log=FALSE) {
		  if(log) theta * copbase@iPsi(t, thetabase, log=TRUE)
		  else copbase@iPsi(t, thetabase)^theta
	      },
	      ## parameter interval
	      paraInterval = interval("[1,Inf)"),
	      ## nesting constraint
	      nestConstr = function(theta0,theta1) {
		  copbase@paraConstr(theta0) &&
		  copbase@paraConstr(theta1) && theta1 >= theta0
	      },
	      ## absolute value of generator derivatives
	      absdPsi = function(t, theta, degree=1, n.MC=0,
				 method=c("stirling", "binomial.coeff"), log=FALSE)
	  {
	      if(theta == 1) return(copbase@absdPsi(t, theta, degree=degree,
		 n.MC=n.MC, log=log))	# copbase case
	      is0 <- t == 0
	      isInf <- is.infinite(t)
	      res <- numeric(n <- length(t))
	      res[is0] <- Inf # Note: absdPsi(0, ...) is correct (even for n.MC > 0)
	      res[isInf] <- -Inf
	      n0Inf <- !(is0 | isInf)
	      if(all(!n0Inf)) return(if(log) res else exp(res))
	      t. <- t[n0Inf]
	      if(n.MC > 0) {
		  res[n0Inf] <- absdPsiMC(t, family=C., theta=theta,
					  degree=degree, n.MC=n.MC, log=TRUE)
	      } else {
		  k <- 1:degree # compute everything once for 1:degree
		  beta <- 1/theta
		  t.beta <- t.^beta   # beta = 1/theta for psi(t^beta)
		  ## in principle, it would be efficient to vectorize
		  ## the absdPsi slots also in the parameter "degree",
		  ## but that would be even more complicated
		  labsdPsi <- do.call(rbind,
				      lapply(k, function(k.)
					     copbase@absdPsi(t.beta,
							     theta=thetabase,
							     degree=k.,
							     log=TRUE))) # (degree,n)-matrix
		  bklt <- k %*% t(log(t.beta)) # (degree,n)-matrix
		  method <- match.arg(method)
		  res[n0Inf] <-
		      switch(method,
			     "stirling" = { # much faster & less prone to errors
				 s <- Stirling1.all(degree)
				 b.one.j <- function(j.) {
				     k <- 1:j.
				     signs <- (-1)^k
				     lS <- log(Stirling2.all(j.))
				     a <- lS + labsdPsi[k,, drop=FALSE] + bklt[k,, drop=FALSE]
				     (-beta)^j. * abs(s[j.]) * colSums(signs*exp(a))
				 }
				 ## => returns a vector of length n containing the values for one j. and all t
				 j <- 1:degree
				 b <- do.call(rbind, lapply(j, FUN=b.one.j)) # (degree, n)-matrix
				 -degree*log(t.)+log(colSums(b))
			     },
			     "binomial.coeff" = {
				 ## outer sum
				 lfac <- lfactorial(0:degree) # log(0!), log(1!), .., log(degree!)
				 log.b.one.j <- function(j.) {
				     k <- j.:degree
				     a <- labsdPsi[k,, drop=FALSE] + bklt[k,, drop=FALSE] - lfac[j.+1] - lfac[k-j.+1] # (degree-j.+1, n)-matrix
				     ls. <- lsum(a) # length = n
				     lchoose(beta*j., degree) + ls. # note: the lchoose() can be non-finite with this approach!
				 }
				 ## => returns a vector of length n containing the values for one j. and all t
				 j <- 1:degree
				 b <- do.call(rbind, lapply(j, FUN=log.b.one.j)) # (degree, n)-matrix
				 signs <- signFF(beta, j, degree)
				 lfac[degree+1] - degree*log(t.) + lssum(b, signs, strict=FALSE)
			     }, stop(gettextf("unsupported method '%s' in absdPsi",
					     method), domain=NA)) # end{switch}
	      }
	      if(log) res else exp(res)
	  },
              ## derivatives of the generator inverse
              absdiPsi = function(t, theta, log=FALSE) {
                  if(theta == 1) return(copbase@absdiPsi(t, theta, log=log)) # copbase case
                  if(log) {
                      log(theta)+(theta-1)*log(copbase@iPsi(t,thetabase))+
                          copbase@absdiPsi(t, thetabase,log=TRUE)
                  } else {
                      theta*copbase@iPsi(t,thetabase)^(theta-1)*
                          copbase@absdiPsi(t, thetabase,log=FALSE)
                  }
              },
              ## density
              dacopula = function(u, theta, n.MC=0, log=FALSE)
          {
              stopifnot(C.@paraConstr(theta))
              if(theta == 1) return(copbase@dacopula(u, theta, n.MC=n.MC, log=log)) # copbase case
              if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
              if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
              ## f() := NaN outside and on the boundary of the unit hypercube
              res <- rep.int(NaN, n <- nrow(u))
              n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
              if(!any(n01)) return(res)
              ## auxiliary results
              u. <- u[n01,, drop=FALSE]
              psiI <- rowSums(C.@iPsi(u.,theta))
              res[n01] <- C.@absdPsi(psiI, theta, degree=d, n.MC=n.MC, log=TRUE) +
                  rowSums(C.@absdiPsi(u., theta, log=TRUE))
              if(log) res else exp(res)
          },
              ## score function
              score = function(u, theta) {
	          stopifnot(C.@paraConstr(theta))
                  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
                  if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                  stop("The score function is currently not implemented for outer power copulas")
              },
              ## V0 and V01
              V0 = function(n, theta) {
	          V0base <- copbase@V0(n, thetabase)
                  if(theta == 1) return(V0base) # the copula is copbase with thetabase
                  alpha <- 1/theta
                  ## Sample from S(alpha,1,(cos(alpha*pi/2))^(1/alpha),0;1)
                  ## with Laplace-Stieltjes transform exp(-t^alpha)
                  S <- rstable1(n, alpha, beta=1,
                                gamma = cospi2(alpha)^(1/alpha))
                  S*V0base^theta
              },
              dV0 = function(x, theta, log=FALSE) {
                  stop("not implemented; it's the density of SV^theta, where V ~ F with LS[F] = copbase@psi and S ~ S(1/theta, 1, cos^theta(pi/(2*theta)), I_{theta==1}; 1)")
              },
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
              dV01 = function(x, V0, theta0, theta1, log=FALSE) {
                  copGumbel@dV01(x, V0, theta0, theta1, log=log)
              },
              ## Kendall's tau
              tau = function(theta) {
                  1-(1-copbase@tau(thetabase))/theta
              },
              iTau = function(tau) {
                  taubase <- copbase@tau(thetabase)
                  if(tau >= taubase) (1-taubase)/(1-tau)
                  else {
                      stop("The provided tau has to be >= taubase")
                      NA * tau
                  }
              },
              ## lower tail dependence coefficient lambda_l
              lambdaL = function(theta) {
                  if(copbase@name=="Clayton") 2^(-1/(thetabase*theta)) else 0*theta
              },
              lambdaLInv = function(lambda) {
                  if(copbase@name=="Clayton") {
                      if(lambda >= 2^(-1/thetabase)) -1/(thetabase*log2(lambda))
                      else {
                          stop("The provided lambda has to be >= 2^(-1/thetabase)")
                          NA * lambda
                      }
                  } else {
                      if(any(lambda != 0))
                          stop("Any parameter for this outer power Archimedean copula gives lambdaL = 0")
                      NA * lambda
                  }
              },
              ## upper tail dependence coefficient lambda_u
              lambdaU = function(theta) {
                  2 - 2^(1/if(copbase@name %in% c("Gumbel", "Joe")) thetabase*theta else theta)
              },
              lambdaUInv = function(lambda) {
                  if(copbase@name %in% c("Gumbel", "Joe")) {
                      if(lambda >= 2-2^(1/thetabase)) 1/(thetabase*log2(2-lambda))
                      else {
                          stop("The provided lambda has to be >= 2-2^(1/thetabase)")
                          NA * lambda
                      }
                  } else 1/log2(2-lambda)
              }
              )
    C.
}
