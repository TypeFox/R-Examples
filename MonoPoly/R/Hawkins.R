###  Copyright (C) 2012 Berwin A. Turlach <Berwin.Turlach@gmail.com>
###
###  This program is free software; you can redistribute it and/or modify
###  it under the terms of the GNU General Public License as published by
###  the Free Software Foundation; either version 2 of the License, or
###  (at your option) any later version.
###
###  This program is distributed in the hope that it will be useful,
###  but WITHOUT ANY WARRANTY; without even the implied warranty of
###  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###  GNU General Public License for more details.
###
###  You should have received a copy of the GNU General Public License
###  along with this program; if not, write to the Free Software
###  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
###  USA.
###
qrpolymat <- function(x, w, q, TolQR){
  n <- length(x)
  Xmat <- matrix(0, nrow=n, ncol=q+1)

  Xmat[,1] <- 1
  for(i in 2:(q+1))
    Xmat[,i] <- Xmat[,i-1]*x
  Xmat <- w*Xmat

  Xmatqr <- qr(Xmat, tol=TolQR)
  if(Xmatqr$rank != q+1){
    cat("Rank of QR factorisation says rank is", Xmatqr$rank, "not", q+1, "\n")
    XQ <- qr.Q(Xmatqr)
    XR <- qr.R(Xmatqr)
    tt <- crossprod(XQ)
    cat("Is Q matrix orthogonal:\n")
    print(max(abs(diag(tt)-1)))
    diag(tt) <- 0
    print(max(abs(tt)))
    cat("Is QR equal to design matrix:\tmax(abs(Xmat-Q%*%R))\n")
    print(max(abs(Xmat-XQ%*%XR)))
  }

  list(q=qr.Q(Xmatqr), r=qr.R(Xmatqr))
} 

hawkins <- function(x, y, w, K, trace, plot.it, control){

  if(plot.it)
    plot(x,y)

  TOL <- control$tol1
  MAXITER <- control$maxiter
  TOL1 <- control$tol2
  TOL2 <- control$tolqr
  
  q <- 2*K+1

  if(missing(w) || is.null(w)){
    w <- rep.int(1, length(x))
  }
  wsq <- sqrt(w)
  tt <- qrpolymat(x, wsq, q, TOL2)
  beta.u <- as.vector(backsolve(tt$r, crossprod(tt$q, wsq*y)))
  names(beta.u) <- paste("beta",0:q, sep="")
  if(trace){
    cat("Coefficients of unconstrained fit:\n")
    print(beta.u)
    cat("RSS:\t", sum(w*(y-evalPol(x,beta.u))^2), "\n\n")
  }
  p1dbeta.u <- beta.u[-1]*(1:q)

  if(plot.it){
    xgr <- seq(-1,1,length=401)
    lines(xgr, evalPol(xgr, beta.u), col="green")
  }
  
  Dmat <- solve(tt$r)
  dvec <- crossprod(tt$r, crossprod(tt$q,wsq*y))
  
  lambda <- xstar.set <- NULL
  beta <- beta.u
  iter <- 0
  if(K>0){
    converged <- FALSE
    while(TRUE){
      if(iter > MAXITER) break
      p1d <- beta[-1]*(1:q)
      p2d <- p1d[-1]*1:(q-1)
      
      p1d.roots <- polyroot(p1d)
      p2d.roots <- polyroot(p2d)

      p1d.troots <- c(-1,1)*max(Mod(p1d.roots))*1.1
      evtr <- evalPol(p1d.troots,p1d)
      p2d.rroots <- which(abs(Im(p2d.roots))<TOL1)
      p2d.rroots <- Re(p2d.roots[p2d.rroots])
      evrr <- evalPol(p2d.rroots,p1d)
      iter <- iter+1
      if(trace){
        cat("\nIteration", iter, "\n")
        cat("\tRoots of the first derivative:\n")
        print(p1d.roots)
        cat("\tFirst derivative at points beyond largest root:\n")
        print(evtr)
        cat("\tRoots of the second derivative:\n")
        print(p2d.roots)
        cat("\tReal roots of the second derivative:\n")
        print(p2d.rroots)
        cat("\tFirst derivative at real roots of the second derivative:\n")
        print(evrr)
      }
      
      if(any(evrr < -TOL) | any(evtr < -TOL)){
        if(any(evtr < -TOL)){
          ind <- which.min(evtr)
          xstar <- p1d.troots[ind]
        }else{
          ind <- which.min(evrr)
          xstar <- p2d.rroots[ind]
        }
        xstar.set <- union(xstar.set, xstar)
        if(trace)
          cat("\tEnforcing hip at:", xstar.set, "\n")
        
        Hmat <- matrix(0, nrow=length(xstar.set), ncol=q+1)
        Hmat[,2] <- 1
        for(l1 in 3:(q+1))
          Hmat[,l1] <- (l1-1)*xstar.set^(l1-2)
        
        bvec <- rep(0,NROW(Hmat))
        res <- solve.QP(Dmat, dvec, Amat=t(Hmat), bvec, factorized=TRUE)
        lambda <- res$Lagrangian
        if(trace)
          cat("\tLagrangian:", lambda,"\n")
        
        if(any(lambda<=0)){
          if(trace)
            cat("\t**Removing xstars with non-positive Lagrangian**\n")
          ind <- which(lambda<=0)
          xstar.set <- xstar.set[-ind]
          lambda <- lambda[-ind]
        }
        beta <- res$solution
        
        HmatBeta <- evalPol(xstar.set, beta[-1]*(1:q))
        stopifnot(max(abs(HmatBeta*lambda))< TOL1)
        if(plot.it)
          lines(xgr, evalPol(xgr, beta), col="blue")
      }else{
        converged <- TRUE
        break
      }
    }
  }else{
    if(beta[2] < -TOL){
      iter <- iter+1
      if(is.null(w)){
        beta <- c(mean(y),0)
      }else{
        beta <- c(weighted.mean(y, w), 0)
      }
    }
    converged <- TRUE
  }
  names(beta) <- paste("beta",0:q, sep="")
  RSS <-  sum(w*(y-evalPol(x,beta))^2)
  if(trace){
    cat("\n\nConstrained beta:\n")
    print(beta)
    cat("RSS:\t", RSS, "\n")
    if(length(xstar.set)>0){
      cat("Hips:\n")
      print(xstar.set)
      cat("Lagrange multiplier:\n")
      print(lambda)
    }
  }
  if(plot.it)
    lines(xgr, evalPol(xgr, beta), col="black")
  list(beta=beta, beta.unconstrained=beta.u, RSS=RSS,
       hips=xstar.set, Lagrangian=lambda,
       niter=iter, converged=converged)
}
