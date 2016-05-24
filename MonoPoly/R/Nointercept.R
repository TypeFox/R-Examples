###  Copyright (C) 2015 Berwin A. Turlach <Berwin.Turlach@gmail.com>
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
monpol.noint <- function(x=x, y=y, w=w, deg.is.odd, K, start,
                         a, b, trace,
                         plot.it=plot.it, control=control){

  no.w <- FALSE
  if(missing(w) || is.null(w)){
    no.w <- TRUE
    w <- rep.int(1, length(x))
  }

  minx <- 0
  sclx <- max(abs(x))
  x <- x/sclx
  a <- a/sclx
  b <- b/sclx

  ## other code (e.g. fitted.monpol) requires that y is sclaed
  ## according to the formula 
  ##            y <- 2*(y-miny)/sly - 1
  ## as we do not want to shift y (only rescale), miny is defined so
  ## that the shift becomes zero, i.e. miny <- -scly/2
  scly <- max(abs(y))*2
  miny <- -scly/2
  y <- 2*y/scly
  
  if( K == 0){
    tmp <- lm.wfit(x, y, w)
    beta <- c(0,coef(tmp))
    names(beta) <- c("beta0", "beta1")
    res <- list(beta=beta, niter=0, converged=TRUE)

    res$RSS <- sum(w*tmp$residuals^2)
    if(plot.it){
      plot(x, y, xlab="x", ylab="y")
      abline(beta, col="green")
    }
  }else{
    xgr <- seq(-1,1,length=401)

    n <- length(x)
    q <- 2*K + deg.is.odd
    Xmat <- matrix(0, nrow=n, ncol=q)
    Xmat[,1] <- x
    for(i in 2:q)
      Xmat[,i] <- Xmat[,i-1]*x

    tmp <- Xmat * Xmat
    scnd.deriv <- colSums(w*tmp)
    Xmat <- w*Xmat
    
    maxiter <- control$maxiter
    beta <- rep(0, q+1)
    names(beta) <- paste0("beta", 0:q)

##    ii <- seq(from=2, to=q+1, by=2)
    ii <- (start+1):2
    converged <- FALSE
    niter <- 0
    while(!converged){
      beta.old <- beta
      RSS.old <- sum(w*(y-evalPol(x,beta.old))^2)
      if(trace){
        cat("Iteration", niter, "\n")
        cat("beta:", beta, "\n")
        cat("RSS:", RSS.old, "\n")
      }
      if(plot.it){
        plot(x,y, xlab="x", ylab="y")
        lines(xgr, evalPol(xgr, beta), col="green")
      }
      
      jjMaxReached <- FALSE
      for(ind in ii){
        if(trace){
          cat("\tUpdating component", ind, "\t")
        }
        frst.deriv <- crossprod(Xmat[,ind-1],y-evalPol(x,beta))
        beta.try <- beta
        step.len <- 1
        step <- frst.deriv/scnd.deriv[ind-1]
        jj <- 1
        while( jj <= 20 ){
          beta.try[ind] <- beta[ind] + step.len * step
          if(ismonotone(beta.try, a, b)){
            beta <- beta.try
            break()
          }
          jj <- jj+1
          step.len <- step.len/2
        }
        if(jj == 20 ) {
          warning(paste("'step.len' reduced to ", step.len))
          jjMaxReached <- TRUE
        }
        if(trace) {
          rr <- crossprod(Xmat[,ind-1],y-evalPol(x,beta))
          if(jj == 1){
            cat("full step taken, X^Tr=", rr, "\n")
          }else{
            cat("part step taken, X^Tr=", rr, "\n")
          }
        }
        if(plot.it){
          lines(xgr, evalPol(xgr, beta), col="green")
        }
      }
      ii <- (q+1):2
      niter <- niter+1
      RSS <- sum(w*(y-evalPol(x,beta))^2)
      if(trace){
        cat("beta:", beta,"\n")
        cat("RSS:", RSS,"\n")
        cat("jjMaxReached", jjMaxReached,"\n")
      }
      if( !jjMaxReached && 
         sum((beta.old-beta)^2)/(1 + sum(beta^2)) < 1e-7 &&
         abs(RSS-RSS.old)/abs(RSS) < 1e-7)
        converged <- TRUE
      if(trace) cat("\n\n")
      if(niter == maxiter){
        warning("Maximum number of iterations reached")
        break()
      }
    }
    res <- list(beta=beta, niter=niter, converged=converged)
    resid <- y-evalPol(x,beta)
    res$RSS <- sum(w*resid^2)
  }

  tmp <- diag(rev(res$beta))
  storage.mode(tmp) <- "double"
  tmp <- .Fortran(.MP_transf,
                  c = tmp,
                  as.integer(NROW(tmp)), as.integer(NCOL(tmp)),
                  scl = as.double(1/sclx),
                  mid = as.double(0))$c
  tmp <- rev(tmp[,1]) * scly/2
  names(tmp) <- names(res$beta)
  res$beta.raw <- tmp
  if(!no.w)
    res$weights <- w

  res$fitted.values <- evalPol(x, res$beta)
  attributes(res$fitted.values) <- attributes(y)
  res$residuals <- y - res$fitted.values
  
  c(res,
    list(K=K, minx=minx, sclx=sclx, miny=miny, scly=scly,
         algorithm="NoIntercept"))
}
