###  Copyright (C) 2011-2012 Berwin A. Turlach <Berwin.Turlach@gmail.com>
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
Full <- function(x, y, w, K, par, trace, plot.it,
                 control, ptype, ctype){

  if(missing(w) || is.null(w)){
    RSS <- makeRSS(x, y, ptype, ctype)
    Grad <- makeGrad(x, y, ptype, ctype)
    Hess <- makeGradHessFull(x, y, ptype, ctype)
  }else{
    RSS <- makewRSS(x, y, w, ptype, ctype)
    Grad <- makewGrad(x, y, w, ptype, ctype)
    Hess <- makewGradHessFull(x, y, w, ptype, ctype)
  }
  
  if(plot.it)
    xgr <- seq(from=-1, to=1, length=401)
  if(trace){
    cat("\n----------------------\n")
    cat("RSS:\t", RSS(par), "\n")
    cat("\n")
    print(zapsmall(rbind(par,grad=Grad(par))))
    cat("\n")
    print(zapsmall(evalCoef(par, ptype, ctype)))
  }
  
  maxiter <- control$maxiter
  tol <- control$tol
  
  converged <- FALSE
  terminate <- FALSE
  lam <- 0.1
  for(ol in 1:maxiter){
    cRSS <- RSS(par)

    tt <- Hess(par)

    j <- 1
    repeat{
##      yy <- pmax(diag(tt$Hess), 0) + 1
      yy <- diag(tt$Hess)
      yy[yy<1e-7] <- 1
      zz <- try(tmp <- par + solve(tt$Hess+lam*diag(yy), -tt$Grad),
                silent=TRUE)
      if(inherits(zz, "try-error")){
        nRSS <- cRSS+1
      }else{
        nRSS <- RSS(tmp)
      }
      
      if(nRSS <= cRSS){
        par <- tmp
        if(j==1) lam <- lam/10
        break
      }else{
        j <- j+1
        lam <- lam*10
      }
      if(lam > 1e10){
        warning("lamda is getting too large")
        terminate <- TRUE
        break
      }
    }
    grad <- Grad(par)

    if(terminate)
      break
    
    if(trace && (ol %% trace == 0)){
      cat("\n----------------------\n")
      cat("Iteration:\t", ol, "\nRSS:\t", RSS(par), "\n")
      cat("\n")
      print(zapsmall(rbind(par,grad=grad)))
      cat("\n")
      print(zapsmall(evalCoef(par, ptype, ctype)))
    }
    if(plot.it && (ol %% plot.it == 0)){
      ygr <- evalPolMonPar(par, xgr, ptype, ctype)
      lines(xgr,ygr, col="blue")
    }
    if(max(abs(grad)) < tol){
      converged <- TRUE
      break
    }
  }
  names(grad) <- names(par)
  list(par=par, grad=grad, beta=evalCoef(par, ptype, ctype),
       RSS=RSS(par), niter=ol, converged=converged)
}

FullConstr <- function(x, y, w, K, par, trace, plot.it,
                       control, ptype, ctype){


  if(missing(w) || is.null(w)){
    RSS <- makeRSS(x, y, ptype, ctype)
    Grad <- makeGrad(x, y, ptype, ctype)
    Hess <- makeGradHessFull(x, y, ptype, ctype)
  }else{
    RSS <- makewRSS(x, y, w, ptype, ctype)
    Grad <- makewGrad(x, y, w, ptype, ctype)
    Hess <- makewGradHessFull(x, y, w, ptype, ctype)
  }
  
  if(plot.it)
    xgr <- seq(from=-1, to=1, length=401)
  if(trace){
    cat("\n----------------------\n")
    cat("RSS:\t", RSS(par), "\n")
    cat("\n")
    print(zapsmall(rbind(par,grad=Grad(par))))
    cat("\n")
    print(zapsmall(evalCoef(par, ptype, ctype)))
  }

  Amat <- matrix(rep(1,K), nrow=1)
  ii <- 2 + 2*(1:K)
  Aind <- rbind(Amat, ii)
  
  maxiter <- control$maxiter
  tol <- control$tol
  TOL <- .Machine$double.eps*1000
  isncpar <- c(TRUE, TRUE, rep(c(TRUE, FALSE), K))

  converged <- FALSE
  terminate <- FALSE
  lam <- 0.1
  for(ol in 1:maxiter){
    cRSS <- RSS(par)

    tt <- Hess(par)

    j <- 1
    repeat{
##      yy <- pmax(diag(tt$Hess), 0) + 1
      yy <- diag(tt$Hess)
      yy[yy<1e-7] <- 1
      bvec <- -par[ii]
      zz <- try(tmp <- solve.QP.compact(tt$Hess+lam*diag(yy), -tt$Grad,
                                        Amat, Aind, bvec),
                silent=TRUE)
      if(inherits(zz, "try-error")){
        nRSS <- cRSS+1
      }else{
        tmp <- par+tmp$solution
        nRSS <- RSS(tmp)
      }
      
      if(nRSS <= cRSS){
        par <- tmp
        if(j==1) lam <- lam/10
        break
      }else{
        j <- j+1
        lam <- lam*10
      }
      if(lam > 1e10){
        warning("lamda is getting too large")
        terminate <- TRUE
        break
      }
    }
    grad <- Grad(par)

    if(terminate)
      break
    
    if(trace && (ol %% trace == 0)){
      cat("\n----------------------\n")
      cat("Iteration:\t", ol, "\nRSS:\t", RSS(par), "\n")
      cat("\n")
      print(zapsmall(rbind(par,grad=grad)))
      cat("\n")
      print(zapsmall(evalCoef(par, ptype, ctype)))
    }
    if(plot.it && (ol %% plot.it == 0)){
      ygr <- evalPolMonPar(par, xgr, ptype, ctype)
      lines(xgr,ygr, col="blue")
    }
    if(max(abs(grad[isncpar | abs(par)>TOL])) < tol){
      converged <- TRUE
      break
    }
  }
  names(grad) <- names(par)
  list(par=par, grad=grad, beta=evalCoef(par, ptype, ctype),
       RSS=RSS(par), niter=ol, converged=converged)
}
