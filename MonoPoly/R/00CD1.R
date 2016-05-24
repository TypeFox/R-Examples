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
CD1 <- function(x, y, w, K, par, trace, plot.it,
                control, ptype, ctype){

  if(missing(w) || is.null(w)){
    RSS <- makeRSS(x, y, ptype, ctype)
    Grad <- makeGrad(x, y, ptype, ctype)
    i <- 1
    GradHess <- makeGradHessBlk1(x, y, par, i, ptype, ctype)
    envGradHess <- environment(GradHess)
  }else{
    RSS <- makewRSS(x, y, w, ptype, ctype)
    Grad <- makewGrad(x, y, w, ptype, ctype)
    i <- 1
    GradHess <- makewGradHessBlk1(x, y, w, par, i, ptype, ctype)
    envGradHess <- environment(GradHess)
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
  for(ol in 1:maxiter){
    for(i in 1:(2*K)){
      j <- i+2L
      assign("i", i, envir=envGradHess)
      
      tt <- GradHess(par)
      par[j]  <- par[j] - tt[1L]/tt[2L]
      tt <- GradHess(par)
      par[j]  <- par[j] - tt[1L]/tt[2L]
    }
    tmp <- lm.fit(cbind(1, evalPolMonPartilde(par, x, ptype, ctype)), y)
    par[1:2] <- as.vector(coef(tmp))
    assign("a", par[2L], envir=envGradHess)
    grad <- Grad(par)
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

CD1Constr <- function(x, y, w, K, par, trace, plot.it,
                      control, ptype, ctype){

  if(missing(w) || is.null(w)){
    RSS <- makeRSS(x, y, ptype, ctype)
    Grad <- makeGrad(x, y, ptype, ctype)
    i <- 1
    GradHess <- makeGradHessBlk1(x, y, par, i, ptype, ctype)
    envGradHess <- environment(GradHess)
  }else{
    RSS <- makewRSS(x, y, w, ptype, ctype)
    Grad <- makewGrad(x, y, w, ptype, ctype)
    i <- 1
    GradHess <- makewGradHessBlk1(x, y, w, par, i, ptype, ctype)
    envGradHess <- environment(GradHess)
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
  TOL <- .Machine$double.eps*1000
  isncpar <- c(TRUE, TRUE, rep(c(TRUE, FALSE), K))
  
  converged <- FALSE
  for(ol in 1:maxiter){
    odd <- FALSE
    for(i in 1:(2*K)){
      j <- i+2L
      odd <- !odd
      assign("i", i, envir=envGradHess)
      
      if( odd ){
        tt <- GradHess(par)
        par[j]  <- par[j] - tt[1L]/tt[2L]
        tt <- GradHess(par)
        par[j]  <- par[j] - tt[1L]/tt[2L]
      }else{
        tt <- GradHess(par)
        tmp <- par[j] - tt[1L]/tt[2L]
        if(tmp >= 0){
          par[j] <- tmp
        }else{
          par[j] <- 0
        }
      }
    }
    tmp <- lm.fit(cbind(1, evalPolMonPartilde(par, x, ptype, ctype)), y)
    par[1:2] <- as.vector(coef(tmp))
    assign("a", par[2L], envir=envGradHess)
    grad <- Grad(par)
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
