###  Copyright (C) 2013 Berwin A. Turlach <Berwin.Turlach@gmail.com>
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
FullSOS <- function(x, y, w, deg.is.odd, type, K, par,
                    trace, plot.it, control){

  if(missing(w) || is.null(w)){
    RSS <- makeRSSSOS(x, y, type, deg.is.odd, K)
    Grad <- makeGradSOS(x, y, type, deg.is.odd, K)
    Hess <- makeGradHessSOS(x, y, type, deg.is.odd, K)
  }else{
    RSS <- makewRSSSOS(x, y, w, type, deg.is.odd, K)
    Grad <- makewGradSOS(x, y, w, type, deg.is.odd, K)
    Hess <- makewGradHessSOS(x, y, w, type, deg.is.odd, K)
  }
  
  if(plot.it)
    xgr <- seq(from=min(x), to=max(x), length=401)
  if(trace){
    cat("\n----------------------\n")
    cat("RSS:\t", RSS(par), "\n")
    cat("\n")
    print(zapsmall(rbind(par[-2L],grad=Grad(par))))
    cat("\n")
    print(zapsmall(evalCoefSOS(par, type, deg.is.odd, K)))
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
      yy <- diag(tt$Hess)
      yy[yy<1e-7] <- 1
      zz <- try(tmp <- par[-2L] + solve(tt$Hess+lam*diag(yy), -tt$Grad),
                silent=TRUE)
      if(inherits(zz, "try-error")){
        nRSS <- cRSS+1
      }else{
        nRSS <- RSS(c(tmp[1], par[2L], tmp[-1]))
      }
      
      if(nRSS <= cRSS){
        par[-2L] <- tmp
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
      print(zapsmall(rbind(par[-2L],grad=grad)))
      cat("\n")
      print(zapsmall(evalCoefSOS(par, type, deg.is.odd, K)))
    }
    if(plot.it && (ol %% plot.it == 0)){
      ygr <- evalPolSOS(par, xgr, type, deg.is.odd, K)
      lines(xgr,ygr, col="blue")
    }
    if(max(abs(grad)) < tol){
      converged <- TRUE
      break
    }
  }
  names(grad) <- names(par)
  list(par=par, grad=grad, beta=evalCoefSOS(par, type, deg.is.odd, K),
       RSS=RSS(par), niter=ol, converged=converged)
}
