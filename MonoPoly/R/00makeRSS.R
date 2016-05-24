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
makeRSS <- function(x, y, ptype=c("Elphinstone", "EHH", "Penttila"),
                    ctype=c("c2", "cge0")){
  force(x)
  force(y)
  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  function(par){
    fit <- evalPolMonPar(par, x, ptype, ctype)
    sum((y-fit)^2)
  }
}

makewRSS <- function(x, y, w, ptype=c("Elphinstone", "EHH", "Penttila"),
                    ctype=c("c2", "cge0")){
  force(x)
  force(y)
  force(w)
  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  function(par){
    fit <- evalPolMonPar(par, x, ptype, ctype)
    sum(w*(y-fit)^2)
  }
}

evalPolMonPar <- function(par, x, ptype=c("Elphinstone", "EHH", "Penttila"),
                          ctype=c("c2", "cge0")){

  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  polynom.coef <- evalCoef(par, ptype, ctype)
  order <- length(polynom.coef)
  
  pol.eval <- polynom.coef[order]
  for(i in rev(seq_len(order-1)))
    pol.eval <- pol.eval*x + polynom.coef[i]
  
  pol.eval
}

evalCoef <- function(par, ptype=c("Elphinstone", "EHH", "Penttila"),
                     ctype=c("c2", "cge0")){

  d <- par[1L]
  a <- par[2L]

  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  funName <- paste("evalCoef", ptype, ctype, sep="")
  res <- c(d, a*do.call(funName, list(par=par)))
  names(res) <- paste("beta", seq_along(res)-1, sep="")
  res
}

evalPolMonPartilde <- function(par, x,
                               ptype=c("Elphinstone", "EHH", "Penttila"),
                               ctype=c("c2", "cge0")){

  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  polynom.coef <- evalCoeftilde(par, ptype, ctype)
  order <- length(polynom.coef)
  
  pol.eval <- polynom.coef[order]*x
  for(i in rev(seq_len(order-1)))
    pol.eval <- (pol.eval + polynom.coef[i])*x
  
  pol.eval
}

evalCoeftilde <- function(par, ptype=c("Elphinstone", "EHH", "Penttila"),
                          ctype=c("c2", "cge0")){

  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  funName <- paste("evalCoef", ptype, ctype, sep="")
  do.call(funName, list(par=par))
}

evalCoefElphinstonec2 <- function(par){
  b1 <- par[3L]
  c1 <- par[4L]
  
  K <- (length(par)-2)/2
  
  integrand.coef <- c(b1^2+c1^2, 2*b1, 1)
  ii <- 5:6
  while(K > 1){
    next.coef <- c(1, 2*par[ii[1L]], par[ii[1L]]^2+par[ii[2L]]^2)
    integrand.coef <- convolve(integrand.coef, next.coef, type="o")
    
    ii <- ii + 2L
    K <- K-1
  }
  
  integrand.coef/(1:length(integrand.coef))
}

evalCoefEHHc2 <- function(par){
  b1 <- par[3L]
  c1 <- par[4L]
  
  K <- (length(par)-2)/2
  
  integrand.coef <- c(1, 2*b1, b1^2+c1^2)
  ii <- 5:6
  while(K > 1){
    next.coef <- c(par[ii[1L]]^2+par[ii[2L]]^2, 2*par[ii[1L]], 1)
    integrand.coef <- convolve(integrand.coef, next.coef, type="o")
    
    ii <- ii + 2L
    K <- K-1
  }
  
  integrand.coef/(1:length(integrand.coef))
}

evalCoefPenttilac2 <- function(par){
  b1 <- par[3L]
  c1 <- par[4L]
  
  K <- (length(par)-2)/2
  
  integrand.coef <- c(b1^2, 2*b1, 1+c1^2)
  ii <- 5:6
  while(K > 1){
    next.coef <- c(1+par[ii[2L]]^2, 2*par[ii[1L]], par[ii[1L]]^2)
    integrand.coef <- convolve(integrand.coef, next.coef, type="o")
    
    ii <- ii + 2L
    K <- K-1
  }
  
  integrand.coef/(1:length(integrand.coef))
}

evalCoefElphinstonecge0 <- function(par){
  b1 <- par[3L]
  c1 <- par[4L]
  
  K <- (length(par)-2)/2
  
  integrand.coef <- c(b1^2+c1, 2*b1, 1)
  ii <- 5:6
  while(K > 1){
    next.coef <- c(1, 2*par[ii[1L]], par[ii[1L]]^2+par[ii[2L]])
    integrand.coef <- convolve(integrand.coef, next.coef, type="o")
    
    ii <- ii + 2L
    K <- K-1
  }
  
  integrand.coef/(1:length(integrand.coef))
}

evalCoefEHHcge0 <- function(par){
  b1 <- par[3L]
  c1 <- par[4L]
  
  K <- (length(par)-2)/2
  
  integrand.coef <- c(1, 2*b1, b1^2+c1)
  ii <- 5:6
  while(K > 1){
    next.coef <- c(par[ii[1L]]^2+par[ii[2L]], 2*par[ii[1L]], 1)
    integrand.coef <- convolve(integrand.coef, next.coef, type="o")
    
    ii <- ii + 2L
    K <- K-1
  }
  
  integrand.coef/(1:length(integrand.coef))
}

evalCoefPenttilacge0 <- function(par){
  b1 <- par[3L]
  c1 <- par[4L]
  
  K <- (length(par)-2)/2
  
  integrand.coef <- c(b1^2, 2*b1, 1+c1)
  ii <- 5:6
  while(K > 1){
    next.coef <- c(1+par[ii[2L]], 2*par[ii[1L]], par[ii[1L]]^2)
    integrand.coef <- convolve(integrand.coef, next.coef, type="o")
    
    ii <- ii + 2L
    K <- K-1
  }
  
  integrand.coef/(1:length(integrand.coef))
}
