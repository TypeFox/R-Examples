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
makeGrad <- function(x, y, ptype=c("Elphinstone", "EHH", "Penttila"),
                     ctype=c("c2", "cge0")){
  force(x)
  force(y)
  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  function(par){
    grad <- evalGradPol(par, x, ptype, ctype)
    fit <- par[1L] + par[2L] * grad[,2L]
    tmp <- -2*(y-fit)
    colSums(tmp*grad)
  }
}

makewGrad <- function(x, y, w, ptype=c("Elphinstone", "EHH", "Penttila"),
                     ctype=c("c2", "cge0")){
  force(x)
  force(y)
  force(w)
  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  function(par){
    grad <- evalGradPol(par, x, ptype, ctype)
    fit <- par[1L] + par[2L] * grad[,2L]
    tmp <- -2*w*(y-fit)
    colSums(tmp*grad)
  }
}

evalGradPol <- function(par, x, ptype=c("Elphinstone", "EHH", "Penttila"),
                        ctype=c("c2", "cge0")){

  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  polynom.coef <- evalGradCoef(par, ptype, ctype)
  res <- matrix(0, nrow=length(x), ncol=length(par))
  res[,1L] <- 1
  
  for(i in 1:length(polynom.coef)){
    pc <- polynom.coef[[i]]
    order <- length(pc)
    tmp <- pc[order]*x
    for(j in rev(seq_len(order-1)))
      tmp <- (tmp + pc[j])*x
    res[,i+1] <- tmp
  }
  res[,-(1:2)] <- res[,-(1:2)] * par[2L]
  
  res
}

evalGradCoef <- function(par, ptype=c("Elphinstone", "EHH", "Penttila"),
                         ctype=c("c2", "cge0")){

  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  funName <- paste("evalGradCoef", ptype, ctype, sep="")
  do.call(funName, list(par=par))
}

evalGradCoefElphinstonec2 <- function(par){
  b1 <- par[3L]
  c1 <- par[4L]
  
  K <- (length(par)-2)/2
  order <- 2*K+1

  integrand.coef <- c(b1^2+c1^2, 2*b1, 1)
  res <- replicate(order, integrand.coef, simplify=FALSE)
  res[[2L]] <- c(2*b1, 2)
  res[[3L]] <- 2*c1
  ii <- 5:6
  jj <- 1:order
  while(K > 1){
    nc  <- c(1, 2*par[ii[1L]], par[ii[1L]]^2+par[ii[2L]]^2)
    nc1 <- c(2, 2*par[ii[1L]])
    nc2 <- 2*par[ii[2L]]
    ll <- ii-1L
    kk <- setdiff(jj, ll)
    for(i in kk){
      res[[i]] <- convolve(res[[i]], nc, type="o")
    }
    res[[ll[1L]]] <- convolve(res[[ll[1L]]], nc1, type="o")
    res[[ll[2L]]] <- res[[ll[2L]]] * nc2
    
    ii <- ii + 2L
    K <- K-1
  }

  for(j in jj)
    res[[j]] <- res[[j]]/(1:length(res[[j]]))

  res
}

evalGradCoefEHHc2 <- function(par){
  b1 <- par[3L]
  c1 <- par[4L]
  
  K <- (length(par)-2)/2
  order <- 2*K+1

  integrand.coef <- c(1, 2*b1, b1^2+c1^2)
  res <- replicate(order, integrand.coef, simplify=FALSE)
  res[[2L]] <- c(0, 2, 2*b1)
  res[[3L]] <- c(0, 0, 2*c1)
  ii <- 5:6
  jj <- 1:order
  while(K > 1){
    nc  <- c(par[ii[1L]]^2+par[ii[2L]]^2, 2*par[ii[1L]], 1)
    nc1 <- c(2*par[ii[1L]], 2, 0)
    nc2 <- c(2*par[ii[2L]], 0, 0)
    ll <- ii-1L
    kk <- setdiff(jj, ll)
    for(i in kk){
      res[[i]] <- convolve(res[[i]], nc, type="o")
    }
    res[[ll[1L]]] <- convolve(res[[ll[1L]]], nc1, type="o")
    res[[ll[2L]]] <- convolve(res[[ll[2L]]], nc2, type="o")
    
    ii <- ii + 2L
    K <- K-1
  }

  for(j in jj)
    res[[j]] <- res[[j]]/(1:length(res[[j]]))

  res
}

evalGradCoefPenttilac2 <- function(par){
  b1 <- par[3L]
  c1 <- par[4L]
  
  K <- (length(par)-2)/2
  order <- 2*K+1

  integrand.coef <- c(b1^2, 2*b1, 1+c1^2)
  res <- replicate(order, integrand.coef, simplify=FALSE)
  res[[2L]] <- c(2*b1, 2)
  res[[3L]] <- c(0, 0, 2*c1)
  ii <- 5:6
  jj <- 1:order
  while(K > 1){
    nc  <- c(1+par[ii[2L]]^2, 2*par[ii[1L]], par[ii[1L]]^2)
    nc1 <- c(2, 2*par[ii[1L]])
    nc2 <- c(2*par[ii[2L]], 0, 0)
    ll <- ii-1L
    kk <- setdiff(jj, ll)
    for(i in kk){
      res[[i]] <- convolve(res[[i]], nc, type="o")
    }
    res[[ll[1L]]] <- convolve(res[[ll[1L]]], nc1, type="o")
    res[[ll[2L]]] <- convolve(res[[ll[2L]]], nc2, type="o")
    
    ii <- ii + 2L
    K <- K-1
  }

  for(j in jj)
    res[[j]] <- res[[j]]/(1:length(res[[j]]))

  res
}

evalGradCoefElphinstonecge0 <- function(par){
  b1 <- par[3L]
  c1 <- par[4L]
  
  K <- (length(par)-2)/2
  order <- 2*K+1

  integrand.coef <- c(b1^2+c1, 2*b1, 1)
  res <- replicate(order, integrand.coef, simplify=FALSE)
  res[[2L]] <- c(2*b1, 2)
  res[[3L]] <- 1
  ii <- 5:6
  jj <- 1:order
  while(K > 1){
    nc  <- c(1, 2*par[ii[1L]], par[ii[1L]]^2+par[ii[2L]])
    nc1 <- c(2, 2*par[ii[1L]])
    ll <- ii-1L
    kk <- setdiff(jj, ll)
    for(i in kk){
      res[[i]] <- convolve(res[[i]], nc, type="o")
    }
    res[[ll[1L]]] <- convolve(res[[ll[1L]]], nc1, type="o")
    
    ii <- ii + 2L
    K <- K-1
  }

  for(j in jj)
    res[[j]] <- res[[j]]/(1:length(res[[j]]))

  res
}

evalGradCoefEHHcge0 <- function(par){
  b1 <- par[3L]
  c1 <- par[4L]
  
  K <- (length(par)-2)/2
  order <- 2*K+1

  integrand.coef <- c(1, 2*b1, b1^2+c1)
  res <- replicate(order, integrand.coef, simplify=FALSE)
  res[[2L]] <- c(0, 2, 2*b1)
  res[[3L]] <- c(0, 0, 1)
  ii <- 5:6
  jj <- 1:order
  while(K > 1){
    nc  <- c(par[ii[1L]]^2+par[ii[2L]], 2*par[ii[1L]], 1)
    nc1 <- c(2*par[ii[1L]], 2, 0)
    nc2 <- c(1, 0, 0)
    ll <- ii-1L
    kk <- setdiff(jj, ll)
    for(i in kk){
      res[[i]] <- convolve(res[[i]], nc, type="o")
    }
    res[[ll[1L]]] <- convolve(res[[ll[1L]]], nc1, type="o")
    res[[ll[2L]]] <- convolve(res[[ll[2L]]], nc2, type="o")
    
    ii <- ii + 2L
    K <- K-1
  }

  for(j in jj)
    res[[j]] <- res[[j]]/(1:length(res[[j]]))

  res
}

evalGradCoefPenttilacge0 <- function(par){
  b1 <- par[3L]
  c1 <- par[4L]
  
  K <- (length(par)-2)/2
  order <- 2*K+1

  integrand.coef <- c(b1^2, 2*b1, 1+c1)
  res <- replicate(order, integrand.coef, simplify=FALSE)
  res[[2L]] <- c(2*b1, 2)
  res[[3L]] <- c(0, 0, 1)
  ii <- 5:6
  jj <- 1:order
  while(K > 1){
    nc  <- c(1+par[ii[2L]], 2*par[ii[1L]], par[ii[1L]]^2)
    nc1 <- c(2, 2*par[ii[1L]])
    nc2 <- c(1, 0, 0)
    ll <- ii-1L
    kk <- setdiff(jj, ll)
    for(i in kk){
      res[[i]] <- convolve(res[[i]], nc, type="o")
    }
    res[[ll[1L]]] <- convolve(res[[ll[1L]]], nc1, type="o")
    res[[ll[2L]]] <- convolve(res[[ll[2L]]], nc2, type="o")
    
    ii <- ii + 2L
    K <- K-1
  }

  for(j in jj)
    res[[j]] <- res[[j]]/(1:length(res[[j]]))

  res
}
