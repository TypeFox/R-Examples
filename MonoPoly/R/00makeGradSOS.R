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
makeGradSOS <- function(x, y, type, deg.is.odd, K){

  force(x)
  force(y)
  force(type)
  force(deg.is.odd)
  force(K)

  function(par){
    grad <- evalGradPolSOS(par, x, type, deg.is.odd, K)
    fit <- as.vector(evalPolSOS(par, x, type, deg.is.odd, K))
    tmp <- -2*(y-fit)
    colSums(tmp*grad)
  }
}

makewGradSOS <- function(x, y, w, type, deg.is.odd, K){

  force(x)
  force(y)
  force(w)
  force(type)
  force(deg.is.odd)
  force(K)
  
  function(par){
    grad <- evalGradPolSOS(par, x, type, deg.is.odd, K)
    fit <- as.vector(evalPolSOS(par, x, type, deg.is.odd, K))
    tmp <- -2*w*(y-fit)
    colSums(tmp*grad)
  }
}

evalGradPolSOS <- function(par, x, type, deg.is.odd, K){

  polynom.coef <- evalGradCoefSOS(par, type, deg.is.odd, K)
  res <- matrix(0, nrow=length(x), ncol=length(polynom.coef)+1)
  res[,1L] <- 1
  
  for(i in 1:length(polynom.coef)){
    pc <- polynom.coef[[i]]
    order <- length(pc)
    tmp <- pc[order]*x
    for(j in rev(seq_len(order-1)))
      tmp <- (tmp + pc[j])*x
    res[,i+1] <- tmp
  }
  
  res
}

evalGradCoefSOS <- function(par, type, deg.is.odd, K){

  a <- par[2L]
  if(type == 0){
    M <- (length(par)-2L)/2L  ## M == K+1, 2K+1==degree

    res <- vector("list",2L*M)
    res[[1]] <- par[2L + 1L:M]
    res[[M+1L]] <- par[M+2L+ 1L:M]
    if( M > 1L ){
      for(i in 2:M){
        res[[i]] <- c(0,res[[i-1L]])
        res[[M+i]] <- c(0,res[[M+i-1L]])
      }
    }
  }else if(type==1){
    if(deg.is.odd){
      M <- K+1
      res <- vector("list", M+K)
      res[[1]] <- par[2L + 1L:M]
      res[[2]] <- c(0, res[[1]])
      res[[M+1L]] <- c(0, par[M+2L + 1L:K])
      if( M > 2L ){
        for(i in 3:M){
          res[[i]] <- c(0, res[[i-1L]])
          res[[M+i-1L]] <- c(0, res[[M+i-2L]])
        }
      }
    }else{
      res <- vector("list", 2L*K)
      res[[1]] <- par[2L + 1L:K]
      res[[K+1L]] <- c(0, par[K+2L + 1L:K])
      if( K > 1L ){
        for(i in 2:K){
          res[[i]] <- c(0,res[[i-1L]])
          res[[K+i]] <- c(0,res[[K+i-1L]])
        }
      }
    }
  }else if(type==2){
    if(deg.is.odd){
      M <- K+1
      res <- vector("list", M+K)
      res[[1]] <- par[2L + 1L:M]
      res[[2]] <- c(0, res[[1]])
      res[[M+1L]] <- convolve(par[M+2L + 1L:K], c(-1, 1, 0), type="o")
      if( M > 2L ){
        for(i in 3:M){
          res[[i]] <- c(0, res[[i-1L]])
          res[[M+i-1L]] <- c(0, res[[M+i-2L]])
        }
      }
    }else{
      res <- vector("list", 2L*K)
      tmp <- par[2L + 1L:K]
      res[[1]] <- c(tmp, 0) - c(0, tmp)
      res[[K+1L]] <- c(0, par[K+2L + 1L:K])
      if( K > 1L ){
        for(i in 2:K){
          res[[i]] <- c(0,res[[i-1L]])
          res[[K+i]] <- c(0,res[[K+i-1L]])
        }
      }
    }
  }else{
    stop("How did we get here?")
  }

  for(i in seq_along(res))
    res[[i]] <- 2*a*res[[i]] / (1:length(res[[i]]))
  res
}

makeGradHessSOS <- function(x, y, type, deg.is.odd, K){

  force(x)
  force(y)
  force(type)
  force(deg.is.odd)
  force(K)
  
  function(par){
    grad <- evalGradPolSOS(par, x, type, deg.is.odd, K)
    fit <- evalPolSOS(par, x, type, deg.is.odd, K)
    resid <- as.vector(y-fit)
    res1 <- -colSums(resid*grad)

    res2 <- crossprod(grad, grad)
    res2 <- .finishHess(par, res2, resid, x, type, deg.is.odd, K)
          
    list(Grad=2*res1, Hess=2*res2)
  }
}

makewGradHessSOS <- function(x, y, w, type, deg.is.odd, K){

  force(x)
  force(y)
  force(w)
  force(type)
  force(deg.is.odd)
  force(K)
  
  function(par){
    grad <- evalGradPolSOS(par, x, type, deg.is.odd, K)
    fit <- evalPolSOS(par, x, type, deg.is.odd, K)
    resid <- w*as.vector(y-fit)
    res1 <- -colSums(resid*grad)

    res2 <- crossprod(w*grad, grad)
    res2 <- .finishHess(par, res2, x, type, deg.is.odd, K)
          
    list(Grad=2*res1, Hess=2*res2)
  }
}

.finishHess <- function(par, res2, resid, x, type, deg.is.odd, K){
  tmp <- 2*par[2L]*resid*x
  if(type==0){
    pm1 <- (NROW(res2)-1)/2 
    p <- pm1 + 1
    for(i in 2:p){
      ind <- cbind(2:i, i:2)
      st <- sum(tmp)
      res2[ind] <- res2[ind] - st
      res2[ind+pm1] <- res2[ind+pm1] - st
      tmp <- tmp*x/i*(i-1)
    }
    for(i in (p+1):(2*p-2)){
      ind <- cbind((i+2-p):p, p:(i+2-p))
      st <- sum(tmp)
      res2[ind] <- res2[ind] - st
      res2[ind+pm1] <- res2[ind+pm1] - st
      tmp <- tmp*x/i*(i-1)
    }
  }else if(type==1){
    if(deg.is.odd){
      p <- NROW(res2)/2
      pm1 <- p - 1
      st <- rep(0, 2*p-1)
      for(i in 1:(2*p-1)){
        st[i] <- sum(tmp)
        tmp <- tmp*x/(i+1)*i
      }
      for(i in 2:(p+1)){
        ind <- cbind(2:i, i:2)
        res2[ind] <- res2[ind] - st[i-1]
      }
      for(i in (p+2):(2*p)){
        ind <- cbind((i+1-p):(p+1), (p+1):(i+1-p))
        res2[ind] <- res2[ind] - st[i-1]
      }
      for(i in 2:p){
        ind <- cbind(2:i, i:2)+p
        res2[ind] <- res2[ind] - st[i]
      }
      if(p > 2){
        for(i in (p+1):(2*p-2)){
          ind <- cbind((i+2-p):p, p:(i+2-p)) + p
          res2[ind] <- res2[ind] - st[i]
        }
      }
    }else{
      pm1 <- (NROW(res2)-1)/2 
      p <- pm1 + 1
      st <- sum(tmp)
      for(i in 2:p){
        ind <- cbind(2:i, i:2)
        res2[ind] <- res2[ind] - st
        tmp <- tmp*x/i*(i-1)
        st <- sum(tmp)
        res2[ind+pm1] <- res2[ind+pm1] - st
      }
      if(pm1>1){
        for(i in (p+1):(2*p-2)){
          ind <- cbind((i+2-p):p, p:(i+2-p))
          res2[ind] <- res2[ind] - st
          tmp <- tmp*x/i*(i-1)
          st <- sum(tmp)
          res2[ind+pm1] <- res2[ind+pm1] - st
        }
      }
    }
  }else if(type==2){
    if(deg.is.odd){
      p <- NROW(res2)/2
      pm1 <- p - 1
      st <- rep(0, 2*p-1)
      for(i in 1:(2*p-1)){
        st[i] <- sum(tmp)
        tmp <- tmp*x/(i+1)*i
      }
      for(i in 2:(p+1)){
        ind <- cbind(2:i, i:2)
        res2[ind] <- res2[ind] - st[i-1]
      }
      for(i in (p+2):(2*p)){
        ind <- cbind((i+1-p):(p+1), (p+1):(i+1-p))
        res2[ind] <- res2[ind] - st[i-1]
      }
      for(i in 2:p){
        ind <- cbind(2:i, i:2)+p
        res2[ind] <- res2[ind] - (st[i] - st[i+1])
      }
      if(p > 2){
        for(i in (p+1):(2*p-2)){
          ind <- cbind((i+2-p):p, p:(i+2-p)) + p
          res2[ind] <- res2[ind] - (st[i] - st[i+1])
        }
      }
    }else{
      pm1 <- (NROW(res2)-1)/2 
      p <- pm1 + 1
      st <- rep(0, 2*p-2)
      for(i in 1:(2*p-2)){
        st[i] <- sum(tmp)
        tmp <- tmp*x/(i+1)*i
      }
      for(i in 2:p){
        ind <- cbind(2:i, i:2)
        res2[ind] <- res2[ind] - (st[i-1]-st[i])
        res2[ind+pm1] <- res2[ind+pm1] - st[i]
      }
      if(pm1>1){
        for(i in (p+1):(2*p-2)){
          ind <- cbind((i+2-p):p, p:(i+2-p))
          res2[ind] <- res2[ind] - (st[i-1]-st[i])
          res2[ind+pm1] <- res2[ind+pm1] - st[i]
        }
      }
    }
  }else{
    stop("How did we get here?")
  }
  res2
}
