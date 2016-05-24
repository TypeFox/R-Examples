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
makeGradHessFull <- function(x, y, ptype=c("Elphinstone", "EHH", "Penttila"),
                             ctype=c("c2", "cge0")){
  force(x)
  force(y)
  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  function(par){
    grad <- evalGradPol(par, x, ptype, ctype)
    fit <- par[1L] + par[2L] * grad[,2L]
    resid <- y-fit
    res1 <- -colSums(resid*grad)

    res2 <- crossprod(grad, grad)
    p <- NROW(res2)
    tmp <- res1[-(1:2)]/par[2L]
    res2[2L, 3:p] <- res2[3:p, 2L] <- res2[2L, 3:p] + tmp
    for(i in 3:p){
      hess <- evalHessPol(par, x, i, p, ptype, ctype)
      tmp <- -colSums(resid*hess) * par[2L]
      if(i%%2)
        tmp <- c(tmp[1], 0, tmp[-1])
      res2[i, i:p] <- res2[i:p, i] <- res2[i, i:p] + tmp
    }

    list(Grad=2*res1, Hess=2*res2)
  }
}

makewGradHessFull <- function(x, y, w,
                              ptype=c("Elphinstone", "EHH", "Penttila"),
                              ctype=c("c2", "cge0")){
  force(x)
  force(y)
  force(w)
  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  function(par){
    grad <- evalGradPol(par, x, ptype, ctype)
    fit <- par[1L] + par[2L] * grad[,2L]
    resid <- w*(y-fit)
    res1 <- -colSums(resid*grad)

    res2 <- crossprod(w*grad, grad)
    p <- NROW(res2)
    tmp <- res1[-(1:2)]/par[2L]
    res2[2L, 3:p] <- res2[3:p, 2L] <- res2[2L, 3:p] + tmp
    for(i in 3:p){
      hess <- evalHessPol(par, x, i, p, ptype, ctype)
      tmp <- -colSums(resid*hess) * par[2L]
      if(i%%2)
        tmp <- c(tmp[1], 0, tmp[-1])
      res2[i, i:p] <- res2[i:p, i] <- res2[i, i:p] + tmp
    }

    list(Grad=2*res1, Hess=2*res2)
  }
}

evalHessPol <- function(par, x, i, p,
                        ptype=c("Elphinstone", "EHH", "Penttila"),
                        ctype=c("c2", "cge0")){

  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  par <- par[-(1:2)]
  i <- i-2
  p <- p-2
  polynom.coef <- evalHessCoef(par, i, p, ptype, ctype)
  if(i%%2){
    res <- matrix(0, nrow=length(x), ncol=p-i)
  }else{
    res <- matrix(0, nrow=length(x), ncol=p-i+1)
  }
  for(i in 1:NCOL(res)){
    pc <- polynom.coef[[i]]
    order <- length(pc)
    tmp <- pc[order]*x
    for(j in rev(seq_len(order-1)))
      tmp <- (tmp + pc[j])*x
    res[,i] <- tmp
  }
  
  res
}

evalHessCoef <- function(par, i, p,
                         ptype=c("Elphinstone", "EHH", "Penttila"),
                         ctype=c("c2", "cge0")){

  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  funName <- paste("evalHessCoef", ptype, ctype, sep="")
  do.call(funName, list(par=par, i=i, len=p))
}

evalHessCoefElphinstonec2 <- function(par, i, len){

  K   <- len/2
  if(i>len)
    stop("i is too large.")

  if( i%% 2){
    ## i is odd, par[i] is a b_j
    ii <- i:(i+1)
    res <- replicate(len-i+1, c(2*par[i], 2), simplify=FALSE)
    res[[1L]] <- c(2)
  }else{
    ## i is even, par[i] is a c_j
    ii <- (i-1):i
    res <- replicate(len-i+1, c(2*par[i]), simplify=FALSE)
    res[[1L]] <- c(2)
  }

  rr <- seq_along(res)
  jj <- setdiff(1:len, ii)

  jj1 <- jj[jj<i]
  if(length(jj1)>1){
    K1 <- length(jj1)/2
    k <- 1L
    nc  <- c(par[jj1[k]]^2+par[jj1[k+1L]]^2, 2*par[jj1[k]], 1)
    k <- k+2L
    K1 <- K1-1
    while(K1 > 0){
      nc1  <- c(1, 2*par[jj1[k]], par[jj1[k]]^2+par[jj1[k+1L]]^2)
      nc <- convolve(nc, nc1, type="o")
      k <- k+2L
      K1 <- K1-1L
    }
    nc <- rev(nc)
    for(k in rr)
      res[[k]] <- convolve(res[[k]], nc, type="o")
  }

  jj1 <- jj[jj>i]
  if(length(jj1)>1){
    K1 <- length(jj1)/2
    k <- 1L
    while(K1 > 0){
      nc <- c(1, 2*par[jj1[k]], par[jj1[k]]^2+par[jj1[k+1L]]^2)
      nc1 <- c(0, 2, 2*par[jj1[k]])
      nc2 <- c(2*par[jj1[k+1L]])

      res[[k+1L]] <- convolve(res[[k+1L]], nc1, type="o")
      res[[k+2L]] <- convolve(res[[k+2L]], nc2, type="o")
      ss <- rr[-c(k+1L, k+2L)]
      for(i in ss)
        res[[i]] <- convolve(res[[i]], nc, type="o")
      
      k <- k+2L
      K1 <- K1-1L
    }
  }
  
  for(i in rr)
    res[[i]] <- res[[i]]/(1:length(res[[i]]))

  res
}

evalHessCoefEHHc2 <- function(par, i, len){

  K   <- len/2
  if(i>len)
    stop("i is too large.")

  if( i%% 2){
    ## i is odd, par[i] is a b_j
    ii <- i:(i+1)
    res <- replicate(len-i+1, c(0, 2, 2*par[i]), simplify=FALSE)
    res[[1L]] <- c(0, 0, 2)
  }else{
    ## i is even, par[i] is a c_j
    ii <- (i-1):i
    res <- replicate(len-i+1, c(0, 0, 2*par[i]), simplify=FALSE)
    res[[1L]] <- c(0, 0, 2)
  }

  rr <- seq_along(res)
  jj <- setdiff(1:len, ii)

  jj1 <- jj[jj<i]
  if(length(jj1)>1){
    K1 <- length(jj1)/2
    k <- 1L
    nc  <- c(1, 2*par[jj1[k]], par[jj1[k]]^2+par[jj1[k+1L]]^2)
    k <- k+2L
    K1 <- K1-1
    while(K1 > 0){
      nc1  <- c(par[jj1[k]]^2+par[jj1[k+1L]]^2, 2*par[jj1[k]], 1)
      nc <- convolve(nc, nc1, type="o")
      k <- k+2L
      K1 <- K1-1L
    }
    nc <- rev(nc)
    for(k in rr)
      res[[k]] <- convolve(res[[k]], nc, type="o")
  }

  jj1 <- jj[jj>i]
  if(length(jj1)>1){
    K1 <- length(jj1)/2
    k <- 1L
    while(K1 > 0){
      nc <- c(par[jj1[k]]^2+par[jj1[k+1L]]^2, 2*par[jj1[k]], 1)
      nc1 <- c(2*par[jj1[k]], 2, 0)
      nc2 <- c(2*par[jj1[k+1L]], 0, 0)

      res[[k+1L]] <- convolve(res[[k+1L]], nc1, type="o")
      res[[k+2L]] <- convolve(res[[k+2L]], nc2, type="o")
      ss <- rr[-c(k+1L, k+2L)]
      for(i in ss)
        res[[i]] <- convolve(res[[i]], nc, type="o")
      
      k <- k+2L
      K1 <- K1-1L
    }
  }
  
  for(i in rr)
    res[[i]] <- res[[i]]/(1:length(res[[i]]))

  res
}

evalHessCoefPenttilac2 <- function(par, i, len){

  K   <- len/2
  if(i>len)
    stop("i is too large.")

  if( i%% 2){
    ## i is odd, par[i] is a b_j
    ii <- i:(i+1)
    res <- replicate(len-i+1, c(2*par[i], 2), simplify=FALSE)
    res[[1L]] <- c(2)
  }else{
    ## i is even, par[i] is a c_j
    ii <- (i-1):i
    res <- replicate(len-i+1, c(0, 0, 2*par[i]), simplify=FALSE)
    res[[1L]] <- c(0, 0, 2)
  }

  rr <- seq_along(res)
  jj <- setdiff(1:len, ii)

  jj1 <- jj[jj<i]
  if(length(jj1)>1){
    K1 <- length(jj1)/2
    k <- 1L
    nc  <- c(par[jj1[k]]^2, 2*par[jj1[k]], 1+par[jj1[k+1L]]^2)
    k <- k+2L
    K1 <- K1-1
    while(K1 > 0){
      nc1  <- c(1+par[jj1[k+1L]]^2, 2*par[jj1[k]], par[jj1[k]]^2)
      nc <- convolve(nc, nc1, type="o")
      k <- k+2L
      K1 <- K1-1L
    }
    nc <- rev(nc)
    for(k in rr)
      res[[k]] <- convolve(res[[k]], nc, type="o")
  }

  jj1 <- jj[jj>i]
  if(length(jj1)>1){
    K1 <- length(jj1)/2
    k <- 1L
    while(K1 > 0){
      nc <- c(1+par[jj1[k+1L]]^2, 2*par[jj1[k]], par[jj1[k]]^2)
      nc1 <- c(0, 2, 2*par[jj1[k]])
      nc2 <- c(2*par[jj1[k+1L]], 0, 0)

      res[[k+1L]] <- convolve(res[[k+1L]], nc1, type="o")
      res[[k+2L]] <- convolve(res[[k+2L]], nc2, type="o")
      ss <- rr[-c(k+1L, k+2L)]
      for(i in ss)
        res[[i]] <- convolve(res[[i]], nc, type="o")
      
      k <- k+2L
      K1 <- K1-1L
    }
  }
  
  for(i in rr)
    res[[i]] <- res[[i]]/(1:length(res[[i]]))

  res
}

evalHessCoefElphinstonecge0 <- function(par, i, len){

  K   <- len/2
  if(i>len)
    stop("i is too large.")

  if( i%% 2){
    ## i is odd, par[i] is a b_j
    ii <- i:(i+1)
    res <- replicate(len-i+1, c(2*par[i], 2), simplify=FALSE)
    res[[1L]] <- c(2)
  }else{
    ## i is even, par[i] is a c_j
    ii <- (i-1):i
    res <- replicate(len-i+1, c(1), simplify=FALSE)
    res[[1L]] <- c(0)
  }

  rr <- seq_along(res)
  jj <- setdiff(1:len, ii)

  jj1 <- jj[jj<i]
  if(length(jj1)>1){
    K1 <- length(jj1)/2
    k <- 1L
    nc  <- c(par[jj1[k]]^2+par[jj1[k+1L]], 2*par[jj1[k]], 1)
    k <- k+2L
    K1 <- K1-1
    while(K1 > 0){
      nc1  <- c(1, 2*par[jj1[k]], par[jj1[k]]^2+par[jj1[k+1L]])
      nc <- convolve(nc, nc1, type="o")
      k <- k+2L
      K1 <- K1-1L
    }
    nc <- rev(nc)
    for(k in rr)
      res[[k]] <- convolve(res[[k]], nc, type="o")
  }

  jj1 <- jj[jj>i]
  if(length(jj1)>1){
    K1 <- length(jj1)/2
    k <- 1L
    while(K1 > 0){
      nc <- c(1, 2*par[jj1[k]], par[jj1[k]]^2+par[jj1[k+1L]])
      nc1 <- c(0, 2, 2*par[jj1[k]])
      nc2 <- c(1)

      res[[k+1L]] <- convolve(res[[k+1L]], nc1, type="o")
      res[[k+2L]] <- convolve(res[[k+2L]], nc2, type="o")
      ss <- rr[-c(k+1L, k+2L)]
      for(i in ss)
        res[[i]] <- convolve(res[[i]], nc, type="o")
      
      k <- k+2L
      K1 <- K1-1L
    }
  }
  
  for(i in rr)
    res[[i]] <- res[[i]]/(1:length(res[[i]]))

  res
}

evalHessCoefEHHcge0 <- function(par, i, len){

  K   <- len/2
  if(i>len)
    stop("i is too large.")

  if( i%% 2){
    ## i is odd, par[i] is a b_j
    ii <- i:(i+1)
    res <- replicate(len-i+1, c(0, 2, 2*par[i]), simplify=FALSE)
    res[[1L]] <- c(0, 0, 2)
  }else{
    ## i is even, par[i] is a c_j
    ii <- (i-1):i
    res <- replicate(len-i+1, c(0, 0, 1), simplify=FALSE)
    res[[1L]] <- c(0)
  }

  rr <- seq_along(res)
  jj <- setdiff(1:len, ii)

  jj1 <- jj[jj<i]
  if(length(jj1)>1){
    K1 <- length(jj1)/2
    k <- 1L
    nc  <- c(1, 2*par[jj1[k]], par[jj1[k]]^2+par[jj1[k+1L]])
    k <- k+2L
    K1 <- K1-1
    while(K1 > 0){
      nc1  <- c(par[jj1[k]]^2+par[jj1[k+1L]], 2*par[jj1[k]], 1)
      nc <- convolve(nc, nc1, type="o")
      k <- k+2L
      K1 <- K1-1L
    }
    nc <- rev(nc)
    for(k in rr)
      res[[k]] <- convolve(res[[k]], nc, type="o")
  }

  jj1 <- jj[jj>i]
  if(length(jj1)>1){
    K1 <- length(jj1)/2
    k <- 1L
    while(K1 > 0){
      nc <- c(par[jj1[k]]^2+par[jj1[k+1L]], 2*par[jj1[k]], 1)
      nc1 <- c(2*par[jj1[k]], 2, 0)
      nc2 <- c(1, 0, 0)

      res[[k+1L]] <- convolve(res[[k+1L]], nc1, type="o")
      res[[k+2L]] <- convolve(res[[k+2L]], nc2, type="o")
      ss <- rr[-c(k+1L, k+2L)]
      for(i in ss)
        res[[i]] <- convolve(res[[i]], nc, type="o")
      
      k <- k+2L
      K1 <- K1-1L
    }
  }
  
  for(i in rr)
    res[[i]] <- res[[i]]/(1:length(res[[i]]))

  res
}

evalHessCoefPenttilacge0 <- function(par, i, len){

  K   <- len/2
  if(i>len)
    stop("i is too large.")

  if( i%% 2){
    ## i is odd, par[i] is a b_j
    ii <- i:(i+1)
    res <- replicate(len-i+1, c(2*par[i], 2), simplify=FALSE)
    res[[1L]] <- c(2)
  }else{
    ## i is even, par[i] is a c_j
    ii <- (i-1):i
    res <- replicate(len-i+1, c(0, 0, 1), simplify=FALSE)
    res[[1L]] <- c(0)
  }

  rr <- seq_along(res)
  jj <- setdiff(1:len, ii)

  jj1 <- jj[jj<i]
  if(length(jj1)>1){
    K1 <- length(jj1)/2
    k <- 1L
    nc  <- c(par[jj1[k]]^2, 2*par[jj1[k]], 1+par[jj1[k+1L]])
    k <- k+2L
    K1 <- K1-1
    while(K1 > 0){
      nc1  <- c(1+par[jj1[k+1L]], 2*par[jj1[k]], par[jj1[k]]^2)
      nc <- convolve(nc, nc1, type="o")
      k <- k+2L
      K1 <- K1-1L
    }
    nc <- rev(nc)
    for(k in rr)
      res[[k]] <- convolve(res[[k]], nc, type="o")
  }

  jj1 <- jj[jj>i]
  if(length(jj1)>1){
    K1 <- length(jj1)/2
    k <- 1L
    while(K1 > 0){
      nc <- c(1+par[jj1[k+1L]], 2*par[jj1[k]], par[jj1[k]]^2)
      nc1 <- c(0, 2, 2*par[jj1[k]])
      nc2 <- c(1, 0, 0)

      res[[k+1L]] <- convolve(res[[k+1L]], nc1, type="o")
      res[[k+2L]] <- convolve(res[[k+2L]], nc2, type="o")
      ss <- rr[-c(k+1L, k+2L)]
      for(i in ss)
        res[[i]] <- convolve(res[[i]], nc, type="o")
      
      k <- k+2L
      K1 <- K1-1L
    }
  }
  
  for(i in rr)
    res[[i]] <- res[[i]]/(1:length(res[[i]]))

  res
}
