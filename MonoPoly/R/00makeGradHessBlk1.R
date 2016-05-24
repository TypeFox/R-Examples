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
makeGradHessBlk1 <- function(x, y, par, i,
                                ptype=c("Elphinstone", "EHH", "Penttila"),
                                ctype=c("c2", "cge0")){
  force(x)
  force(y)
  force(i)
  a <- par[2L]
  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  function(par){
    grhe <- evalGradHessPolBlk1(par, x, i, ptype, ctype)
    grad <- a * grhe[,1L]
    fit <- evalPolMonPar(par, x, ptype, ctype)
    tmp <- -2*(y-fit)
    res1 <- sum(tmp*grad)
    res2 <- sum(2*grad^2 + tmp*a*grhe[,2L] )
    c(Grad=res1, Hess=res2)
  }
}

makewGradHessBlk1 <- function(x, y, w, par, i,
                                ptype=c("Elphinstone", "EHH", "Penttila"),
                                ctype=c("c2", "cge0")){
  force(x)
  force(y)
  force(w)
  force(i)
  a <- par[2L]
  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  function(par){
    grhe <- evalGradHessPolBlk1(par, x, i, ptype, ctype)
    grad <- a * grhe[,1L]
    fit <- evalPolMonPar(par, x, ptype, ctype)
    tmp <- -2*w*(y-fit)
    res1 <- sum(tmp*grad)
    res2 <- sum(2*w*grad^2 + tmp*a*grhe[,2L] )
    c(Grad=res1, Hess=res2)
  }
}

evalGradHessPolBlk1 <- function(par, x, i,
                                ptype=c("Elphinstone", "EHH", "Penttila"),
                                ctype=c("c2", "cge0")){

  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  par <- par[-(1:2)]
  polynom.coef <- evalGradHessCoefBlk1(par, i, ptype, ctype)
  res <- matrix(0, nrow=length(x), ncol=2)
  
  for(i in 1:2){
    pc <- polynom.coef[[i]]
    order <- length(pc)
    tmp <- pc[order]*x
    for(j in rev(seq_len(order-1)))
      tmp <- (tmp + pc[j])*x
    res[,i] <- tmp
  }
  
  res
}

evalGradHessCoefBlk1 <- function(par, i,
                                 ptype=c("Elphinstone", "EHH", "Penttila"),
                                 ctype=c("c2", "cge0")){

  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  funName <- paste("evalGradHessCoefBlk1", ptype, ctype, sep="")
  do.call(funName, list(par=par, i=i))
}

evalGradHessCoefBlk1Elphinstonec2 <- function(par, i){
  
  len <- length(par)
  K   <- len/2
  if(i>len)
    stop("i is too large.")

  jj <- 1:len
  res <- vector("list", 2)
  if( i%% 2){
    ## i is odd, par[i] is a b_j
    ii <- i:(i+1)
    res[[1L]] <- c(2*par[i], 2)
    res[[2L]] <- c(2)
  }else{
    ## i is even, par[i] is a c_j
    ii <- (i-1):i
    res[[1L]] <- c(2*par[i])
    res[[2L]] <- c(2)
  }

  if(K > 1){
    jj <- setdiff(jj, ii)
    i <- 1L
    nc  <- c(par[jj[i]]^2+par[jj[i+1L]]^2, 2*par[jj[i]], 1)
    i <- i+2L
    K <- K-1
    while(K > 1){
      nc1  <- c(1, 2*par[jj[i]], par[jj[i]]^2+par[jj[i+1L]]^2)
      nc <- convolve(nc, nc1, type="o")
      
      i <- i+2L
      K <- K-1
    }
    nc <- rev(nc)
    res[[1L]] <- convolve(res[[1L]], nc, type="o")
    res[[2L]] <- convolve(res[[2L]], nc, type="o")
  }

  res[[1L]] <- res[[1L]]/(1:length(res[[1L]]))
  res[[2L]] <- res[[2L]]/(1:length(res[[2L]]))

  res
}

evalGradHessCoefBlk1EHHc2 <- function(par, i){
  
  len <- length(par)
  K   <- len/2
  if(i>len)
    stop("i is too large.")

  jj <- 1:len
  res <- vector("list", 2)
  if( i%% 2){
    ## i is odd, par[i] is a b_j
    ii <- i:(i+1)
    res[[1L]] <- c(0, 2, 2*par[i])
    res[[2L]] <- c(0, 0, 2)
  }else{
    ## i is even, par[i] is a c_j
    ii <- (i-1):i
    res[[1L]] <- c(0, 0, 2*par[i])
    res[[2L]] <- c(0, 0, 2)
  }

  if(K > 1){
    jj <- setdiff(jj, ii)
    i <- 1L
    nc  <- c(1, 2*par[jj[i]], par[jj[i]]^2+par[jj[i+1L]]^2)
    i <- i+2L
    K <- K-1
    while(K > 1){
      nc1  <- c(par[jj[i]]^2+par[jj[i+1L]]^2, 2*par[jj[i]], 1)
      nc <- convolve(nc, nc1, type="o")
      
      i <- i+2L
      K <- K-1
    }
    nc <- rev(nc)
    res[[1L]] <- convolve(res[[1L]], nc, type="o")
    res[[2L]] <- convolve(res[[2L]], nc, type="o")
  }

  res[[1L]] <- res[[1L]]/(1:length(res[[1L]]))
  res[[2L]] <- res[[2L]]/(1:length(res[[2L]]))

  res
}

evalGradHessCoefBlk1Penttilac2 <- function(par, i){
  
  len <- length(par)
  K   <- len/2
  if(i>len)
    stop("i is too large.")

  jj <- 1:len
  res <- vector("list", 2)
  if( i%% 2){
    ## i is odd, par[i] is a b_j
    ii <- i:(i+1)
    res[[1L]] <- c(2*par[i], 2)
    res[[2L]] <- c(2)
  }else{
    ## i is even, par[i] is a c_j
    ii <- (i-1):i
    res[[1L]] <- c(0, 0, 2*par[i])
    res[[2L]] <- c(0, 0, 2)
  }

  if(K > 1){
    jj <- setdiff(jj, ii)
    i <- 1L
    nc  <- c(par[jj[i]]^2, 2*par[jj[i]], 1+par[jj[i+1L]]^2)
    i <- i+2L
    K <- K-1
    while(K > 1){
      nc1  <- c(1+par[jj[i+1L]]^2, 2*par[jj[i]], par[jj[i]]^2)
      nc <- convolve(nc, nc1, type="o")
      
      i <- i+2L
      K <- K-1
    }
    nc <- rev(nc)
    res[[1L]] <- convolve(res[[1L]], nc, type="o")
    res[[2L]] <- convolve(res[[2L]], nc, type="o")
  }

  res[[1L]] <- res[[1L]]/(1:length(res[[1L]]))
  res[[2L]] <- res[[2L]]/(1:length(res[[2L]]))

  res
}

evalGradHessCoefBlk1Elphinstonecge0 <- function(par, i){
  
  len <- length(par)
  K   <- len/2
  if(i>len)
    stop("i is too large.")

  jj <- 1:len
  res <- vector("list", 2)
  if( i%% 2){
    ## i is odd, par[i] is a b_j
    odd <- TRUE
    ii <- i:(i+1)
    res[[1L]] <- c(2*par[i], 2)
    res[[2L]] <- c(2)
  }else{
    ## i is even, par[i] is a c_j
    odd <- FALSE
    ii <- (i-1):i
    res[[1L]] <- c(1)
    res[[2L]] <- c(0)
  }

  if(K > 1){
    jj <- setdiff(jj, ii)
    i <- 1L
    nc  <- c(par[jj[i]]^2+par[jj[i+1L]], 2*par[jj[i]], 1)
    i <- i+2L
    K <- K-1
    while(K > 1){
      nc1  <- c(1, 2*par[jj[i]], par[jj[i]]^2+par[jj[i+1L]])
      nc <- convolve(nc, nc1, type="o")
      
      i <- i+2L
      K <- K-1
    }
    nc <- rev(nc)
    res[[1L]] <- convolve(res[[1L]], nc, type="o")
    if(odd)
      res[[2L]] <- convolve(res[[2L]], nc, type="o")
  }

  res[[1L]] <- res[[1L]]/(1:length(res[[1L]]))
  res[[2L]] <- res[[2L]]/(1:length(res[[2L]]))

  res
}

evalGradHessCoefBlk1EHHcge0 <- function(par, i){
  
  len <- length(par)
  K   <- len/2
  if(i>len)
    stop("i is too large.")

  jj <- 1:len
  res <- vector("list", 2)
  if( i%% 2){
    ## i is odd, par[i] is a b_j
    odd <- TRUE
    ii <- i:(i+1)
    res[[1L]] <- c(0, 2, 2*par[i])
    res[[2L]] <- c(0, 0, 2)
  }else{
    ## i is even, par[i] is a c_j
    odd <- FALSE
    ii <- (i-1):i
    res[[1L]] <- c(0, 0, 1)
    res[[2L]] <- c(0)
  }

  if(K > 1){
    jj <- setdiff(jj, ii)
    i <- 1L
    nc  <- c(1, 2*par[jj[i]], par[jj[i]]^2+par[jj[i+1L]])
    i <- i+2L
    K <- K-1
    while(K > 1){
      nc1  <- c(par[jj[i]]^2+par[jj[i+1L]], 2*par[jj[i]], 1)
      nc <- convolve(nc, nc1, type="o")
      
      i <- i+2L
      K <- K-1
    }
    nc <- rev(nc)
    res[[1L]] <- convolve(res[[1L]], nc, type="o")
    if(odd)
      res[[2L]] <- convolve(res[[2L]], nc, type="o")
  }
    
  res[[1L]] <- res[[1L]]/(1:length(res[[1L]]))
  res[[2L]] <- res[[2L]]/(1:length(res[[2L]]))

  res
}

evalGradHessCoefBlk1Penttilacge0 <- function(par, i){
  
  len <- length(par)
  K   <- len/2
  if(i>len)
    stop("i is too large.")

  jj <- 1:len
  res <- vector("list", 2)
  if( i%% 2){
    ## i is odd, par[i] is a b_j
    odd <- TRUE
    ii <- i:(i+1)
    res[[1L]] <- c(2*par[i], 2)
    res[[2L]] <- c(2)
  }else{
    ## i is even, par[i] is a c_j
    odd <- FALSE
    ii <- (i-1):i
    res[[1L]] <- c(0, 0, 1)
    res[[2L]] <- c(0)
  }

  if(K > 1){
    jj <- setdiff(jj, ii)
    i <- 1L
    nc  <- c(par[jj[i]]^2, 2*par[jj[i]], 1+par[jj[i+1L]])
    i <- i+2L
    K <- K-1
    while(K > 1){
      nc1  <- c(1+par[jj[i+1L]], 2*par[jj[i]], par[jj[i]]^2)
      nc <- convolve(nc, nc1, type="o")
      
      i <- i+2L
      K <- K-1
    }
    nc <- rev(nc)
    res[[1L]] <- convolve(res[[1L]], nc, type="o")
    if(odd)
      res[[2L]] <- convolve(res[[2L]], nc, type="o")
  }

  res[[1L]] <- res[[1L]]/(1:length(res[[1L]]))
  res[[2L]] <- res[[2L]]/(1:length(res[[2L]]))

  res
}
