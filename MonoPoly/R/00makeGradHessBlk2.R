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
makeGradHessBlk2 <- function(x, y, par, i,
                                ptype=c("Elphinstone", "EHH", "Penttila"),
                                ctype=c("c2", "cge0")){
  force(x)
  force(y)
  force(i)
  a <- par[2L]
  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  function(par){
    grhe <- evalGradHessPolBlk2(par, x, i, ptype, ctype)
    grad <- a * grhe[,c(1L,3L)]
    fit <- evalPolMonPar(par, x, ptype, ctype)
    tmp <- -2*(y-fit)
    res1 <- colSums(tmp*grad)
    res2 <- matrix(0, ncol=2, nrow=2)
    res2[1L, 1L] <- sum(2*grad[,1L]^2 + tmp*a*grhe[,2L])
    res2[2L, 2L] <- sum(2*grad[,2L]^2 + tmp*a*grhe[,4L])
    res2[1L, 2L] <- res2[2L, 1L] <- 2*sum(grad[,1L]*grad[,2L])
    
    list(Grad=res1, Hess=res2)
  }
}

makewGradHessBlk2 <- function(x, y, w, par, i,
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
    grhe <- evalGradHessPolBlk2(par, x, i, ptype, ctype)
    grad <- a * grhe[,c(1L,3L)]
    fit <- evalPolMonPar(par, x, ptype, ctype)
    tmp <- -2*w*(y-fit)
    res1 <- colSums(tmp*grad)
    res2 <- matrix(0, ncol=2, nrow=2)
    res2[1L, 1L] <- sum(2*w*grad[,1L]^2 + tmp*a*grhe[,2L])
    res2[2L, 2L] <- sum(2*w*grad[,2L]^2 + tmp*a*grhe[,4L])
    res2[1L, 2L] <- res2[2L, 1L] <- 2*sum(w*grad[,1L]*grad[,2L])
    
    list(Grad=res1, Hess=res2)
  }
}

evalGradHessPolBlk2 <- function(par, x, i,
                                ptype=c("Elphinstone", "EHH", "Penttila"),
                                ctype=c("c2", "cge0")){

  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  par <- par[-(1:2)]
  polynom.coef <- evalGradHessCoefBlk2(par, i, ptype, ctype)
  res <- matrix(0, nrow=length(x), ncol=4)
  
  for(i in 1:4){
    pc <- polynom.coef[[i]]
    order <- length(pc)
    tmp <- pc[order]*x
    for(j in rev(seq_len(order-1)))
      tmp <- (tmp + pc[j])*x
    res[,i] <- tmp
  }
  
  res
}

evalGradHessCoefBlk2 <- function(par, i,
                                 ptype=c("Elphinstone", "EHH", "Penttila"),
                                 ctype=c("c2", "cge0")){

  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  funName <- paste("evalGradHessCoefBlk2", ptype, ctype, sep="")
  do.call(funName, list(par=par, i=i))
}

evalGradHessCoefBlk2Elphinstonec2 <- function(par, i){
  
  len <- length(par)
  K   <- len/2
  if(i>K)
    stop("i is too large.")

  i <- 2*(i-1)+1

  jj <- 1:len
  res <- vector("list", 4)
  ii <- i:(i+1)
  res[[1L]] <- c(2*par[i], 2)
  res[[2L]] <- c(2)
  res[[3L]] <- c(2*par[i+1])
  res[[4L]] <- c(2)

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
    res[[3L]] <- convolve(res[[3L]], nc, type="o")
    res[[4L]] <- convolve(res[[4L]], nc, type="o")
  }

  res[[1L]] <- res[[1L]]/(1:length(res[[1L]]))
  res[[2L]] <- res[[2L]]/(1:length(res[[2L]]))
  res[[3L]] <- res[[3L]]/(1:length(res[[3L]]))
  res[[4L]] <- res[[4L]]/(1:length(res[[4L]]))

  res
}

evalGradHessCoefBlk2EHHc2 <- function(par, i){
  
  len <- length(par)
  K   <- len/2
  if(i>K)
    stop("i is too large.")

  i <- 2*(i-1)+1

  jj <- 1:len
  res <- vector("list", 4)
  ii <- i:(i+1)
  res[[1L]] <- c(0, 2, 2*par[i])
  res[[2L]] <- c(0, 0, 2)
  res[[3L]] <- c(0, 0, 2*par[i+1])
  res[[4L]] <- c(0, 0, 2)

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
    res[[3L]] <- convolve(res[[3L]], nc, type="o")
    res[[4L]] <- convolve(res[[4L]], nc, type="o")
  }

  res[[1L]] <- res[[1L]]/(1:length(res[[1L]]))
  res[[2L]] <- res[[2L]]/(1:length(res[[2L]]))
  res[[3L]] <- res[[3L]]/(1:length(res[[3L]]))
  res[[4L]] <- res[[4L]]/(1:length(res[[4L]]))

  res
}

evalGradHessCoefBlk2Penttilac2 <- function(par, i){
  
  len <- length(par)
  K   <- len/2
  if(i>K)
    stop("i is too large.")

  i <- 2*(i-1)+1

  jj <- 1:len
  res <- vector("list", 4)
  ii <- i:(i+1)
  res[[1L]] <- c(2*par[i], 2)
  res[[2L]] <- c(2)
  res[[3L]] <- c(0, 0, 2*par[i+1])
  res[[4L]] <- c(0, 0, 2)

  if(K > 1){
    jj <- setdiff(jj, ii)
    i <- 1L
    nc  <- c(par[jj[i]]^2, 2*par[jj[i]], 1+par[jj[i+1L]]^2)
    i <- i+2L
    K <- K-1
    while(K > 1){
      nc1  <- c(1+par[jj[i+1L]]^2, 2*par[jj[i]],  par[jj[i]]^2)
      nc <- convolve(nc, nc1, type="o")
      
      i <- i+2L
      K <- K-1
    }
    nc <- rev(nc)
    res[[1L]] <- convolve(res[[1L]], nc, type="o")
    res[[2L]] <- convolve(res[[2L]], nc, type="o")
    res[[3L]] <- convolve(res[[3L]], nc, type="o")
    res[[4L]] <- convolve(res[[4L]], nc, type="o")
  }

  res[[1L]] <- res[[1L]]/(1:length(res[[1L]]))
  res[[2L]] <- res[[2L]]/(1:length(res[[2L]]))
  res[[3L]] <- res[[3L]]/(1:length(res[[3L]]))
  res[[4L]] <- res[[4L]]/(1:length(res[[4L]]))

  res
}

evalGradHessCoefBlk2Elphinstonecge0 <- function(par, i){
  
  len <- length(par)
  K   <- len/2
  if(i>K)
    stop("i is too large.")

  i <- 2*(i-1)+1

  jj <- 1:len
  res <- vector("list", 4)
  ii <- i:(i+1)
  res[[1L]] <- c(2*par[i], 2)
  res[[2L]] <- c(2)
  res[[3L]] <- c(1)
  res[[4L]] <- c(0)

  if(K > 1){
    jj <- setdiff(jj, ii)
    i <- 1L
    nc  <- c(par[jj[i]]^2+par[jj[i+1L]], 2*par[jj[i]], 1)
    i <- i+2L
    K <- K-1
    while(K > 1){
      nc1  <- c(1, 2*par[jj[i]],  par[jj[i]]^2+par[jj[i+1L]])
      nc <- convolve(nc, nc1, type="o")
      
      i <- i+2L
      K <- K-1
    }
    nc <- rev(nc)
    res[[1L]] <- convolve(res[[1L]], nc, type="o")
    res[[2L]] <- convolve(res[[2L]], nc, type="o")
    res[[3L]] <- convolve(res[[3L]], nc, type="o")
  }
  
  res[[1L]] <- res[[1L]]/(1:length(res[[1L]]))
  res[[2L]] <- res[[2L]]/(1:length(res[[2L]]))
  res[[3L]] <- res[[3L]]/(1:length(res[[3L]]))

  res
}

evalGradHessCoefBlk2EHHcge0 <- function(par, i){
  
  len <- length(par)
  K   <- len/2
  if(i>K)
    stop("i is too large.")

  i <- 2*(i-1)+1

  jj <- 1:len
  res <- vector("list", 4)
  ii <- i:(i+1)
  res[[1L]] <- c(0, 2, 2*par[i])
  res[[2L]] <- c(0, 0, 2)
  res[[3L]] <- c(0, 0, 1)
  res[[4L]] <- c(0)

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
    res[[2L]] <- convolve(res[[2L]], nc, type="o")
    res[[3L]] <- convolve(res[[3L]], nc, type="o")
  }

  res[[1L]] <- res[[1L]]/(1:length(res[[1L]]))
  res[[2L]] <- res[[2L]]/(1:length(res[[2L]]))
  res[[3L]] <- res[[3L]]/(1:length(res[[3L]]))

  res
}

evalGradHessCoefBlk2Penttilacge0 <- function(par, i){
  
  len <- length(par)
  K   <- len/2
  if(i>K)
    stop("i is too large.")

  i <- 2*(i-1)+1

  jj <- 1:len
  res <- vector("list", 4)
  ii <- i:(i+1)
  res[[1L]] <- c(2*par[i], 2)
  res[[2L]] <- c(2)
  res[[3L]] <- c(0, 0, 1)
  res[[4L]] <- c(0)

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
    res[[2L]] <- convolve(res[[2L]], nc, type="o")
    res[[3L]] <- convolve(res[[3L]], nc, type="o")
  }

  res[[1L]] <- res[[1L]]/(1:length(res[[1L]]))
  res[[2L]] <- res[[2L]]/(1:length(res[[2L]]))
  res[[3L]] <- res[[3L]]/(1:length(res[[3L]]))

  res
}
