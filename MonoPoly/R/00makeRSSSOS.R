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

makeRSSSOS <- function(x, y, type, deg.is.odd, K){

  force(x)
  force(y)
  force(type)
  force(deg.is.odd)
  force(K)
  
  function(par){
    fit <- evalPolSOS(par, x, type, deg.is.odd, K)
    sum((y-fit)^2)
  }
}

makewRSSSOS <- function(x, y, w, type, deg.is.odd, K){

  force(x)
  force(y)
  force(w)
  force(type)
  force(deg.is.odd)
  force(K)
  
  function(par){
    fit <- evalPolSOS(par, x, type, deg.is.odd, K)
    sum(w*(y-fit)^2)
  }
}


evalPolSOS <- function(par, x, type, deg.is.odd, K){

  beta <- evalCoefSOS(par, type, deg.is.odd, K)
  evalPol(x, beta)
}

##
## type: 0 unbounded, 1 half open, 2 bounded
evalCoefSOS <- function(par, type, deg.is.odd, K){
  
  d <- par[1L]
  a <- par[2L]
  if(K==0){
    tmp <- c(d, a*par[3L])
  }else{
    if(type == 0){
      M <- (length(par)-2L)/2L  ## M == K+1, 2K+1 == degree
      gamma <- par[2L + 1L:M]
      delta <- par[M+2L + 1L:M]
      
      tmp <- convolve(gamma, rev(gamma), type="o")
      tmp <- tmp + convolve(delta, rev(delta), type="o")
    }else if(type == 1){
      if(deg.is.odd){
        M <- K+1
        gamma <- par[2L + 1L:M]
        delta <- par[M+2L + 1L:K]
        tmp <- convolve(gamma, rev(gamma), type="o")
        tmp <- tmp + c(0, convolve(delta, rev(delta), type="o"), 0)
      }else{
        gamma <- par[2L + 1L:K]
        delta <- par[K+2L + 1L:K]
        tmp <- c(convolve(gamma, rev(gamma), type="o"), 0)
        tmp <- tmp + c(0, convolve(delta, rev(delta), type="o"))
      }
    }else if(type == 2){
      if(deg.is.odd){
        M <- K+1
        gamma <- par[2L + 1L:M]
        delta <- par[M+2L + 1L:K]
        tmp <- convolve(gamma, rev(gamma), type="o")
        tmp <- tmp + convolve(convolve(delta, rev(delta), type="o"),
                              c(-1,1,0),type="o")
        
      }else{
        gamma <- par[2L + 1L:K]
        delta <- par[K+2L + 1L:K]
        tmp <- convolve(gamma, rev(gamma), type="o")
        tmp <- c(tmp, 0) - c(0, tmp)
        tmp <- tmp + c(0, convolve(delta, rev(delta), type="o"))
      }
    }else{
      stop("Why are we here.")
    }
    res <- c(d, a*tmp/(1:length(tmp)))
  }

  names(res) <- paste("beta", seq_along(res)-1, sep="")
  res
}
