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

SOSpol.fit <- function(x, y, w=NULL, deg.is.odd, K, start,
                       a, b, monotone=c("increasing", "decreasing"),
                       trace=FALSE, plot.it=FALSE,
                       type,
                       control=monpol.control()){

  x <- switch(type+1,
              ## type == 0, rescale to [-1,1]
              {minx <- min(x)
               sclx <- max(x)-minx
               2*(x-minx)/sclx - 1},
              ## type == 1, rescale [a,Inf) to [0,Inf)
              {minx <- a
               sclx <- max(x)-minx
               (x-a)/sclx},
              ## type == 2, rescale [a,b] to [0,1]
              {minx <- a
               sclx <- b-minx
               (x-a)/sclx})
  if(is.null(x))
    stop("Why is rescaled x null?")

  miny <- min(y)
  scly <- max(y)-miny
  y <- 2*(y-miny)/scly - 1
  
  if(K==0){
    if(is.null(w)){ 
      tmp <- lm.fit(cbind(1, x), y)
    }else{
      tmp <- lm.wfit(cbind(1,x), y, w)
    }
    beta <- coef(tmp)
    par <- rep(0, 4)
    par[1L] <- beta[1L]
    par[2L] <- ifelse(beta[2L]>0, 1, -1)
    par[3L] <- abs(beta[2L])
    names(par) <- c("d", "a", "beta10", "beta20")
    names(beta) <- c("beta0", "beta1")
    res <- list(par=par, grad=rep(0,4), beta=beta,
                niter=0, converged=TRUE,
                type=type)
    if(is.null(w)){
      res$RSS <- sum(tmp$residuals^2)
    }else{
      res$RSS <- sum(w*tmp$residuals^2)
    }
    if(plot.it){
      plot(x, y, xlab="x", ylab="y")
      abline(beta, col="green")
    }
  }else{
    
    if(missing(start)){
      if(type==0){
        cfm <- coef(lm.fit(cbind(1,x,x^3),y))
        par <- rep(0, 2*K+4L)
        par[1L] <- cfm[1L]
        par[2L] <- ifelse(cor(x,y)>0, 1, -1)
        par[3L] <- sqrt(abs(cfm[2L]))
        par[K+5L] <- sqrt(3*abs(cfm[3L]))
      }else if(type==1){
        ind <- x>=0
        xx <- x[ind]
        yy <- y[ind]
        if(deg.is.odd){
          cfm <- coef(lm.fit(cbind(1,xx,xx^3),yy))
        }else{
          cfm <- coef(lm.fit(cbind(1,xx,xx^2),yy))
        }
        par <- rep(0, 2*K+2L+deg.is.odd)
        par[1L] <- cfm[1L]
        par[2L] <- ifelse(cor(xx,yy)>0, 1, -1)
        par[3L] <- sqrt(abs(cfm[2L]))
        if(deg.is.odd){
          par[4L] <- -sqrt(3*abs(cfm[3L]))
          par[K+4L] <- sqrt(2*abs(par[3L]*par[4L]))
        }else{
          par[K+3L] <- sqrt(2*abs(cfm[3L]))
        }
      }else if(type==2){
        ind <- x>=0 & x<=1
        xx <- x[ind]
        yy <- y[ind]
        if(deg.is.odd){
          cfm <- coef(lm.fit(cbind(1,xx,xx^3),yy))
        }else{
          cfm <- coef(lm.fit(cbind(1,xx,xx^2),yy))
        }
        par <- rep(0, 2*K+2L+deg.is.odd)
        par[1L] <- cfm[1L]
        par[2L] <- ifelse(cor(xx,yy)>0, 1, -1)
        par[3L] <- sqrt(abs(cfm[2L]))
        if(deg.is.odd){
          par[4L] <- -par[3L] - sqrt(par[3L]^2 + 3*abs(cfm[3L]))
          par[K+4L] <- sqrt(abs(par[4L]^2 - 3*abs(cfm[3L])))
        }else{
          par[K+3L] <- sqrt(abs(cfm[2L]) + 2*abs(cfm[3L]))
        }
      }else{
        stop("How did we get here?")
      }
      if(!missing(monotone)){
        monotone <- match.arg(monotone)
        if(monotone=="increasing")
          par[2L] <- 1
        else
          par[2L] <- -1
      }
    }else{
      par <- start
    }
    
    if(plot.it){
      xgr <- seq(from=min(x), to=max(x), length=401)
      plot(x, y, xlab="x", ylab="y")
      ygr <- evalPolSOS(par, xgr, type, deg.is.odd, K)
      lines(xgr, ygr, col="green")
    }

    res <- FullSOS(x=x, y=y, w=w, deg.is.odd=deg.is.odd,
                   type=type, K=K, par=par,
                   trace=trace, plot.it=plot.it, control=control)

    res <- c(res, list(type=type))
  }
    
  
  
  tmp <- diag(rev(res$beta))
  MN <- NROW(tmp)
  storage.mode(tmp) <- "double"
  if(type==0){
    tmp <- .Fortran(.MP_transf,
                    c = tmp,
                    as.integer(MN), as.integer(MN),
                    scl = as.double(2/sclx),
                    mid = as.double(-(minx + sclx/2)))$c
  }else{
    tmp <- .Fortran(.MP_transf,
                    c = tmp,
                    as.integer(MN), as.integer(MN),
                    scl = as.double(1/sclx),
                    mid = as.double(-minx))$c
  }
  tmp <- rev(tmp[,1])
  tmp[1] <- tmp[1]+1
  tmp <- tmp * scly / 2
  tmp[1] <- tmp[1] + miny
  names(tmp) <- names(res$beta)
  res$beta.raw <- tmp
  if(!is.null(w))
    res$weights <- w
  
  res$fitted.values <- evalPol(x, res$beta)
  attributes(res$fitted.values) <- attributes(y)
  res$residuals <- y - res$fitted.values
  
  c(res,
    list(K=K, minx=minx, sclx=sclx, miny=miny, scly=scly,
         algorithm="Full", ptype="SOS"))
}

