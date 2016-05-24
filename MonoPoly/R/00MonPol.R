###  Copyright (C) 2011-2015 Berwin A. Turlach <Berwin.Turlach@gmail.com>
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
monpol.control <- function(maxiter=1000, tol=1e-05,
                           tol1=1e-10, tol2=1e-07, tolqr=1e-07){
  list(maxiter=maxiter, tol=tol, tol1=tol1, tol2=tol2, tolqr=tolqr)
}

monpol.fit <- function(x, y, w=NULL, K=1, start,
                       trace=FALSE, plot.it=FALSE,
                       control=monpol.control(),
                       algorithm=c("Full", "Hawkins", "BCD", "CD1", "CD2"),
                       ptype=c("Elphinstone", "EHH", "Penttila"),
                       ctype=c("cge0", "c2")){
  
  algorithm <- match.arg(algorithm)
  ptype <- match.arg(ptype)
  ctype <- match.arg(ctype)
  
  minx <- min(x)
  sclx <- max(x)-minx
  x <- 2*(x-minx)/sclx - 1

  miny <- min(y)
  scly <- max(y)-miny
  y <- 2*(y-miny)/scly - 1
  
  if(algorithm=="Hawkins"){
    res <- hawkins(x, y, w, K, trace, plot.it, control)
  }else{
    if(K==0){
      if(is.null(w)){ 
        tmp <- lm.fit(cbind(1, x), y)
      }else{
        tmp <- lm.wfit(cbind(1,x), y, w)
      }
      beta <- par <- coef(tmp)
      names(par) <- c("d", "a")
      names(beta) <- c("beta0", "beta1")
      res <- list(par=par, grad=c(0,0), beta=beta,
                  niter=0, converged=TRUE,
                  ptype=ptype, ctype=ctype)
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
        if(ptype=="Penttila"){
          cjs <- seq(from=1, to=2, length=K)
        }else{
          cjs <- seq(from=0.1, to=1, length=K)
        }
        if(ctype=="c2")
          cjs <- sqrt(cjs)
        par <- c(0, 0, as.numeric(rbind(0, cjs)))
      }else{
        par <- start
      }
      tmp <- lm.fit(cbind(1, evalPolMonPartilde(par, x, ptype, ctype)), y)
      par[1:2] <- as.vector(coef(tmp))
      names(par) <- c("d", "a", paste(rep(c("b","c"),K), rep(1:K, each=2), sep=""))
      
      if(plot.it){
        xgr <- seq(from=-1, to=1, length=401)
        plot(x, y, xlab="x", ylab="y")
        ygr <- evalPolMonPar(par, xgr, ptype=ptype, ctype=ctype)
        lines(xgr, ygr, col="green")
      }
      
      if(ctype=="cge0")
        algorithm <- paste(algorithm, "Constr", sep="")
      res <- do.call(algorithm,
                     list(x=x, y=y, w=w, K=K, par=par,
                          trace=trace, plot.it=plot.it, control=control,
                          ptype=ptype, ctype=ctype))
      res <- c(res,
               list(ptype=ptype, ctype=ctype))
    }    
  }
  
  tmp <- diag(rev(res$beta))
  storage.mode(tmp) <- "double"
  tmp <- .Fortran(.MP_transf,
                  c = tmp,
                  as.integer(2L*K+2L), as.integer(2L*K+2L),
                  scl = as.double(2/sclx),
                  mid = as.double(-(minx + sclx/2)))$c
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
         algorithm=algorithm))
}

monpol <- function (formula, data, subset, weights, na.action,
                    degree = 3, K, start,
                    a = -Inf, b = Inf,
                    trace=FALSE, plot.it=FALSE,
                    control=monpol.control(),
                    algorithm=c("Full", "Hawkins", "BCD", "CD1", "CD2"),
                    ptype=c("SOS", "Elphinstone", "EHH", "Penttila"),
                    ctype=c("cge0", "c2"),
                    monotone,
                    model = FALSE, x = FALSE, y = FALSE) 
{
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights",
                 "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) 
      stop("'weights' must be a numeric vector")
    if (is.empty.model(mt)) 
      stop("You should specify a regressor variable.")

    NoIntercept <- attr(mt, "intercept") == 0
    attr(mt, "intercept") <- 0
    x <- model.matrix(mt, mf)
    if(NCOL(x) != 1)
      stop("Regressor variable should be univariate.")

    algorithm <- match.arg(algorithm)
    ptype <- match.arg(ptype)
    ctype <- match.arg(ctype)
    
    if(!missing(K)){
      if(!missing(degree))
        stop("Either 'degree' or 'K' should be specified, not both.")
      if( trunc(K) != K || K < 0 )
        stop("'K' should be a nonnegative integer.")
      degree <- 2*K+1
    }else{
      if( trunc(degree) != degree || degree <= 0 )
        stop("'degree' should be a positive integer.")
      K <- (degree-1)/2
    }

    ## if 'degree' is specified, K might not be integer at this point
    ## but we also want to define 'degree.is.odd' without having the
    ## code in both arms of the above if(){...}else{..}
    if( !(degree.is.odd <- trunc(K) == K) ){
      K <- trunc(K)+1
    }

    type <- 3-(is.infinite(a)+2*is.infinite(b))
    if(type==2)
      stop("'a' must be finite if 'b' is finite.")
    if(type==0 && !degree.is.odd)
      stop("'degree' must be odd if 'a=-Inf' and 'b=Inf'.")
    if(type==3)
      type <- type-1
    if(type!=0 && ptype!="SOS" && !NoIntercept)
      stop("Parameterisation 'SOS' has to be used, or a no-intercept model, if 'a' or 'b', or both, are finite.")
    
    if(NoIntercept){
      z <- monpol.noint(x=x, y=y, w=w, deg.is.odd=degree.is.odd,
                        a=a, b=b, K=K, start=start,
                        trace=trace, plot.it=plot.it, control=control)
    }else{
      if(ptype=="SOS"){
        z <- SOSpol.fit(x=x, y=y, w=w, deg.is.odd=degree.is.odd,
                        a=a, b=b, monotone=monotone,
                        K=K, start=start,
                        trace=trace, plot.it=plot.it, type=type,
                        control=control)
      }else{
        z <- monpol.fit(x=x, y=y, w=w, K=K, start=start, trace=trace,
                        plot.it=plot.it, control=control,
                        algorithm=algorithm, ptype=ptype, ctype=ctype)
      }
    }
    
    z$na.action <- attr(mf, "na.action")
    z$call <- cl
    z$terms <- mt
    if (model) 
        z$model <- mf
    if (ret.x) 
        z$x <- x
    if (ret.y) 
        z$y <- y
    structure(z, class="monpol")
}
