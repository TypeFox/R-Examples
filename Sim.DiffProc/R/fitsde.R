## Wed Dec 30 22:05:18 2015
## Original file Copyright © 2016 A.C. Guidoum, K. Boukhetala
## This file is part of the R package Sim.DiffProc
## Department of Probabilities & Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algiers
## Algeria

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## A copy of the GNU General Public License is available at
## http://www.r-project.org/Licenses/
## Unlimited use and distribution (see LICENCE).
###################################################################################################

###
### fitsde

fitsde <- function(data, ...)
              { 
    call <- match.call()
    UseMethod("fitsde")
}

fitsde.default <- function(data,drift,diffusion,start = list(),
                  pmle=c("euler","kessler","ozaki","shoji"),
                  optim.method = "L-BFGS-B",lower = -Inf, upper = Inf,...)
               {
    if (is.ts(data) == FALSE) stop("data must be time-series objects")
    if (!missing(start) && (!is.list(start) || is.null(names(start)))) stop("'start' must be a named list")
    if (any(!is.expression(drift) || !is.expression(diffusion) )) stop(" coefficient of 'drift' and 'diffusion' must be expressions")   
    pmle <- match.arg(pmle)
    .Euler.lik <- function(theta)
              {
     n <- length(data)
     -sum(.dcEuler(x=data[2:n],t=time(data)[2:n],x0=data[1:(n-1)],t0=time(data)[1:(n-1)],
          theta,drift,diffusion,log=TRUE),na.rm=TRUE)
       }
    .Ozaki.lik <- function(theta)
              {
     n <- length(data)
     -sum(.dcOzaki(x=data[2:n],t=time(data)[2:n],x0=data[1:(n-1)],t0=time(data)[1:(n-1)],
          theta,drift,diffusion,log=TRUE),na.rm=TRUE)
       }
    .Shoji.lik <- function(theta)
              {
     n <- length(data)
     -sum(.dcShoji(x=data[2:n],t=time(data)[2:n],x0=data[1:(n-1)],t0=time(data)[1:(n-1)],
          theta,drift,diffusion,log=TRUE),na.rm=TRUE)
       }
    #.Elerian.lik <- function(theta)
    #          {
    # n <- length(data)
    # -sum(.dcElerian(x=data[2:n],t=time(data)[2:n],x0=data[1:(n-1)],t0=time(data)[1:(n-1)],
    #      theta,drift,diffusion,log=TRUE),na.rm=TRUE)
    #   }
    .Kessler.lik <- function(theta)
              {
     n <- length(data)
     -sum(.dcKessler(x=data[2:n],t=time(data)[2:n],x0=data[1:(n-1)],t0=time(data)[1:(n-1)],
          theta,drift,diffusion,log=TRUE),na.rm=TRUE)
       }
    if (pmle == 'euler')       {f <- .Euler.lik}
    else if (pmle == 'ozaki')  {f <- .Ozaki.lik}
    else if (pmle== 'shoji')   {f <- .Shoji.lik}
    else if (pmle== 'kessler') {f <- .Kessler.lik}
	#else if (pmle== 'elerian') {f <- .Elerian.lik}
    out <- optim(start, f, method = optim.method, hessian = TRUE,lower=lower,upper=upper,...)
    coef=out$par
    structure(list(call=call,data=data,drift=drift, diffusion=diffusion,method=pmle,value=out$value,
                   optim.method=optim.method,coef=coef,hessian=out$hessian),class="fitsde")
}
############

print.fitsde <- function(x, digits=NULL, ...)
           {
    class(x) <- "fitsde"
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coef)
    invisible(x)
}

summary.fitsde <- function(object, ...)
           {
    x <- object
    class(x) <- "fitsde"
    if (x$method=="euler")         {meth <- "Euler"}
    else if (x$method=="kessler")  {meth <- "Kessler"}
    else if (x$method=="ozaki")    {meth <- "Ozaki"}
    else if (x$method=="shoji")    {meth <- "Shoji"}
    #else if (x$method=="elerian")  {meth <- "Elerian"}	
    vcov <- if (length(x$coef)) 
            solve(x$hessian)
    m2logL <- 2*x$value
    cat("Pseudo maximum likelihood estimation","\n",sep="")
    cat("\nMethod: ",meth)
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    cmat <- cbind(Estimate = x$coef,`Std. Error` = sqrt(diag(vcov)))
    print(cmat)
    cat("\n-2 log L:", m2logL, "\n")
    invisible(x)
}

vcov.fitsde <- function(object, ...)
       {
    x <- object
    class(x) <- "fitsde"
    vcov <- if (length(x$coef)) 
            solve(x$hessian)
    return(vcov)
    invisible(x)
}

AIC.fitsde <- function(object, ...)
       {
    x <- object
    class(x) <- "fitsde"
    n <- length(x$coef)
    aic <- 2*x$value + 2*n
    return(aic)
    invisible(x)
}

BIC.fitsde <- function(object, ...)
       {
    x <- object
    class(x) <- "fitsde"
    n <- length(x$data)
    bic <- 2*x$value + 2*log(n)
    return(bic)
    invisible(x)
}

logLik.fitsde <- function(object, ...)
         {
    x <- object
    class(x) <- "fitsde"
    return(-x$value)  
    invisible(x)
}

coef.fitsde <- function(object,...)
            {
    x <- object
    class(x) <- "fitsde"
    return(x$coef)
    invisible(x)
}

confint.fitsde <- function(object,parm,level=0.95, ...)
       {
    x <- object
    class(x) <- "fitsde"
    n <- length(x$coef)
    pnames <- names(x$coef)
    if (missing(parm)) parm <- pnames
    vcov <- if (length(x$coef)) 
            solve(x$hessian)
    cf <- (1 - level)/2
    a <- c(cf, 1 - cf)
    conf <- data.frame(do.call("rbind",lapply(1:n, function(i) 
                      c(x$coef[[i]]+ sqrt(diag(vcov))[[i]] * qnorm(cf),
                        x$coef[[i]]- sqrt(diag(vcov))[[i]] * qnorm(cf)) ) ) )
    rownames(conf) <- names(x$coef)
    names(conf)     <- paste(round(100 * a, 1), "%")
    conf <- conf[parm,]
    return(conf)
}
