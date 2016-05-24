## Mon Apr 15 04:04:20 2013
## Original file Copyright 2013 A.C. Guidoum
## This file is part of the R package kedd.
## Arsalane Chouaib GUIDOUM <acguidoum@usthb.dz> and <starsalane@gmail.com> 
## Department of Probabilities-Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algeris
## Algeria
##############################################################################


## Maximum-Likelihood (Kullback-Leibler information) Cross-Validation (MLCV) 

h.mlcv <- function(x, ...)  UseMethod("h.mlcv")

h.mlcv.default  <- function(x,lower=0.1,upper=5,tol=0.1 * lower,kernel=c("gaussian",
                            "epanechnikov","uniform","triangular","triweight",
                            "tricube","biweight","cosine"),...)
        {
     if (!is.numeric(x) || length(dim(x)) >=1 || length(x) < 2L) 
              stop("argument 'x' must be numeric and need at least 3 data points") 
     if (!is.numeric(lower) || lower < 0)
              stop("invalid 'lower'")
     if (!is.numeric(upper))
              stop("invalid 'upper'")
     if (!is.numeric(tol) || tol < 0)
              stop("invalid 'tol'")
     if (lower >= upper ) 
              stop("the boundaries must be positive and 'lower' must be smaller than 'upper'. Default boundaries were used")
     if (missing(kernel)) kernel <- "gaussian"
     name <- deparse(substitute(x))
     x <- x[!is.na(x)]
     x <- sort(x)
     n <- length(x)
     fmlcv <- function(h)
            {
      D <- kernel_fun_der(kernel, outer(x,x,"-")/h,deriv.order=0)
      diag(D) <- 0
      D <- (1/((n-1)*h))* colSums(D)
      mean(log(D))
      }
     obj <- optimize(fmlcv ,c(lower,upper),tol=tol,maximum = TRUE) 
     structure(list(x=x, data.name=name,n=n, kernel=kernel, h = obj$maximum, 
                    mlcv=obj$objective),class="h.mlcv")
}


###### 

print.h.mlcv <- function(x, digits=NULL, ...)
              {
    class(x) <- "h.mlcv"
    cat("\nCall:\t","\tMaximum-Likelihood Cross-Validation","\n",
        "\nData: ",x$data.name," (",x$n," obs.);","\tKernel: ",x$kernel, 
	    "\nMax CV = ",formatC(x$mlcv,digits=digits),";","\tBandwidth 'h' = ",formatC(x$h,digits=digits), "\n\n",sep="")
    invisible(x)
}

######

plot.mlcv <- function(f,seq.bws=NULL,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                      type="l",las=1,lwd=1,...)
                    {
    class(f) <- "h.mlcv"
    n <- f$n
    r <- 0
    kernel <- f$kernel
    x <- sort(f$x)
    if(is.null(xlab)) xlab <- "Bandwidths"
    if(is.null(ylab)) ylab <- bquote(MLCV~(h))                
    if(is.null(main)) main <- "Maximum-Likelihood Cross-Validation function for \nBandwidth Choice for Density Function"                
    if(is.null(sub)) sub <- paste("Kernel",kernel)					
    if(is.null(seq.bws)){
       hos <- ((243 *(2*r+1)*A3_kMr(kernel,r))/(35* A2_kM(kernel)^2))^(1/(2*r+5)) * sd(x,na.rm = TRUE) * n^(-1/(2*r+5))
       seq.bws <- seq(0.15*hos,2*hos,length=50)
                         }
    fmlcv <- function(h)
            {
        D <- kernel_fun_der(kernel, outer(x,x,"-")/h,deriv.order=0)
        diag(D) <- 0
        D <- (1/((n-1)*h))* colSums(D)
        mean(log(D))
        }
    D <- lapply(1:length(seq.bws), function(i) fmlcv(seq.bws[i]))
    Maxf <- c(do.call("rbind",D))
    plot.default(seq.bws,Maxf,type=type,las=las,lwd=lwd,xlab=xlab,ylab=ylab,
		         main=main,sub=sub,font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,...)
    return(list(kernel=kernel,seq.bws=seq.bws, mlcv=Maxf)) 
}


plot.h.mlcv <- function(x,seq.bws=NULL,...) plot.mlcv(x,seq.bws,...)

lines.mlcv <- function(f,seq.bws=NULL,...)
                    {
    class(f) <- "h.mlcv"
    r <- 0
    n <- f$n
    kernel <- f$kernel
    x <- sort(f$x)
    if(is.null(seq.bws)){
       hos <- ((243 *(2*r+1)*A3_kMr(kernel,r))/(35* A2_kM(kernel)^2))^(1/(2*r+5)) * sd(x,na.rm = TRUE) * n^(-1/(2*r+5))
       seq.bws <- seq(0.15*hos,2*hos,length=50)
                         }
    fmlcv <- function(h)
            {
        D <- kernel_fun_der(kernel, outer(x,x,"-")/h,deriv.order=0)
        diag(D) <- 0
        D <- (1/((n-1)*h))* colSums(D)
        mean(log(D))
        }
    D <- lapply(1:length(seq.bws), function(i) fmlcv(seq.bws[i]))
    Minf <- c(do.call("rbind",D))
    lines.default(seq.bws,Minf,...)
    invisible(NULL)
}

lines.h.mlcv <- function(x,seq.bws=NULL,...) lines.mlcv(x,seq.bws,...) 

