## Fri Jun 21 01:27:08 2013
## Original file Copyright 2013 A.C. Guidoum
## This file is part of the R package kedd.
## Arsalane Chouaib GUIDOUM <acguidoum@usthb.dz> and <starsalane@gmail.com> 
## Department of Probabilities-Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algeris
## Algeria
##############################################################################



## Trimmed Cross-Validation (TCV)


h.tcv <- function(x, ...)  UseMethod("h.tcv")

h.tcv.default <- function(x,deriv.order=0,lower=0.1*hos,upper=2*hos,tol=0.1 * lower,
                          kernel=c("gaussian","epanechnikov", "uniform",
                                   "triangular","triweight","tricube", "biweight",
                                   "cosine"),...)
   {
     if (!is.numeric(x) || length(dim(x)) >=1 || length(x) < 3L) 
           stop("argument 'x' must be numeric and need at least 3 data points") 
     if (any(deriv.order < 0 || deriv.order != round(deriv.order))) 
           stop("argument 'deriv.order' is non-negative integers")
     if (missing(kernel)) kernel <- "gaussian"
     r <- deriv.order
     name <- deparse(substitute(x))
     x <- x[!is.na(x)]
     x <- sort(x)
     n <- length(x)
     if (kernel=="epanechnikov" && 2*r >= 3)    return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = NA ,min.tcv=NA),class="h.tcv"))
     else if (kernel=="uniform" && 2*r >= 1)    return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = NA ,min.tcv=NA),class="h.tcv"))
     else if (kernel=="triweight" && 2*r >= 7)  return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = NA ,min.tcv=NA),class="h.tcv"))
     else if (kernel=="biweight" && 2*r >= 5)   return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = NA ,min.tcv=NA),class="h.tcv"))
     else if (kernel=="triangular" && 2*r >= 2) return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = NA ,min.tcv=NA),class="h.tcv"))
     else if (kernel=="tricube" && 2*r >= 10)   return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = NA ,min.tcv=NA),class="h.tcv"))
     hos <- ((243 *(2*r+1)*A3_kMr(kernel,r))/(35* A2_kM(kernel)^2))^(1/(2*r+5)) * sd(x,na.rm = TRUE) * n^(-1/(2*r+5))
     if (!is.numeric(upper)){ 
		stop("argument 'upper' must be numeric. Default 2*hos (Oversmoothing) boundary was used")
		upper= 2*hos
	}
	if (!is.numeric(lower)){
		stop("argument 'lower' must be numeric. Default 0.1*hos boundary was used")
		lower=0.1*hos
	}
	if (lower < 0 | lower >= upper){
      	stop("the boundaries must be positive and 'lower' must be smaller than 'upper'. Default boundaries were used")
		upper=2*hos
		lower=0.1*hos
	}
    R_Kr1 <- A3_kMr(kernel,r)
     ftcv <- function(h)
              {
	  L1 <- kernel_fun_conv(kernel,outer(x,x,"-")/h,deriv.order=r)
      diag(L1) <- 0
      L2 <- ((-1)^(r)/((n-1)*h^(2*r+1)))* colSums(L1)
      Q1 <- mean(L2)	  
      D1 <- kernel_fun_der(kernel, outer(x,x,"-")/h,deriv.order=2*r)*(abs(outer(x,x,"-")/h)> 1/(n*h^(2*r+1))) 
      diag(D1) <- 0
      D2 <- ((-1)^r / ((n-1)*h^(2*r+1)))* colSums(D1)
      Q2 <- mean(D2)
      (1/(n*h^(2*r+1)))* R_Kr1 + Q1 - 2 * Q2  
              }
     obj <- optimize(ftcv , c(lower, upper),tol=tol)
     structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r,h = obj$minimum , 
                    min.tcv=obj$objective),class="h.tcv")  
}

###### 

print.h.tcv <- function(x, digits=NULL, ...)
              {
    class(x) <- "h.tcv"
    cat("\nCall:\t","\tTrimmed Cross-Validation","\n",
	   "\nDerivative order = ",x$deriv.order,
        "\nData: ",x$data.name," (",x$n," obs.);","\tKernel: ",x$kernel, 
	    "\nMin TCV = ",format(x$min.tcv,digits=digits),";","\tBandwidth 'h' = ",format(x$h,digits=digits), "\n\n",sep="")
    invisible(x)
}

######

plot.tcv <- function(f,seq.bws=NULL,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                      type="l",las=1,lwd=1,...)
                    {
    class(f) <- "h.tcv"
    r <- f$deriv.order
    n <- f$n
    kernel <- f$kernel
    x <- sort(f$x)
    if(is.null(xlab)) xlab <- "Bandwidths"
    if(is.null(ylab)) ylab <- bquote(TCV~(h[(.(r))]))                
    if(is.null(main)){ 
	     if(r !=0) {main <- "Trimmed Cross-Validation function for \nBandwidth Choice for Density Derivative"}else{
	                main <- "Trimmed Cross-Validation function for \nBandwidth Choice for Density Function"}
	                }
    if(is.null(sub)) sub <- paste("Kernel",kernel,";","Derivative order = ",r)					
    if(is.null(seq.bws)){
       hos <- ((243 *(2*r+1)*A3_kMr(kernel,r))/(35* A2_kM(kernel)^2))^(1/(2*r+5)) * sd(x,na.rm = TRUE) * n^(-1/(2*r+5))
       seq.bws <- seq(0.15*hos,2*hos,length=50)
                         }
    if (kernel=="epanechnikov" && 2*r >= 3)    stop(" 'epanechnikov kernel derivative = 0' for '2*order >= 3' ")
    else if (kernel=="uniform" && 2*r >= 1)    stop(" 'uniform kernel derivative = 0' for '2*order >= 1' ")
    else if (kernel=="triweight" && 2*r >= 7)  stop(" 'triweight kernel derivative = 0' for '2*order >= 7' ")
    else if (kernel=="biweight" && 2*r >= 5)   stop(" 'biweight kernel derivative = 0' for '2*order >= 5' ")
    else if (kernel=="triangular" && 2*r >= 2) stop(" 'triangular kernel derivative = 0' for '2*order >= 2' ")
    else if (kernel=="tricube" && 2*r >= 10)   stop(" 'tricube kernel derivative = 0' for '2*order >= 10' ")
    R_Kr1 <- A3_kMr(kernel,r)
     ftcv <- function(h)
              {
	  L1 <- kernel_fun_conv(kernel,outer(x,x,"-")/h,deriv.order=r)
      diag(L1) <- 0
      L2 <- ((-1)^(r)/((n-1)*h^(2*r+1)))* colSums(L1)
      Q1 <- mean(L2)	  
      D1 <- kernel_fun_der(kernel, outer(x,x,"-")/h,deriv.order=2*r)*(abs(outer(x,x,"-")/h)> 1/(n*h^(2*r+1))) 
      diag(D1) <- 0
      D2 <- ((-1)^r / ((n-1)*h^(2*r+1)))* colSums(D1)
      Q2 <- mean(D2)
      (1/(n*h^(2*r+1)))* R_Kr1 + Q1 - 2 * Q2  
              }
    D <- lapply(1:length(seq.bws), function(i) ftcv(seq.bws[i]))
    Minf <- c(do.call("rbind",D))
    plot.default(seq.bws,Minf,type=type,las=las,lwd=lwd,xlab=xlab,ylab=ylab,
		         main=main,sub=sub,font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,...)
    return(list(kernel=kernel,deriv.order=r,seq.bws=seq.bws, tcv=Minf))
}

plot.h.tcv <- function(x,seq.bws=NULL,...) plot.tcv(x,seq.bws,...)


#####

lines.tcv <- function(f,seq.bws=NULL,...)
                    {
    class(f) <- "h.tcv"
    r <- f$deriv.order
    n <- f$n
    kernel <- f$kernel
    x <- sort(f$x)
    if (kernel=="epanechnikov" && 2*r >= 3)    stop(" 'epanechnikov kernel derivative = 0' for '2*order >= 3' ")
    else if (kernel=="uniform" && 2*r >= 1)    stop(" 'uniform kernel derivative = 0' for '2*order >= 1' ")
    else if (kernel=="triweight" && 2*r >= 7)  stop(" 'triweight kernel derivative = 0' for '2*order >= 7' ")
    else if (kernel=="biweight" && 2*r >= 5)   stop(" 'biweight kernel derivative = 0' for '2*order >= 5' ")
    else if (kernel=="triangular" && 2*r >= 2) stop(" 'triangular kernel derivative = 0' for '2*order >= 2' ")
    else if (kernel=="tricube" && 2*r >= 10)   stop(" 'tricube kernel derivative = 0' for '2*order >= 10' ")
    if(is.null(seq.bws)){
       hos <- ((243 *(2*r+1)*A3_kMr(kernel,r))/(35* A2_kM(kernel)^2))^(1/(2*r+5)) * sd(x,na.rm = TRUE) * n^(-1/(2*r+5))
       seq.bws <- seq(0.15*hos,2*hos,length=50)
                         }
    R_Kr1 <- A3_kMr(kernel,r)
     ftcv <- function(h)
              {
	  L1 <- kernel_fun_conv(kernel,outer(x,x,"-")/h,deriv.order=r)
      diag(L1) <- 0
      L2 <- ((-1)^(r)/((n-1)*h^(2*r+1)))* colSums(L1)
      Q1 <- mean(L2)	  
      D1 <- kernel_fun_der(kernel, outer(x,x,"-")/h,deriv.order=2*r)*(abs(outer(x,x,"-")/h)> 1/(n*h^(2*r+1))) 
      diag(D1) <- 0
      D2 <- ((-1)^r / ((n-1)*h^(2*r+1)))* colSums(D1)
      Q2 <- mean(D2)
      (1/(n*h^(2*r+1)))* R_Kr1 + Q1 - 2 * Q2  
              }
    D <- lapply(1:length(seq.bws), function(i) ftcv(seq.bws[i]))
    Minf <- c(do.call("rbind",D))
    lines.default(seq.bws,Minf,...)
    invisible(NULL)
}

lines.h.tcv <- function(x,seq.bws=NULL,...) lines.tcv(x,seq.bws,...) 
