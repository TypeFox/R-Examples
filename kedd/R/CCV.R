## Mon Jun 10 02:18:40 2013
## Original file Copyright 2013 A.C. Guidoum
## This file is part of the R package kedd.
## Arsalane Chouaib GUIDOUM <acguidoum@usthb.dz> and <starsalane@gmail.com> 
## Department of Probabilities-Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algeris
## Algeria
##############################################################################


## Complete Cross-Validation (CCV)


h.ccv <- function(x, ...)  UseMethod("h.ccv")

h.ccv.default <- function(x,deriv.order=0,lower=0.1*hos,upper=hos,tol=0.1 * lower,
                          kernel=c("gaussian","triweight","tricube",
                                   "biweight","cosine"),...)
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
     if (kernel=="triweight" && 2*(r+2) >= 7)       return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = NA ,min.ccv=NA),class="h.ccv"))
     else if (kernel=="biweight" && 2*(r+2) >= 5)   return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = NA ,min.ccv=NA),class="h.ccv"))
     else if (kernel=="tricube" && 2*(r+2) >= 10)   return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = NA ,min.ccv=NA),class="h.ccv"))  
     hos <- ((243 *(2*r+1)*A3_kMr(kernel,r))/(35* A2_kM(kernel)^2))^(1/(2*r+5)) * sd(x,na.rm = TRUE) * n^(-1/(2*r+5))
     if (!is.numeric(upper)){ 
		stop("argument 'upper' must be numeric. Default 2*hos (Oversmoothing) boundary was used")
		upper= hos
	}
	if (!is.numeric(lower)){
		stop("argument 'lower' must be numeric. Default 0.1*hos boundary was used")
		lower=0.1*hos
	}
	if (lower < 0 | lower >= upper){
      	stop("the boundaries must be positive and 'lower' must be smaller than 'upper'. Default boundaries were used")
		upper=hos
		lower=0.1*hos
	}
	R_Kr1 <- A3_kMr(kernel,r)
     fccv <- function(h)
              {		
		L1 <- kernel_fun_conv(kernel,outer(x,x,"-")/h,deriv.order=r)
        diag(L1) <- 0
        L2 <- ((-1)^(r)/((n-1)*h^(2*r+1)))* colSums(L1)
        Q1 <- mean(L2)	  
        D2 <- kernel_fun_der(kernel,outer(x,x,"-")/h,deriv.order=2*r)
        diag(D2) <- 0
        D2 <- ((-1)^r /((n-1)*h^(2*r+1)))* colSums(D2)
        Q2 <- mean(D2)
        D3 <- kernel_fun_der(kernel,outer(x,x,"-")/h,deriv.order=2*(r+1))
        diag(D3) <- 0
        D3 <- ((-1)^(r+1) /((n-1)*h^(2*r+3)))* colSums(D3)
        Q3 <- mean(D3)
        D4 <- kernel_fun_der(kernel,outer(x,x,"-")/h,deriv.order=2*(r+2))
        diag(D4) <- 0
        D4 <- ((-1)^(r+2) /((n-1)*h^(2*r+5)))* colSums(D4)
        Q4 <- mean(D4)
        (1/(n*h^(2*r+1)))* R_Kr1 + Q1 - Q2 + 0.5 * h^2 *A2_kM(kernel) * Q3 + (h^4 / 24) *(6*(A2_kM(kernel))^2 - A5_kM(kernel)) * Q4
              }
    obj <- optimize(fccv ,c(lower, upper),tol=tol)
    structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r,h = obj$minimum , 
                   min.ccv=obj$objective),class="h.ccv")  
}

###### 

print.h.ccv <- function(x, digits=NULL, ...)
              {
    class(x) <- "h.ccv"
    cat("\nCall:\t","\tComplete Cross-Validation","\n",
	   "\nDerivative order = ",x$deriv.order,
        "\nData: ",x$data.name," (",x$n," obs.);","\tKernel: ",x$kernel, 
	    "\nMin CCV = ",format(x$min.ccv,digits=digits),";","\tBandwidth 'h' = ",format(x$h,digits=digits), "\n\n",sep="")
    invisible(x)
}

######

plot.ccv <- function(f,seq.bws=NULL,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                      type="l",las=1,lwd=1,...)
                    {
    class(f) <- "h.ccv"
    r <- f$deriv.order
    n <- f$n
    kernel <- f$kernel
    x <- sort(f$x)
    if (kernel=="triweight" && 2*(r+2) >= 7)       stop(" 'triweight kernel derivative = 0' for '2 * (order + 2) >= 7' ")
    else if (kernel=="biweight" && 2*(r+2) >= 5)   stop(" 'biweight kernel derivative = 0' for '2 * (order + 2) >= 5' ")
    else if (kernel=="tricube" && 2*(r+2) >= 10)   stop(" 'tricube kernel derivative = 0' for '2 * (order + 2) >= 10' ")
    if(is.null(xlab)) xlab <- "Bandwidths"
    if(is.null(ylab)) ylab <- bquote(CCV~(h[(.(r))]))                
    if(is.null(main)){ 
	     if(r !=0) {main <- "Complete Cross-Validation function for \nBandwidth Choice for Density Derivative"}else{
	                main <- "Complete Cross-Validation function for \nBandwidth Choice for Density Function"}
	                }
    if(is.null(sub)) sub <- paste("Kernel",kernel,";","Derivative order = ",r)					
    if(is.null(seq.bws)){
       hos <- ((243 *(2*r+1)*A3_kMr(kernel,r))/(35* A2_kM(kernel)^2))^(1/(2*r+5)) * sd(x,na.rm = TRUE) * n^(-1/(2*r+5))
       seq.bws <- seq(0.15*hos,2*hos,length=50)
                         }
	R_Kr1 <- A3_kMr(kernel,r)
     fccv <- function(h)
              {		
		L1 <- kernel_fun_conv(kernel,outer(x,x,"-")/h,deriv.order=r)
        diag(L1) <- 0
        L2 <- ((-1)^(r)/((n-1)*h^(2*r+1)))* colSums(L1)
        Q1 <- mean(L2)	  
        D2 <- kernel_fun_der(kernel,outer(x,x,"-")/h,deriv.order=2*r)
        diag(D2) <- 0
        D2 <- ((-1)^r /((n-1)*h^(2*r+1)))* colSums(D2)
        Q2 <- mean(D2)
        D3 <- kernel_fun_der(kernel,outer(x,x,"-")/h,deriv.order=2*(r+1))
        diag(D3) <- 0
        D3 <- ((-1)^(r+1) /((n-1)*h^(2*r+3)))* colSums(D3)
        Q3 <- mean(D3)
        D4 <- kernel_fun_der(kernel,outer(x,x,"-")/h,deriv.order=2*(r+2))
        diag(D4) <- 0
        D4 <- ((-1)^(r+2) /((n-1)*h^(2*r+5)))* colSums(D4)
        Q4 <- mean(D4)
        (1/(n*h^(2*r+1)))* R_Kr1 + Q1 - Q2 + 0.5 * h^2 *A2_kM(kernel) * Q3 + (h^4 / 24) *(6*(A2_kM(kernel))^2 - A5_kM(kernel)) * Q4
              }
    D <- lapply(1:length(seq.bws), function(i) fccv(seq.bws[i]))
    Minf <- c(do.call("rbind",D))
    plot.default(seq.bws,Minf,type=type,las=las,lwd=lwd,xlab=xlab,ylab=ylab,
		         main=main,sub=sub,font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,...)
    return(list(kernel=kernel,deriv.order=r,seq.bws=seq.bws, ccv=Minf))
}

plot.h.ccv <- function(x,seq.bws=NULL,...) plot.ccv(x,seq.bws,...)

#####

lines.ccv <- function(f,seq.bws=NULL,...)
                    {
    class(f) <- "h.ccv"
    r <- f$deriv.order
    n <- f$n
    kernel <- f$kernel
    x <- sort(f$x)
    if (kernel=="triweight" && 2*(r+2) >= 7)       stop(" 'triweight kernel derivative = 0' for '2 * (order + 2) >= 7' ")
    else if (kernel=="biweight" && 2*(r+2) >= 5)   stop(" 'biweight kernel derivative = 0' for '2 * (order + 2) >= 5' ")
    else if (kernel=="tricube" && 2*(r+2) >= 10)   stop(" 'tricube kernel derivative = 0' for '2 * (order + 2) >= 10' ")
    if(is.null(seq.bws)){
       hos <- ((243 *(2*r+1)*A3_kMr(kernel,r))/(35* A2_kM(kernel)^2))^(1/(2*r+5)) * sd(x,na.rm = TRUE) * n^(-1/(2*r+5))
       seq.bws <- seq(0.15*hos,2*hos,length=50)
                         }
	R_Kr1 <- A3_kMr(kernel,r)
     fccv <- function(h)
              {		
		L1 <- kernel_fun_conv(kernel,outer(x,x,"-")/h,deriv.order=r)
        diag(L1) <- 0
        L2 <- ((-1)^(r)/((n-1)*h^(2*r+1)))* colSums(L1)
        Q1 <- mean(L2)  
        D2 <- kernel_fun_der(kernel,outer(x,x,"-")/h,deriv.order=2*r)
        diag(D2) <- 0
        D2 <- ((-1)^r /((n-1)*h^(2*r+1)))* colSums(D2)
        Q2 <- mean(D2)
        D3 <- kernel_fun_der(kernel,outer(x,x,"-")/h,deriv.order=2*(r+1))
        diag(D3) <- 0
        D3 <- ((-1)^(r+1) /((n-1)*h^(2*r+3)))* colSums(D3)
        Q3 <- mean(D3)
        D4 <- kernel_fun_der(kernel,outer(x,x,"-")/h,deriv.order=2*(r+2))
        diag(D4) <- 0
        D4 <- ((-1)^(r+2) /((n-1)*h^(2*r+5)))* colSums(D4)
        Q4 <- mean(D4)
        (1/(n*h^(2*r+1)))* R_Kr1 + Q1 - Q2 + 0.5 * h^2 *A2_kM(kernel) * Q3 + (h^4 / 24) *(6*(A2_kM(kernel))^2 - A5_kM(kernel)) * Q4
              }
    D <- lapply(1:length(seq.bws), function(i) fccv(seq.bws[i]))
    Minf <- c(do.call("rbind",D))
    lines.default(seq.bws,Minf,...)
    invisible(NULL)
}

lines.h.ccv <- function(x,seq.bws=NULL,...) lines.ccv(x,seq.bws,...) 
