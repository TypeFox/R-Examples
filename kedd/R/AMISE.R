## Mon May 20 13:31:03 2013
## Original file Copyright 2013 A.C. Guidoum
## This file is part of the R package kedd.
## Arsalane Chouaib GUIDOUM <acguidoum@usthb.dz> and <starsalane@gmail.com> 
## Department of Probabilities-Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algeris
## Algeria
##############################################################################


## Aymptotic Mean Integrated Squared Error (AMISE)



h.amise <- function(x, ...)  UseMethod("h.amise")

h.amise.default <- function(x,deriv.order=0,lower=0.1*hos,upper=2*hos,tol=0.1 * lower,
                       kernel=c("gaussian","epanechnikov","triweight",
                                "tricube", "biweight","cosine"),...)
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
     if (kernel=="epanechnikov" && r+2 >= 3)    return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = NA ,amise=NA),class="h.amise"))
     else if (kernel=="triweight" && r+2 >= 7)  return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = NA ,amise=NA),class="h.amise"))
     else if (kernel=="biweight" && r+2 >= 5)   return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = NA ,amise=NA),class="h.amise"))
     else if (kernel=="tricube" && r+2 >= 10)   return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = NA ,amise=NA),class="h.amise"))
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
     R_Kr <- A3_kMr(kernel,r)
     famise <- function(h)
              {
       D <- ((-1)^(r+2)/h^(2*r+5)) * mean(kernel_fun_conv(kernel,outer(x,x,"-")/h,deriv.order=r+2))
       0.25*(2*r+5)* R_Kr^(4/(2*r+5)) * ((A2_kM(kernel)^2 * D) / (2*r+1))^((2*r+1)/(2*r+5)) * n^(-4/(2*r+5))    
              }
     obj <- optimize(famise , c(lower, upper),tol=tol)
     structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r, h = obj$minimum , 
                   amise=obj$objective),class="h.amise")
}

###### 

print.h.amise <- function(x, digits=NULL, ...)
              {
    class(x) <- "h.amise"
    cat("\nCall:\t","\tAymptotic Mean Integrated Squared Error","\n",
	   "\nDerivative order = ",x$deriv.order,
        "\nData: ",x$data.name," (",x$n," obs.);","\tKernel: ",x$kernel, 
	    "\nAMISE = ",format(x$amise,digits=digits),";","\tBandwidth 'h' = ",format(x$h,digits=digits), "\n\n",sep="")
    invisible(x)
}

######

plot.amise <- function(f,seq.bws=NULL,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                      type="l",las=1,lwd=1,...)
                    {
    class(f) <- "h.amise"
    r <- f$deriv.order
    n <- f$n
    kernel <- f$kernel
    x <- sort(f$x)
    if (kernel=="epanechnikov" && r+2 >= 3)    stop(" 'epanechnikov kernel derivative = 0' for 'order + 2 >= 3' ")
    else if (kernel=="triweight" && r+2 >= 7)  stop(" 'triweight kernel derivative = 0' for 'order + 2 >= 7' ")
    else if (kernel=="biweight" && r+2 >= 5)   stop(" 'biweight kernel derivative = 0' for 'order + 2 >= 5' ")
    else if (kernel=="tricube" && r+2 >= 10)   stop(" 'tricube kernel derivative = 0' for 'order + 2 >= 10' ")
    if(is.null(xlab)) xlab <- "Bandwidths"
    if(is.null(ylab)) ylab <-  bquote(AMISE~(h[(.(r))]))                
    if(is.null(main)){ 
	     if(r !=0) {main <- "Aymptotic Mean Integrated Squared Error for \nBandwidth Choice for Density Derivative"}else{
	                main <- "Aymptotic Mean Integrated Squared Error for \nBandwidth Choice for Density Function"}
	                }
    if(is.null(sub)) sub <- paste("Kernel",kernel,";","Derivative order = ",r)					
    if(is.null(seq.bws)){
       hos <- ((243 *(2*r+1)*A3_kMr(kernel,r))/(35* A2_kM(kernel)^2))^(1/(2*r+5)) * sd(x,na.rm = TRUE) * n^(-1/(2*r+5))
       seq.bws <- seq(0.15*hos,2*hos,length=50)
                         }
    R_Kr <- A3_kMr(kernel,r)
     famise <- function(h)
              {
       D <- ((-1)^(r+2)/h^(2*r+5)) * mean(kernel_fun_conv(kernel,outer(x,x,"-")/h,deriv.order=r+2))
       0.25*(2*r+5)* R_Kr^(4/(2*r+5)) * ((A2_kM(kernel)^2 * D) / (2*r+1))^((2*r+1)/(2*r+5)) * n^(-4/(2*r+5))    
              }
    D <- lapply(1:length(seq.bws), function(i) famise(seq.bws[i]))
    Minf <- c(do.call("rbind",D))
    plot.default(seq.bws,Minf,type=type,las=las,lwd=lwd,xlab=xlab,ylab=ylab,
		         main=main,sub=sub,font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,...)
    return(list(kernel=kernel,deriv.order=r,seq.bws=seq.bws, amise=Minf)) 
}


plot.h.amise <- function(x,seq.bws=NULL,...) plot.amise(x,seq.bws,...)

#####

lines.amise <- function(f,seq.bws=NULL,...)
                    {
    class(f) <- "h.amise"
    r <- f$deriv.order
    n <- f$n
    kernel <- f$kernel
	x <- sort(f$x)
    if (kernel=="epanechnikov" && r+2 >= 3)    stop(" 'epanechnikov kernel derivative = 0' for 'order + 2 >= 3' ")
    else if (kernel=="triweight" && r+2 >= 7)  stop(" 'triweight kernel derivative = 0' for 'order + 2 >= 7' ")
    else if (kernel=="biweight" && r+2 >= 5)   stop(" 'biweight kernel derivative = 0' for 'order + 2 >= 5' ")
    else if (kernel=="tricube" && r+2 >= 10)   stop(" 'tricube kernel derivative = 0' for 'order + 2 >= 10' ")
    if(is.null(seq.bws)){
       hos <- ((243 *(2*r+1)*A3_kMr(kernel,r))/(35* A2_kM(kernel)^2))^(1/(2*r+5)) * sd(x,na.rm = TRUE) * n^(-1/(2*r+5))
       seq.bws <- seq(0.15*hos,2*hos,length=50)
                         }
    R_Kr <- A3_kMr(kernel,r)
     famise <- function(h)
              {
       D <- ((-1)^(r+2)/h^(2*r+5)) * mean(kernel_fun_conv(kernel,outer(x,x,"-")/h,deriv.order=r+2))
       0.25*(2*r+5)* R_Kr^(4/(2*r+5)) * ((A2_kM(kernel)^2 * D) / (2*r+1))^((2*r+1)/(2*r+5)) * n^(-4/(2*r+5))    
              }
    D <- lapply(1:length(seq.bws), function(i) famise(seq.bws[i]))
    Minf <- c(do.call("rbind",D))
    lines.default(seq.bws,Minf,...)
    invisible(NULL)
}

lines.h.amise <- function(x,seq.bws=NULL,...) lines.amise(x,seq.bws,...) 
