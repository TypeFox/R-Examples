## Tue May 21 02:04:46 2013
## Original file Copyright 2013 A.C. Guidoum
## This file is part of the R package kedd.
## Arsalane Chouaib GUIDOUM <acguidoum@usthb.dz> and <starsalane@gmail.com> 
## Department of Probabilities-Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algeris
## Algeria
##############################################################################


## Biased Cross-Validation (BCV)


h.bcv <- function(x, ...)  UseMethod("h.bcv")

h.bcv.default <- function(x,whichbcv = 1,deriv.order=0,lower=0.1*hos,upper=2*hos,tol=0.1 * lower,
                          kernel=c("gaussian","epanechnikov","triweight","tricube",
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
     if(whichbcv == 1){
     if (kernel=="epanechnikov" && r+2 >= 3)    return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r,whichbcv=whichbcv, h = NA ,min.bcv1=NA),class="h.bcv"))
     else if (kernel=="triweight" && r+2 >= 7)  return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r,whichbcv=whichbcv, h = NA ,min.bcv1=NA),class="h.bcv"))
     else if (kernel=="biweight" && r+2 >= 5)   return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r,whichbcv=whichbcv, h = NA ,min.bcv1=NA),class="h.bcv"))
     else if (kernel=="tricube" && r+2 >= 10)   return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r,whichbcv=whichbcv, h = NA ,min.bcv1=NA),class="h.bcv"))
     } else if (whichbcv == 2) {
     if (kernel=="epanechnikov" && 2*(r+2) >= 3)    return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r,whichbcv=whichbcv, h = NA ,min.bcv2=NA),class="h.bcv"))
     else if (kernel=="triweight" && 2*(r+2) >= 7)  return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r,whichbcv=whichbcv, h = NA ,min.bcv2=NA),class="h.bcv"))
     else if (kernel=="biweight" && 2*(r+2) >= 5)   return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r,whichbcv=whichbcv, h = NA ,min.bcv2=NA),class="h.bcv"))
     else if (kernel=="tricube" && 2*(r+2) >= 10)   return(structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r,whichbcv=whichbcv, h = NA ,min.bcv2=NA),class="h.bcv"))
     }
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
      R_Kr2 <- A3_kMr(kernel,r+2)
      fbcv1 <- function(h)
              {
         D1 <- kernel_fun_conv(kernel,outer(x,x,"-")/h,deriv.order=r+2)
         diag(D1) <- 0
         D1 <- ((-1)^(r+2)/((n-1)*h^(2*r+5)))* colSums(D1)
         D <- mean(D1)
         ##(1/(n*h^(2*r+1)))* R_Kr1 + (0.25*h^4)*(A2_kM(kernel))^2 * (D - (1/((n-1)*h^(2*r+5))) * R_Kr2 )
		 (1/(n*h^(2*r+1)))* R_Kr1 + (0.25*h^4)*(A2_kM(kernel))^2 * D
              }
      fbcv2 <- function(h)
              {
         D1 <- kernel_fun_der(kernel,outer(x,x,"-")/h,deriv.order=2*(r+2))
         diag(D1) <- 0
         D1 <- ((-1)^(r+2)/((n-1)*h^(2*r+5)))* colSums(D1)
         D <- mean(D1)
         (1/(n*h^(2*r+1)))* R_Kr1 + (0.25*h^4)*(A2_kM(kernel))^2 *  D
              }
      if(whichbcv == 1) {obj <- optimize(fbcv1 ,c(lower, upper),tol=tol)
      structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r,whichbcv=whichbcv, h = obj$minimum , 
                   min.bcv1=obj$objective),class="h.bcv")}else{
      obj <- optimize(fbcv2 ,c(lower, upper),tol=tol)  
      structure(list(x=x, data.name=name,n=n, kernel=kernel, deriv.order=r,whichbcv=whichbcv, h = obj$minimum , 
                   min.bcv2=obj$objective),class="h.bcv")}
}

###### 

print.h.bcv <- function(x, digits=NULL, ...)
              {
    class(x) <- "h.bcv"
    if (x$whichbcv == 1){
     cat("\nCall:\t","\tBiased Cross-Validation 1","\n",
	   "\nDerivative order = ",x$deriv.order,
        "\nData: ",x$data.name," (",x$n," obs.);","\tKernel: ",x$kernel, 
	    "\nMin BCV = ",format(x$min.bcv1,digits=digits),";","\tBandwidth 'h' = ",format(x$h,digits=digits), "\n\n",sep="")}else{
     cat("\nCall:\t","\tBiased Cross-Validation 2","\n",
	   "\nDerivative order = ",x$deriv.order,
        "\nData: ",x$data.name," (",x$n," obs.);","\tKernel: ",x$kernel, 
	    "\nMin BCV = ",format(x$min.bcv2,digits=digits),";","\tBandwidth 'h' = ",format(x$h,digits=digits), "\n\n",sep="")
   }
    invisible(x)
}

######

plot.bcv <- function(f,seq.bws=NULL,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                      type="l",las=1,lwd=1,...)
                    {
    class(f) <- "h.bcv"
    r <- f$deriv.order
    n <- f$n
    kernel <- f$kernel
    x <- sort(f$x)
    if(f$whichbcv == 1){
    if (kernel=="epanechnikov" && r+2 >= 3)    stop(" 'epanechnikov kernel derivative = 0' for 'order + 2 >= 3' ")
    else if (kernel=="triweight" && r+2 >= 7)  stop(" 'triweight kernel derivative = 0' for 'order + 2 >= 7' ")
    else if (kernel=="biweight" && r+2 >= 5)   stop(" 'biweight kernel derivative = 0' for 'order + 2 >= 5' ")
    else if (kernel=="tricube" && r+2 >= 10)   stop(" 'tricube kernel derivative = 0' for 'order + 2 >= 10' ")
    } else if (f$whichbcv == 2){
    if (kernel=="epanechnikov" && 2*(r+2) >= 3)    stop(" 'epanechnikov kernel derivative = 0' for '2 * (order + 2) >= 3' ")
    else if (kernel=="triweight" && 2*(r+2) >= 7)  stop(" 'triweight kernel derivative = 0' for '2 * (order + 2) >= 7' ")
    else if (kernel=="biweight" && 2*(r+2) >= 5)   stop(" 'biweight kernel derivative = 0' for '2 * (order + 2) >= 5' ")
    else if (kernel=="tricube" && 2*(r+2) >= 10)   stop(" 'tricube kernel derivative = 0' for '2 * (order + 2) >= 10' ")
    }
    if(is.null(xlab)) xlab <- "Bandwidths"
    if(is.null(ylab)) ylab <-  bquote(BCV~(h[(.(r))]))               
    if(is.null(main)){
        if(f$whichbcv == 1){ 
	     if(r !=0) {main <- "Biased Cross-Validation (1) function for \nBandwidth Choice for Density Derivative"}else{
	                main <- "Biased Cross-Validation (1) function for \nBandwidth Choice for Density Function"}}else{
          if(r !=0) {main <- "Biased Cross-Validation (2) function for \nBandwidth Choice for Density Derivative"}else{
	                main <- "Biased Cross-Validation (2) function for \nBandwidth Choice for Density Function"}
                    }
	                }
    if(is.null(sub)) sub <- paste("Kernel",kernel,";","Derivative order = ",r)					
    if(is.null(seq.bws)){
       hos <- ((243 *(2*r+1)*A3_kMr(kernel,r))/(35* A2_kM(kernel)^2))^(1/(2*r+5)) * sd(x,na.rm = TRUE) * n^(-1/(2*r+5))
       seq.bws <- seq(0.15*hos,2*hos,length=50)
                         }
    R_Kr1 <- A3_kMr(kernel,r)
    R_Kr2 <- A3_kMr(kernel,r+2)
      fbcv1 <- function(h)
              {
         D1 <- kernel_fun_conv(kernel,outer(x,x,"-")/h,deriv.order=r+2)
         diag(D1) <- 0
         D1 <- ((-1)^(r+2)/((n-1)*h^(2*r+5)))* colSums(D1)
         D <- mean(D1)
         ##(1/(n*h^(2*r+1)))* R_Kr1 + (0.25*h^4)*(A2_kM(kernel))^2 * (D - (1/((n-1)*h^(2*r+5))) * R_Kr2 )
		 (1/(n*h^(2*r+1)))* R_Kr1 + (0.25*h^4)*(A2_kM(kernel))^2 * D
              }
      fbcv2 <- function(h)
              {
        D1 <- kernel_fun_der(kernel,outer(x,x,"-")/h,deriv.order=2*(r+2))
        diag(D1) <- 0
        D1 <- ((-1)^(r+2)/((n-1)*h^(2*r+5)))* colSums(D1)
        D <- mean(D1)
       (1/(n*h^(2*r+1)))* R_Kr1 + (0.25*h^4)*(A2_kM(kernel))^2 *  D
              }
    if(f$whichbcv == 1){
       D <- lapply(1:length(seq.bws), function(i) fbcv1(seq.bws[i]))}else{
       D <- lapply(1:length(seq.bws), function(i) fbcv2(seq.bws[i]))}
    Minf <- c(do.call("rbind",D))
    plot.default(seq.bws,Minf,type=type,las=las,lwd=lwd,xlab=xlab,ylab=ylab,
		         main=main,sub=sub,font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,...)
    return(list(kernel=kernel,deriv.order=r,seq.bws=seq.bws, bcv=Minf))
}

plot.h.bcv <- function(x,seq.bws=NULL,...) plot.bcv(x,seq.bws,...)

#####

lines.bcv <- function(f,seq.bws=NULL,...)
                    {
    class(f) <- "h.bcv"
    r <- f$deriv.order
    n <- f$n
    kernel <- f$kernel
    x <- sort(f$x)
    if(f$whichbcv == 1){
    if (kernel=="epanechnikov" && r+2 >= 3)    stop(" 'epanechnikov kernel derivative = 0' for 'order + 2 >= 3' ")
    else if (kernel=="triweight" && r+2 >= 7)  stop(" 'triweight kernel derivative = 0' for 'order + 2 >= 7' ")
    else if (kernel=="biweight" && r+2 >= 5)   stop(" 'biweight kernel derivative = 0' for 'order + 2 >= 5' ")
    else if (kernel=="tricube" && r+2 >= 10)   stop(" 'tricube kernel derivative = 0' for 'order + 2 >= 10' ")
    } else if (f$whichbcv == 2){
    if (kernel=="epanechnikov" && 2*(r+2) >= 3)    stop(" 'epanechnikov kernel derivative = 0' for '2 * (order + 2) >= 3' ")
    else if (kernel=="triweight" && 2*(r+2) >= 7)  stop(" 'triweight kernel derivative = 0' for '2 * (order + 2) >= 7' ")
    else if (kernel=="biweight" && 2*(r+2) >= 5)   stop(" 'biweight kernel derivative = 0' for '2 * (order + 2) >= 5' ")
    else if (kernel=="tricube" && 2*(r+2) >= 10)   stop(" 'tricube kernel derivative = 0' for '2 * (order + 2) >= 10' ")
    }
    if(is.null(seq.bws)){
       hos <- ((243 *(2*r+1)*A3_kMr(kernel,r))/(35* A2_kM(kernel)^2))^(1/(2*r+5)) * sd(x,na.rm = TRUE) * n^(-1/(2*r+5))
       seq.bws <- seq(0.15*hos,2*hos,length=50)
                         }
    R_Kr1 <- A3_kMr(kernel,r)
    R_Kr2 <- A3_kMr(kernel,r+2)
      fbcv1 <- function(h)
              {
         D1 <- kernel_fun_conv(kernel,outer(x,x,"-")/h,deriv.order=r+2)
         diag(D1) <- 0
         D1 <- ((-1)^(r+2)/((n-1)*h^(2*r+5)))* colSums(D1)
         D <- mean(D1)
         ##(1/(n*h^(2*r+1)))* R_Kr1 + (0.25*h^4)*(A2_kM(kernel))^2 * (D - (1/((n-1)*h^(2*r+5))) * R_Kr2 )
		 (1/(n*h^(2*r+1)))* R_Kr1 + (0.25*h^4)*(A2_kM(kernel))^2 * D
              }
      fbcv2 <- function(h)
              {
        D1 <- kernel_fun_der(kernel,outer(x,x,"-")/h,deriv.order=2*(r+2))
        diag(D1) <- 0
        D1 <- ((-1)^(r+2)/((n-1)*h^(2*r+5)))* colSums(D1)
        D <- mean(D1)
       (1/(n*h^(2*r+1)))* R_Kr1 + (0.25*h^4)*(A2_kM(kernel))^2 *  D
              }
    if(f$whichbcv == 1){
       D <- lapply(1:length(seq.bws), function(i) fbcv1(seq.bws[i]))}else{
       D <- lapply(1:length(seq.bws), function(i) fbcv2(seq.bws[i]))}
    Minf <- c(do.call("rbind",D))
    lines.default(seq.bws,Minf,...)
    invisible(NULL)
}

lines.h.bcv <- function(x,seq.bws=NULL,...) lines.bcv(x,seq.bws,...) 
