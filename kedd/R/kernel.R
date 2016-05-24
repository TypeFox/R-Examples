## Mon Jun 03 00:13:35 2013
## Original file Copyright 2013 A.C. Guidoum
## This file is part of the R package kedd.
## Arsalane Chouaib GUIDOUM <acguidoum@usthb.dz> and <starsalane@gmail.com> 
## Department of Probabilities-Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algeris
## Algeria
##############################################################################


## Kernels functions 


kernel.fun <- function(x,...) UseMethod("kernel.fun")

kernel.fun.default <- function(x=NULL,deriv.order=0,kernel=c("gaussian","epanechnikov",
                                  "uniform","triangular","triweight","tricube", 
                                  "biweight","cosine","silverman"),...)
              {
     if (any(deriv.order < 0 || deriv.order != round(deriv.order))) 
             stop("argument 'deriv.order' is non-negative integers")
     r <- deriv.order
     if (missing(kernel)) kernel <- "gaussian"
     if (is.null(x)){
         if (kernel == "gaussian"){x <- seq(-5,5,length=1000)} 
            else if (kernel == "silverman"){x <- seq(-10,10,length=1000)}
                  else {x <- seq(-1.5,1.5,length=1000)}
        }
     kx <- kernel_fun_der(kernel,u=x,deriv.order=r)
     structure(list(kernel = kernel,deriv.order=r,x=x,kx=kx),class="kernel.fun")
}

##############
##############

kernel.conv <- function(x,...) UseMethod("kernel.conv")

kernel.conv.default <- function(x=NULL,deriv.order=0,kernel=c("gaussian","epanechnikov",
                                "uniform","triangular","triweight","tricube", 
                                "biweight","cosine","silverman"),...)
              {
     if (any(deriv.order < 0 || deriv.order != round(deriv.order))) 
        stop("argument 'deriv.order' is non-negative integers")
     r <- deriv.order
     if (missing(kernel)) kernel <- "gaussian"
     if (is.null(x)){
         if (kernel == "gaussian"){x <- seq(-8,8,length=1000)} 
            else if (kernel == "silverman"){x <- seq(-10,10,length=1000)}
                  else {x <- seq(-2.5,2.5,length=1000)}
        }
     kx <- kernel_fun_conv(kernel,u=x,deriv.order=r)
     structure(list(kernel = kernel,deriv.order=r,x=x,kx=kx),class="kernel.conv")
}

#############
#############

plot.kernel.fun1d <- function(f,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                      type="l",las=1,lwd=1,...)
                  {
    class(f) <- "kernel.fun"
    r <- f$deriv.order
    kernel <- f$kernel
    if(is.null(xlab)) xlab <- "x"
    if(is.null(ylab)) ylab <- "" 
    if(is.null(main)){ 
    	     if(r != 0) {main <- paste("Derivative of ",kernel,"kernel")}else{
    	                 main <- paste(kernel,"kernel")}
    	                }               
    if(is.null(sub)){ 
    	     if(r != 0) {sub <- paste("Derivative order = ",r)}
    	                }
    plot.default(f$x,f$kx,type=type,las=las,lwd=lwd,xlab=xlab,ylab=ylab,
		       main=main,sub=sub,font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,...)
}

plot.kernel.fun <- function(x,...) plot.kernel.fun1d (x,...)


################################
################################

plot.kernel.conv1d <- function(f,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                             type="l",las=1,lwd=1,...)
                  {
    class(f) <- "kernel.conv"
    r <- f$deriv.order
    kernel <- f$kernel
    if(is.null(xlab)) xlab <- "x"
    if(is.null(ylab)) ylab <- ""             
    if(is.null(main)){ 
    	     if(r != 0) {main <- paste("Convolution of derivative",kernel,"kernel")}else{
    	                 main <- paste("Convolution of",kernel,"kernel")}
    	                }
    if(is.null(sub)){
	     if(r !=0) {sub <- paste("Derivative order = ",r)}
	                } 
    plot.default(f$x,f$kx,type=type,las=las,lwd=lwd,xlab=xlab,ylab=ylab,
		       main=main,sub=sub,font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,...)
}

plot.kernel.conv <- function(x,...) plot.kernel.conv1d (x,...)
