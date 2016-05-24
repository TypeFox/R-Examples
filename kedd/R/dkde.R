## Tue Feb 26 00:39:10 2013
## Original file Copyright 2013 A.C. Guidoum
## This file is part of the R package kedd.
## Arsalane Chouaib GUIDOUM <acguidoum@usthb.dz> and <starsalane@gmail.com> 
## Department of Probabilities-Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algeris
## Algeria
##############################################################################


## Derivatives of Kernel Density Estimator (dkde)


dkde <- function(x, ...)  UseMethod("dkde")


dkde.default <- function(x,y=NULL,deriv.order=0,h,
                  kernel=c("gaussian","epanechnikov","uniform","triangular",
                  "triweight","tricube", "biweight","cosine"),...)
       {
    if (!is.numeric(x) || length(dim(x)) >=1 || length(x) < 3L) 
           stop("argument 'x' must be numeric and need at least 3 data points") 
    if (any(deriv.order < 0 || deriv.order != round(deriv.order))) 
             stop("argument 'deriv.order' is non-negative integers")
    r <- deriv.order
    if (kernel=="epanechnikov" && r >= 3) stop(" 'epanechnikov kernel derivative = 0' for 'order >= 3' ")
    if (kernel=="uniform" && r >= 1)      stop(" 'uniform kernel derivative = 0' for 'order >= 1' ")
    if (kernel=="triweight" && r >= 7)    stop(" 'triweight kernel derivative = 0' for 'order >= 7' ")
    if (kernel=="biweight" && r >= 5)     stop(" 'biweight kernel derivative = 0' for 'order >= 5' ")
    if (kernel=="triangular" && r >= 2)   stop(" 'triangular kernel derivative = 0' for 'order >= 2' ")
    if (kernel=="tricube" && r >= 10)     stop(" 'tricube kernel derivative = 0' for 'order >= 10' ")
    name <- deparse(substitute(x))
    if (missing(r))           r <- 0
    if (missing(kernel)) kernel <- "gaussian"
    if (missing(h))      h <- h.ucv(x,deriv.order=r,kernel=kernel)$h
	x <- x[!is.na(x)]
    if (is.null(y)){ 
               #from <- min(x) - 3 * h
               #to   <- max(x) + 3 * h
               #y <- seq(from - 4 * h, to + 4 * h,length.out=512L)	
               #tau <- if (kernel == "gaussian") 4 else 1
               tau <- 4			   
               range.x <- c(min(x) - tau * h, max(x) + tau * h)
               y <- seq(range.x[1L], range.x[2L],length.out=512L)
                  }
    aux <- outer(y,x,"-")/h
    fn <- kernel_fun_der(kernel, aux,r)/h^(r+1)
    fn <- apply(fn,1,mean)
    structure(list(x=x, data.name=name, n=length(x), kernel=kernel, deriv.order=r, h = h, 
                  eval.points=y, est.fx = fn),class="dkde")
         }
		 
###### 

print.dkde <- function(x, digits=NULL, ...)
           {
    class(x) <- "dkde"
    cat("\nData: ",x$data.name," (",x$n," obs.);", "\tKernel: ",x$kernel,"\n",
	    "\nDerivative order: ",x$deriv.order,";","\tBandwidth 'h' = ",formatC(x$h,digits=digits), "\n\n",sep="")
    print(summary(as.data.frame(x[c("eval.points","est.fx")])), digits=digits, ...)
    invisible(x)
}

######

plot.dkde1d <- function(f,fx=NULL,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                      type="l",las=1,lwd=1,col="red",lty=1,...)
                    {
    class(f) <- "dkde"
	#if(!is.function(fx)) stop("fx must be a function.")
    if(is.null(xlab)) xlab <- "x"
    if(is.null(ylab)){
	     if(f$deriv.order !=0){ylab <- "density derivative function"}else{
	                           ylab <- "density function"}
	                }
    if(is.null(main)){ 
	     if(f$deriv.order !=0) {main <- "Kernel density derivative estimate"}else{
	                            main <- "Kernel density estimation"}
	                }
    if(is.null(sub)){
	     if(f$deriv.order !=0) {sub <- paste("Kernel",f$kernel,";","Derivative order = ",f$deriv.order,";", "Bandwidth = ", formatC(f$h))}else{
	                            sub <- paste("Kernel",f$kernel,";", "Bandwidth = ", formatC(f$h))}
	                }  
    gn <- if (!is.null(fx)) function(par) fx(par, ...)
    if (is.null(fx)){
    plot.default(f$eval.points,f$est.fx, main=main,sub=sub, xlab=xlab, ylab=ylab, type=type,font.lab=2,cex.lab=0.9,
                 font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,las=las,lwd=lwd,...)}else{
    plot.default(f$eval.points,f$est.fx, ylim=c(min(f$est.fx,fx(f$eval.points),na.rm = TRUE),max(f$est.fx,fx(f$eval.points),na.rm = TRUE)),
	             main=main,sub=sub, xlab=xlab, ylab=ylab, type=type,col=col,lty=lty,
                 font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,las=las,lwd=lwd,...)
    curve(fx,xlim = c(min(f$eval.points,na.rm = TRUE), max(f$eval.points,na.rm = TRUE)), n = length(f$eval.points), lty = 8,lwd=lwd,add=TRUE)			 
    legend("topright", c("Estimate","True"),lty=c(lty,8),col=c(col,"black"),lwd=c(lwd,lwd), inset = .015,cex=1)
	}
    invisible(NULL)
}


plot.dkde <- function(x,fx=NULL,...) plot.dkde1d(x,fx,...)

######

lines.dkde1d <- function(f,...)
                    {
    class(f) <- "dkde"
    lines.default(f$eval.points,f$est.fx,...)
    invisible(NULL)
}

lines.dkde <- function(x,...) lines.dkde1d(x,...) 
