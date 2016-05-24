#  *****************************************************************************
#   File : kerfon.R
#         ************************************************************
#   Description : 
#       Kernel funcional modelization of autoregressive processes
#       Modelization, cross-validation and prediction
#   Version : 1.0
#   Date : 2001-07-06
#         ************************************************************
#   Author : Julien Damon <julien.damon@gmail.com>
#   License : LGPL
#   URL: https://github.com/Looping027/far
#  *****************************************************************************

#  *****************************************************************************
#   Title : kerfon
#         ************************************************************
#   Description : 
#       Modelization of Vectorized Functional Processes
#   Version : 1.0
#   Date : 2001-07-06
#  *****************************************************************************
kerfon <- function(data, x, r, hmin, hmax, na.rm=TRUE)
{
    if ( (!is.null(class(data)))
        && (class((data)) != "fdata")) # test the type of data
        stop("data is not of class fdata")
    call <- match.call()

    if (missing(x)) {x <- names(data)}
    if (length(x) > 1) x<-x[1]
    
    # Find dimensions
    p <- nrow(data[[x]])
    n <- ncol(data[[x]])

    # Removing non available data if required
    if (na.rm) {
        listobs <- c((!is.na(data))[x,])
        listobs2 <- c(FALSE,listobs) * c(listobs,FALSE) == 1
    } else {
        listobs <- rep(TRUE,n)
        listobs2 <- c(FALSE,rep(TRUE,n-1),FALSE)
    }
    nbobs <- sum(listobs == TRUE)
    nbobs2 <- sum(listobs2 == TRUE)
    
   	if (missing(r))
      	r<-ceiling(nbobs2/10)
      
    # Create Data
   	xdata <- (data[[x]])[,listobs2[-1],drop=FALSE]
   	ydata <- (data[[x]])[,listobs2[-(n+1)],drop=FALSE]
   	xtemp <- c(xdata)
   	ytemp <- c(ydata)
   	
   	varx<-var(xtemp)
   	if (missing(hmin))
      	hmin<-sqrt(varx)/8
   	if (missing(hmax))
      	hmax<-sqrt(varx)*4

   	f <- function(h,x,y,n,p,r)
   	{
	.C("CVkerfon",
		as.double(x),
		as.double(y),
		as.double(h),
		as.integer(n),
		as.integer(p),
		as.integer(r),
		result=double(1),
		PACKAGE="far")$result
   	}

   	if (hmin==hmax)
   	{
   		hopt <- hmin
   	} else {
   		hopt <- optimize(f,c(hmin,hmax),tol=0.1,x=xtemp,
                                 y=ytemp,n=nbobs2,p=p,r=r)
   	}
   	
   	output <- list(
		        "call" = call,
   				"h" = c(hopt[[1]],hmin,hmax),
   				"x" = x,
   				"xdata" = xdata,
   				"ydata" = ydata)
   class(output) <- "kerfon"
   invisible(output)
}

#  *****************************************************************************
#   Title : print.kerfon
#         ************************************************************
#   Description : 
#       Printing of class model "kerfon"
#   Version : 1.0
#   Date : 2001-06-07
#  *****************************************************************************
print.kerfon <- function(x, ..., digits = max(3, getOption("digits") - 3),
            na.print = "", file="", append=TRUE)
{
    cat("Functional Kernel Model\n",file=file,append=append)
    cat("Call: ", deparse(x$call), "\n\n",file=file,append=append)
    cat("Window: ",format(x$h[1], digits = digits),"optimized between ",
    	format(x$h[2], digits = digits)," and ",
    	format(x$h[3], digits = digits),".\n",file=file,append=append)
}

#  *****************************************************************************
#   Title : predict.kerfon
#         ************************************************************
#   Description : 
#       Computation of prediction for the class model "kerfon"
#   Version : 1.0
#   Date : 2001-06-07
#  *****************************************************************************
predict.kerfon <- function(object, ..., newdata = NULL, label, na.rm=TRUE,
                           positive=FALSE)
{
    if ( (!is.null(class(object)))
        && (class((object)) != "kerfon")) # test the type of data
        stop("object is not of class kerfon")
    if ( (!is.null(class(newdata)))
        && (class((newdata)) != "fdata")) # test the type of data
        stop("newdata is not of class fdata")

   	h <- object$h[1]
   	p <- dim(object$xdata)[1]
   	n <- dim(object$xdata)[2]
    if (missing(label))
        label <- c(colnames(newdata[[object$x]])[-1],
                   paste(ncol(newdata[[object$x]])+1))
    else
        label <- c(colnames(newdata[[object$x]])[-1],label)

    if (na.rm)
        listobs <- c((!is.na(newdata))[object$x,])
    else
        listobs <- rep(TRUE,n)
    nbobs <- sum(listobs==TRUE)

    pred <- matrix(0,nrow=p,ncol=nbobs)
    dimnames(pred) <- list(rownames(newdata[[object$x]]),label[listobs])

	newdata <- (newdata[[object$x]])[,listobs,drop=FALSE]
   	for (i in 1:nbobs)
   	{
		pred[,i]<-.C("prevkerfon",
			as.double(c(object$xdata)),
			as.double(c(object$ydata)),
			as.double(newdata[,i]),
			as.double(h),
			as.integer(n),
			as.integer(p),
			result=double(p),
			PACKAGE="far")$result
	}
    if (positive)
    {
         pred <- (pred+abs(pred))/2
    }
    return(as.fdata(pred,name=object$x))
}
