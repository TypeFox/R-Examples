#############################################################
#                                                           #
#	mle.cp function                                     #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: August, 2, 2001                               #
#	Version: 0.4                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

mle.cp <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, var.full=0, contrasts=NULL, verbose=FALSE) {

    ret.x <- x
    ret.y <- y
    result <- list()	
    mt <- terms(formula, data = data)
    mf <- cl <- match.call()
    mf$var.full <- mf$contrasts <- NULL
    mf$model <- mf$x <- mf$y <- NULL
    mf$verbose <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    xvars <- as.character(attr(mt, "variables"))[-1]
    inter <- attr(mt, "intercept")
    if((yvar <- attr(mt, "response")) > 0) xvars <- xvars[-yvar]
    xlev <-
	if(length(xvars) > 0) {
	    xlev <- lapply(mf[xvars], levels)
	    xlev[!sapply(xlev, is.null)]
	}
    ydata <- model.response(mf, "numeric")
    if (is.empty.model(mt)) 
	stop("The model is empty")
    else 
	xdata <- model.matrix(mt, mf, contrasts)

if (is.null(size <- nrow(xdata)) | is.null(nvar <- ncol(xdata))) stop("'x' must be a matrix")
if (length(ydata)!=size) stop("'y' and 'x' are not compatible")

nrep <- 2^nvar-1

if (size<nvar+1) {
stop("Number of observation must be at least equal to the number of predictors (including intercept) + 1")
}

if (var.full<0) {
    if (verbose) cat("mle.cp: the variance of the full model can not be negative, using default value \n")
    var.full <- 0
}

  z <- .Fortran("mlecp",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(0), 
	as.integer(size),
	as.integer(nvar),
	as.integer(nrep),
	as.double(var.full),
	cp=mat.or.vec(nrep,nvar+1),
	param=mat.or.vec(nrep,nvar),
	var=double(nrep),
	resid=mat.or.vec(nrep,size),
	info=integer(1),
	PACKAGE="wle")


result$cp <- z$cp
result$coefficients <- z$param
result$scale <- sqrt(z$var)
result$residuals <- z$resid
result$call <- cl
result$info <- z$info
result$contrasts <- attr(xdata, "contrasts")
result$xlevels <- xlev
result$terms <- mt

if (model)
    result$model <- mf
if (ret.x)
    result$x <- xdata
if (ret.y)
    result$y <- ydata

dn <- colnames(xdata)
dimnames(result$coefficients) <- list(NULL,dn)
dimnames(result$cp) <- list(NULL,c(dn,"cp"))

class(result) <- "mle.cp" 

return(result)

}

#############################################################
#                                                           #
#	summary.mle.cp function                             #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: October, 28, 2003                             #
#	Version: 0.4-1                                      #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#############################################################

summary.mle.cp <- function (object, num.max=20, verbose=FALSE, ...) {

if (is.null(object$terms)) {
    stop("invalid \'mle.cp\' object")
}

ans <- list()
cp <- object$cp

if (num.max<1) {
    if (verbose) cat("summary.mle.cp: num.max can not less than 1, num.max=1 \n")
    num.max <- 1
}

if(is.null(nmodel <- nrow(cp))) nmodel <- 1
num.max <- min(nmodel,num.max)
if (nmodel!=1) { 
    nvar <- ncol(cp)-1
    nparam <- apply(cp[,(1:nvar)],1,sum)
    cp <- cp[cp[,(nvar+1)]<=(nparam+0.00001),]
    if (!is.null(nrow(cp)) && nrow(cp)>1) {
	num.max <- min(nrow(cp),num.max)
    	cp <- cp[order(cp[,(nvar+1)]),]
    	cp <- cp[1:num.max,]
    } else num.max <- 1
}

ans$cp <- cp
ans$num.max <- num.max
ans$call <- object$call

class(ans) <- "summary.mle.cp"
return(ans)
}

#############################################################
#                                                           #
#	print.mle.cp function                                   #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: October, 27, 2003                                 #
#	Version: 0.4-1                                          #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.mle.cp <- function (x, digits = max(3, getOption("digits") - 3),  num.max=max(1, nrow(x$cp)), ...) {
    res <- summary.mle.cp(object=x, num.max=num.max, ...)
    print.summary.mle.cp(res, digits=digits, ...)
}

#############################################################
#                                                           #
#	print.summary.mle.cp function                       #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: August, 2, 2001                               #
#	Version: 0.4                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

print.summary.mle.cp <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")

    cat("\nMallows Cp:\n")
    if(x$num.max>1) {
    nvar <- ncol(x$cp)-1
    x$cp[,(nvar+1)] <- signif(x$cp[,(nvar+1)],digits)
    } else {
    nvar <- length(x$cp)-1
    x$cp[(nvar+1)] <- signif(x$cp[(nvar+1)],digits)
    }
    print(x$cp)
    cat("\n")

    cat("Printed the first ",x$num.max," best models \n") 
    invisible(x)
}


