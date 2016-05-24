#############################################################
#                                                           #
#	mle.aic function                                    #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: November, 17, 2005                            #
#	Version: 0.4-1                                      #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################

mle.aic <- function(formula, data=list(), model=TRUE, x=FALSE,
                    y=FALSE, var.full=0, alpha=2, contrasts = NULL,
                    se=FALSE, verbose=FALSE) {

    ret.x <- x
    ret.y <- y
    result <- list()	
    mt <- terms(formula, data = data)
    mf <- cl <- match.call()
    mf$var.full <- mf$alpha <- mf$contrasts <- NULL
    mf$model <- mf$x <- mf$y <- NULL
    mf$se <- mf$verbose <- NULL
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

if (size<nvar+1) {stop("Number of observation must be at least equal to the number of predictors (including intercept) + 1")}
if (var.full<0) {
    if (verbose) cat("mle.aic: the variance of the full model can not be negative, using default value \n")
    var.full <- 0
}

  z <- .Fortran("mleaic",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(0), 
	as.integer(size),
	as.integer(nvar),
	as.integer(nrep),
	as.double(var.full),
	as.double(alpha),
	aic=mat.or.vec(nrep,nvar+1),
	param=mat.or.vec(nrep,nvar),
	var=double(nrep),
	resid=mat.or.vec(nrep,size),
	info=integer(1),
	PACKAGE="wle")
	
result$aic <- z$aic
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
dimnames(result$aic) <- list(NULL,c(dn,"aic"))

nrc <- dim(result$coefficients)

if (se){
    semat <- matrix(0, nrow=nrc[1], ncol=nrc[2])
  
    for (i in 1:nrow(result$aic)) {
         pos <- as.logical(result$aic[i,-ncol(result$aic)])
         xtemp <- xdata[,pos]
         setemp <- result$scale[i]*sqrt(diag(solve(t(xtemp)%*%xtemp)))
         semat[i, pos] <- setemp 
    }
    result$se <- semat
    dimnames(result$se) <- list(NULL,dn)
}
    

class(result) <- "mle.aic"

return(result)
}

#############################################################
#                                                           #
#	summary.mle.aic function                            #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 3, 2001                             #
#	Version: 0.4-1                                      #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

summary.mle.aic <- function (object, num.max=20, verbose=FALSE, ...) {

if (is.null(object$terms)) {
    stop("invalid \'mle.aic\' object")
}

if (num.max<1) {
    if (verbose) cat("summary.mle.aic: num.max can not less than 1, num.max=1 \n")
    num.max <- 1
}

ans <- list()
aic <- object$aic
if (is.null(nmodel <- nrow(aic))) nmodel <- 1
num.max <- min(nmodel,num.max)
if (nmodel!=1) { 
    nvar <- ncol(aic)-1
    aic <- aic[order(aic[,(nvar+1)]),]
    aic <- aic[1:num.max,]
}

ans$aic <- aic
ans$num.max <- num.max
ans$call <- object$call

class(ans) <- "summary.mle.aic"
return(ans)
}

#############################################################
#                                                           #
#	print.mle.aic function                                  #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: October, 27, 2003                                 #
#	Version: 0.4-1                                          #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.mle.aic <- function (x, digits = max(3, getOption("digits") - 3), num.max=max(1, nrow(x$aic)), ...) {
   res <- summary.mle.aic(object=x, num.max=num.max, ...)
   print.summary.mle.aic(res, digits=digits, ...)
}

#############################################################
#                                                           #
#	print.summary.mle.aic function                      #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: August, 2, 2001                               #
#	Version: 0.4                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

print.summary.mle.aic <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")

    cat("\nAkaike Information Criterion (AIC):\n")
    if(x$num.max>1) {
    nvar <- ncol(x$aic)-1
    x$aic[,(nvar+1)] <- signif(x$aic[,(nvar+1)],digits)
    } else {
    nvar <- length(x$aic)-1
    x$aic[(nvar+1)] <- signif(x$aic[(nvar+1)],digits)
    }
    print(x$aic)
    cat("\n")

    cat("Printed the first ",x$num.max," best models \n") 
    invisible(x)
}






