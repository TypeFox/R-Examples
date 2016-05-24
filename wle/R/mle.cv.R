#############################################################
#                                                           #
#	mle.cv function                                         #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: April, 02, 2002                                   #
#	Version: 0.4                                            #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli                  #
#                                                           #
#############################################################

mle.cv <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, monte.carlo=500, split, contrasts=NULL, verbose=FALSE) {

if (missing(split)) {
    split <- 0
}

    ret.x <- x
    ret.y <- y
    result <- list()	
    mt <- terms(formula, data = data)
    mf <- cl <- match.call()
    mf$monte.carlo <- mf$split <- mf$contrasts <- NULL
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

if (size<nvar+1) {stop("Number of observation must be at least equal to the number of predictors (including intercept) + 1")}

if (split<nvar+2 | split>(size-2)) {
    split <- max(round(size^(3/4)),nvar+2)
    if (verbose) cat("mle.cv: dimension of the split subsample set to default value = ",split,"\n")
}
maxcarlo <- sum(log(1:size))-(sum(log(1:split))+sum(log(1:(size-split))))
if (monte.carlo<1 | log(monte.carlo) > maxcarlo){
    stop("MonteCarlo replication not in the range")
}

  z <- .Fortran("mlecv",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(0), 
	as.integer(size),
	as.integer(nvar),
	as.integer(nrep),
	as.integer(monte.carlo),
	as.integer(split),
	cv=mat.or.vec(nrep,nvar+1),
	info=integer(1),
	PACKAGE="wle")


result$cv <- z$cv
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
dimnames(result$cv) <- list(NULL,c(dn,"cv"))

class(result) <- "mle.cv" 

return(result)
}

#############################################################
#                                                           #
#	summary.mle.cv function                             #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: August, 2, 2001                               #
#	Version: 0.4-1                                      #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

summary.mle.cv <- function (object, num.max=20, verbose=FALSE, ...) {

    if (is.null(object$terms)) {
        stop("invalid \'mle.cv\' object")
    }

    if (num.max<1) {
        if (verbose) cat("summary.mle.cv: num.max can not less than 1, num.max=1 \n")
        num.max <- 1
    }

    ans <- list()
    cv <- object$cv
    if(is.null(nmodel <- nrow(cv))) nmodel <- 1
    num.max <- min(nmodel,num.max)
    if (nmodel!=1) { 
        nvar <- ncol(cv)-1
        cv <- cv[order(cv[,(nvar+1)]),]
        cv <- cv[1:num.max,]
    }

    ans$cv <- cv
    ans$num.max <- num.max
    ans$call <- object$call

    class(ans) <- "summary.mle.cv"
    return(ans)
}

#############################################################
#                                                           #
#	print.mle.cv function                                   #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: October, 27, 2003                                 #
#	Version: 0.4-1                                          #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.mle.cv <- function (x, digits = max(3, getOption("digits") - 3), num.max=max(1, nrow(x$cv)), ...) {
    res <- summary.mle.cv(object=x, num.max=num.max, ...)
    print.summary.mle.cv(res, digits=digits, ...)
}

#############################################################
#                                                           #
#	print.summary.mle.cv function                       #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: August, 2, 2001                               #
#	Version: 0.4                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

print.summary.mle.cv <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")

    cat("\nCross Validation selection criteria:\n")
    if(x$num.max>1) {
    nvar <- ncol(x$cv)-1
    x$cv[,(nvar+1)] <- signif(x$cv[,(nvar+1)],digits)
    } else {
    nvar <- length(x$cv)-1
    x$cv[(nvar+1)] <- signif(x$cv[(nvar+1)],digits)
    }
    print(x$cv)
    cat("\n")

    cat("Printed the first ",x$num.max," best models \n") 
    invisible(x)
}



