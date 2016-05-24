#############################################################
#                                                           #
#	wle.aic function                                    #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: September, 28, 2010                           #
#	Version: 0.5                                        #
#                                                           #
#	Copyright (C) 2010 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.aic <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, boot=30, group, var.full=0, num.sol=1, raf="HD", smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, min.weight=0.5, method="full", alpha=2, contrasts=NULL, verbose=FALSE) {

raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

if (raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

type <- switch(method,
	full = 0,
	reduced = 1,
	-1)

if (type==-1) stop("Please, choose the method: full=wieghts based on full model, reduced=weights based on the actual model")

if (missing(group)) {
group <- 0
}

    ret.x <- x
    ret.y <- y
    result <- list()	
    mt <- terms(formula, data = data)
    mf <- cl <- match.call()
    mf$boot <- mf$group <- mf$smooth <- NULL
    mf$tol <- mf$equal <- mf$num.sol <- NULL
    mf$min.weight <- mf$max.iter <- mf$raf <- NULL
    mf$var.full <- mf$alpha <- mf$contrasts <- NULL
    mf$model <- mf$x <- mf$y <- mf$method <- NULL
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

if (size<nvar) {
    stop("Number of observations must be at least equal to the number of predictors (including intercept)")
}

if (group<nvar) {
    group <- max(round(size/4),nvar)
    if (verbose) cat("wle.aic: dimension of the subsample set to default value = ",group,"\n")
}

maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))

if (boot<1 | log(boot) > maxboot) {
    stop("bootstrap replication not in the range")
}

if (!(num.sol>=1)) {
    if (verbose) cat("wle.aic: number of solution to report set to 1 \n")
    num.sol <- 1
}

if (max.iter<1) {
    if (verbose) cat("wle.aic: max number of iteration set to 500 \n")
    max.iter <- 500
}

if (smooth<10^(-5)) {
    if (verbose) cat("wle.aic: the smooth parameter seems too small \n")
}

if (tol<=0) {
    if (verbose) cat("wle.aic: the accuracy must be positive, using default value: 10^(-6) \n")
    tol <- 10^(-6)
}

if (equal<=tol) {
    if (verbose) cat("wle.aic: the equal parameter must be greater than tol, using default value: tol+10^(-3) \n")
    equal <- tol+10^(-3)
}

if (var.full<0) {
    if (verbose) cat("wle.aic: the variance of the full model can not be negative, using default value \n")
    var.full <- 0
}

if (min.weight<0) {
    if (verbose) cat("wle.aic: the minimum sum of the weights can not be negative, using default value \n")
    min.weight <- 0.5
}

  z <- .Fortran("wleaic",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(0), 
	as.integer(size),
	as.integer(nvar),
	as.integer(boot),
	as.integer(group),
	as.integer(nrep),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	as.double(var.full),
	as.integer(num.sol),
	as.double(min.weight),
	as.integer(type),
	as.double(alpha),
	waic=mat.or.vec(nrep*num.sol,nvar+1),
	param=mat.or.vec(nrep*num.sol,nvar),
	var=double(nrep*num.sol),
	resid=mat.or.vec(nrep*num.sol,size),
	totweight=double(nrep*num.sol),
	weight=mat.or.vec(nrep*num.sol,size),
 	same=integer(nrep*num.sol),
	info=integer(1),
	PACKAGE="wle")

delnull <- z$same==0

z$waic[delnull,nvar+1] <- NA
z$param[delnull,] <- NA
z$var[delnull] <- NA
z$resid[delnull] <- NA
z$weight[delnull,] <- NA
z$totweight[delnull] <- NA

result$waic <- z$waic
result$coefficients <- z$param
result$scale <- sqrt(z$var)
result$residuals <- z$resid
result$weights <- z$weight
result$tot.weights <- z$totweight
result$freq <- z$same
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
dimnames(result$waic) <- list(NULL,c(dn,"waic"))

class(result) <- "wle.aic"

return(result)
}

#############################################################
#                                                           #
#	summary.wle.aic function                            #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: August, 2, 2001                               #
#	Version: 0.4-1                                      #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

summary.wle.aic <- function (object, num.max=20, verbose=FALSE, ...) {

if (is.null(object$terms)) {
    stop("invalid \'wle.aic\' object")
}

if (num.max<1) {
    if (verbose) cat("summary.wle.aic: num.max can not less than 1, num.max=1 \n")
    num.max <- 1
}

ans <- list()
waic <- object$waic
if (is.null(nmodel <- nrow(waic))) nmodel <- 1
num.max <- min(nmodel,num.max)
if (nmodel!=1) { 
    nvar <- ncol(waic)-1
    waic <- waic[order(waic[,(nvar+1)]),]
    waic <- waic[1:num.max,]
}

ans$waic <- waic
ans$num.max <- num.max
ans$call <- object$call

class(ans) <- "summary.wle.aic"
return(ans)
}

#############################################################
#                                                           #
#	print.wle.aic function                                  #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: October, 27, 2003                                 #
#	Version: 0.4-1                                          #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.wle.aic <- function (x, digits = max(3, getOption("digits") - 3), num.max=max(1, nrow(x$waic)), ...) {
    res <- summary.wle.aic(object=x, num.max=num.max, ...)
    print.summary.wle.aic(res, digits=digits, ...)
}

#############################################################
#                                                           #
#	print.summary.wle.aic function                      #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: August, 2, 2001                               #
#	Version: 0.4                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

print.summary.wle.aic <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")

    cat("\nWeighted Akaike Information Criterion (WAIC):\n")
    if(x$num.max>1) {
    nvar <- ncol(x$waic)-1
    x$waic[,(nvar+1)] <- signif(x$waic[,(nvar+1)],digits)
    } else {
    nvar <- length(x$waic)-1
    x$waic[(nvar+1)] <- signif(x$waic[(nvar+1)],digits)
    }
    print(x$waic)
    cat("\n")

    cat("Printed the first ",x$num.max," best models \n") 
    invisible(x)
}








