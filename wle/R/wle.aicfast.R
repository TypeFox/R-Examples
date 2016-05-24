#############################################################
#                                                           #
#	wle.aicfast function                                #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: July, 7, 2011                                 #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2011 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.aicfast <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, boot=30, group, var.full=0, num.sol=1, raf="HD", smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, min.weight=0.5, alpha=2, contrasts=NULL, verbose=FALSE) {

raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

if (raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

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

  z <- .Fortran("wleaicfast",
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
	as.double(alpha),
	waic=mat.or.vec(nrep,nvar+1),
	param=mat.or.vec(num.sol,nvar),
	var=double(num.sol),
	resid=mat.or.vec(num.sol,size),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
 	same=integer(num.sol),
	info=integer(1),
	PACKAGE="wle")

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
