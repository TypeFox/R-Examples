#############################################################
#                                                           #
#	wle.cv function                                     #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: April, 02, 2002                               #
#	Version: 0.4-1                                      #
#                                                           #
#	Copyright (C) 2002 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.cv <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, monte.carlo=500, split, boot=30, group, num.sol=1, raf="HD", smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, min.weight=0.5, contrasts=NULL, type=c('fast', 'slow'), verbose=FALSE) {

type <- match.arg(type)
  
raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

if (raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

if (missing(group)) {
group <- 0
}

if (missing(split)) {
split <- 0
}

    ret.x <- x
    ret.y <- y
    result <- list()	
    mt <- terms(formula, data = data)
    mf <- cl <- match.call()
    mf$monte.carlo <- mf$split <- NULL
    mf$boot <- mf$group <- mf$smooth <- NULL
    mf$tol <- mf$equal <- mf$num.sol <- NULL
    mf$min.weight <- mf$max.iter <- mf$raf <- NULL
    mf$contrasts <- NULL
    mf$model <- mf$x <- mf$y <- NULL
    mf$type <- NULL
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
    if (verbose) cat("wle.cv: dimension of the subsample set to default value = ",group,"\n")
}

maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))

if (boot<1 | log(boot) > maxboot) {
    stop("Bootstrap replication not in the range")
}

if (split<nvar+2 | split>(size-2)) {
    split <- max(round(size^(3/4)),nvar+2)
    if (verbose) cat("wle.cv: dimension of the split subsample set to default value = ",split,"\n")
}

maxcarlo <- sum(log(1:size))-(sum(log(1:split))+sum(log(1:(size-split))))

if (monte.carlo<1 | log(monte.carlo) > maxcarlo) {
    stop("MonteCarlo replication not in the range")
}

if (!(num.sol>=1)) {
    if (verbose) cat("wle.cv:number of solution to report set to 1 \n")
    num.sol <- 1
}

if (max.iter<1) {
    if (verbose) cat("wle.cv: max number of iteration set to 500 \n")
    max.iter <- 500
}

if (smooth<10^(-5)) {
    if (verbose) cat("wle.cv: the smooth parameter seems too small \n")
}

if (tol<=0) {
    if (verbose) cat("wle.cv: the accuracy must be positive, using default value: 10^(-6) \n")
    tol <- 10^(-6)
}

if (equal<=tol) {
    if (verbose) cat("wle.cv: the equal parameter must be greater than tol, using default value: tol+10^(-3) \n")
    equal <- tol+10^(-3)
}

if (min.weight<0) {
    if (verbose) cat("wle.cv: the minimum sum of the weights can not be negative, using default value \n")
    min.weight <- 0.5
}

if (type=='fast') {
  z <- .Fortran("wlecv",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(0), 
	as.integer(size),
	as.integer(nvar),
	as.integer(boot),
	as.integer(group),
	as.integer(nrep),
	as.integer(monte.carlo),
	as.integer(split),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	as.integer(num.sol),
	as.double(min.weight),
	wcv=mat.or.vec(nrep,nvar+1),
	param=mat.or.vec(num.sol,nvar),
	var=double(num.sol),
	resid=mat.or.vec(num.sol,size),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
	same=integer(num.sol),
	index=integer(1),
	info=integer(1),
	PACKAGE="wle")
} else {
  z <- .Fortran("wwlecv",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(0), 
	as.integer(size),
	as.integer(nvar),
	as.integer(boot),
	as.integer(group),
	as.integer(nrep),
	as.integer(monte.carlo),
	as.integer(split),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	as.integer(num.sol),
	as.double(min.weight),
	wcv=mat.or.vec(nrep,nvar+1),
	param=mat.or.vec(num.sol,nvar),
	var=double(num.sol),
	resid=mat.or.vec(num.sol,size),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
	same=integer(num.sol),
	index=integer(1),
	info=integer(1),
	PACKAGE="wle")
}

delnull <- z$same==0

result$wcv <- z$wcv
result$coefficients <- z$param[!delnull,]
result$scale <- sqrt(z$var[!delnull])
result$residuals <- z$resid[!delnull]
result$weights <- z$weight[!delnull,]
result$tot.weights <- z$totweight[!delnull]
result$freq <- z$same[!delnull]
result$call <- cl
result$info <- z$info
result$index <- z$index
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
if (is.null(nrow(result$coefficients))) {
names(result$coefficients) <- dn
} else {
dimnames(result$coefficients) <- list(NULL,dn)
}
dimnames(result$wcv) <- list(NULL,c(dn,"wcv"))

class(result) <- "wle.cv"

return(result)
}

#############################################################
#                                                           #
#	summary.wle.cv function                             #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 3, 2001                             #
#	Version: 0.4-1                                      #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

summary.wle.cv <- function (object, num.max=20, verbose=FALSE, ...) {

if (is.null(object$terms)) {
    stop("invalid \'wle.cv\' object")
}

if (num.max<1) {
    if (verbose) cat("summary.wle.cv: num.max can not less than 1, num.max=1 \n")
    num.max <- 1
}

ans <- list()
wcv <- object$wcv
if(is.null(nmodel <- nrow(wcv))) nmodel <- 1
num.max <- min(nmodel,num.max)
if (nmodel!=1) { 
nvar <- ncol(wcv)-1
wcv <- wcv[order(wcv[,(nvar+1)]),]
wcv <- wcv[1:num.max,]
}

ans$wcv <- wcv
ans$num.max <- num.max
ans$call <- object$call

class(ans) <- "summary.wle.cv"
return(ans)
}

#############################################################
#                                                           #
#	print.wle.cv function                                   #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: October, 27, 2003                                 #
#	Version: 0.4-1                                          #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.wle.cv <- function (x, digits = max(3, getOption("digits") - 3), num.max=max(1, nrow(x$wcv)), ...) {
    res <- summary.wle.cv(object=x, num.max=num.max, ...)
    print.summary.wle.cv(res, digits=digits, ...)
}

#############################################################
#                                                           #
#	print.summary.wle.cv function                       #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: August, 2, 2001                               #
#	Version: 0.4                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

print.summary.wle.cv <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")

    cat("\nWeighted Cross Validation selection criteria:\n")
    if(x$num.max>1) {
    nvar <- ncol(x$wcv)-1
    x$wcv[,(nvar+1)] <- signif(x$wcv[,(nvar+1)],digits)
    } else {
    nvar <- length(x$wcv)-1
    x$wcv[(nvar+1)] <- signif(x$wcv[(nvar+1)],digits)
    }
    print(x$wcv)
    cat("\n")

    cat("Printed the first ",x$num.max," best models \n") 
    invisible(x)
}



