#############################################################
#                                                           #
#	wle.onestep function                                #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: December, 19, 2000                            #
#	Version: 0.3                                        #
#                                                           #
#	Copyright (C) 2000 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.onestep <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, ini.param, ini.scale, raf="HD", smooth=0.031, num.step=1, contrasts=NULL, verbose=FALSE) {

raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

if (raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

    ret.x <- x
    ret.y <- y
    result <- list()	
    mt <- terms(formula, data = data)
    mf <- cl <- match.call()
    mf$ini.param <- mf$ini.scale <- mf$smooth <- NULL
    mf$num.step <- mf$raf <- mf$contrasts <- NULL
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

if (size<(nvar+1)) {stop("Number of observation must be at least equal to the number of predictors (including intercept) + 1")}

if (!(ini.scale>=0)) {
stop("The initial scale error must be non negative")
}

if (!(num.step>=1)) {
    if (verbose) cat("wle.onestep: number of steps can not be negative, set to 1 \n")
    num.step <- 1
}

if (smooth<10^(-5)) {
    if (verbose) cat("wle.onestep: the smooth parameter seems too small \n")
}

ini.var <- ini.scale^2

  z <- .Fortran("wleonestepfix",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(0), 
	as.integer(size),
	as.integer(nvar),
	as.integer(nvar),
	as.double(ini.param),
	as.double(ini.var),
	as.integer(raf),
	as.double(smooth),
	as.integer(num.step),
	param=double(nvar),
	var=double(1),
	resid=double(size),
	totweight=double(1),
	weight=double(size),
	density=double(size),
	model=double(size),
	delta=double(size),
        PACKAGE = "wle")

if (z$var>0) {

    devparam <- sqrt(z$var*diag(solve(t(xdata)%*%diag(z$weight)%*%xdata)))

    result$coefficients <- z$param
    result$standard.error <- devparam
    result$scale <- sqrt(z$var)
    result$residuals <- z$resid
    result$fitted.values <- as.vector(xdata%*%z$param)
    result$weights <- z$weight
    result$f.density <- z$density
    result$m.density <- z$model
    result$delta <- z$delta
    result$tot.weights <- z$totweight

} else {
    if (verbose) cat("The initial estimates do not seems very good: the total sum of the weights is less than number of independent variables\n")


    result$coefficients <- rep(NA,nvar)
    result$standard.error <- rep(NA,nvar)
    result$scale <- NA
    result$residuals <- rep(NA,size)
    result$fitted.values <- rep(NA,size)
    result$weights <- rep(NA,size)
    result$f.density <- rep(NA,size)
    result$m.density <- rep(NA,size)
    result$delta <- rep(NA,size)
    result$tot.weights <- NA
}

result$call <- cl
result$contrasts <- attr(xdata, "contrasts")
result$xlevels <- xlev
result$terms <- mt

if (model)
    result$model <- mf
if (ret.x)
    result$x <- xdata
if (ret.y)
    result$y <- ydata

if (is.null(names(ini.param))) {
    dn <- colnames(xdata)
} else {
    dn <- names(ini.param)
}

if (is.null(nrow(result$coefficients))) {
    names(result$coefficients) <- dn
} else {
    dimnames(result$coefficients) <- list(NULL,dn)
}

class(result) <- "wle.onestep"
return(result)
}

#############################################################
#                                                           #
#	print.wle.onestep function                          #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: December, 19, 2000                            #
#	Version: 0.3                                        #
#                                                           #
#	Copyright (C) 2000 Claudio Agostinelli              #
#                                                           #
#############################################################

print.wle.onestep <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Coefficients:\n")
    print.default(format(coef(x), digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Scale estimate: ",format(x$scale, digits=digits))
    cat("\n")

    invisible(x)
}
