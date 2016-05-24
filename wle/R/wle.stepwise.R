#############################################################
#                                                           #
#	wle.stepwise function                               #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: February, 10, 2010                            #
#	Version: 0.6                                        #
#                                                           #
#	Copyright (C) 2010 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.stepwise <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, boot=30, group, num.sol=1, raf="HD", smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, min.weight=0.5, type="Forward", f.in=4.0, f.out=4.0, method="WLE", contrasts=NULL, verbose=FALSE)
{

  raf <- switch(raf,
    HD = 1,
    NED = 2,
    SCHI2 = 3,
    -1)

  if (raf==-1)
    stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

  ntype <- switch(type,
    Forward = 1,
    Backward = 2,
    Stepwise = 3,
    -1)

  if (ntype==-1)
    stop("The type must be Forward, Backward or Stepwise")

  nmethod <- switch(method,
    WLE = 0,
    WLS = 1,
    -1)

  if (nmethod==-1)
    stop("The method must be WLE, or WLS the default value is WLE")

  if (missing(group))
    group <- 0

  ret.x <- x
  ret.y <- y
  result <- list()	
  mt <- terms(formula, data = data)
  mf <- cl <- match.call()
  mf$boot <- mf$group <- mf$smooth <- NULL
  mf$tol <- mf$equal <- mf$num.sol <- NULL
  mf$max.iter <- mf$raf <- mf$contrasts <- NULL
  mf$min.weight <- NULL
  mf$type <- mf$f.in <- mf$f.out <- NULL
  mf$model <- mf$x <- mf$y <- mf$method <- NULL
  mf$verbose <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent()))
  xvars <- as.character(attr(mt, "variables"))[-1]
  inter <- attr(mt, "intercept")
  if((yvar <- attr(mt, "response")) > 0) xvars <- xvars[-yvar]
  xlev <- if(length(xvars) > 0) {
    xlev <- lapply(mf[xvars], levels)
    xlev[!sapply(xlev, is.null)]
  }
  ydata <- model.response(mf, "numeric")
  if (is.empty.model(mt)) 
    stop("The model is empty")
  else 
    xdata <- model.matrix(mt, mf, contrasts)
  if (inter>1) {
      xdata <- cbind(xdata[,inter, drop=FALSE], xdata[,-inter])
  } else if(inter==0 & ntype!=2) {
    warning("An intercept term is inserted in the model, please do not insert an intercept in the response yourself")
    xdata <- cbind(rep(1, NROW(xdata)), xdata)
    colnames(xdata)[1] <- "(Intercept)"
  }
  if (is.null(size <- nrow(xdata)) | is.null(nvar <- ncol(xdata)))
    stop("'x' must be a matrix")
  if (length(ydata)!=size)
    stop("'y' and 'x' are not compatible")

  if (size<nvar+1)
    stop("Number of observations must be at least equal to the number of predictors (including intercept) + 1")
  if (f.in<0 | f.out<0)
    stop("f.in and f.out can not be negative")

  if (group < 1) {
    group <- max(round(size/4),nvar+1)
    if (verbose) cat("wle.stepwise: dimension of the subsample set to default value: ",group,"\n")
  }
  maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))
  if (boot<1 | log(boot) > maxboot)
    stop("Bootstrap replication not in the range")
  if (!(num.sol>=1)) {
    if (verbose)
      cat("wle.stepwise: number of solution to report set to 1 \n")
    num.sol <- 1
  }
  if (max.iter<1) {
    if (verbose)
      cat("wle.stepwise: max number of iteration set to 500 \n")
    max.iter <- 500
  }
  if (smooth<10^(-5)) {
    if (verbose)
      cat("wle.stepwise: the smooth parameter seems too small \n")
  }

  if (tol<=0) {
    if (verbose)
      cat("wle.stepwise: the accuracy must be positive, using default value: 10^(-6)\n")
    tol <- 10^(-6)
  }

  if (equal<=tol) {
    if (verbose)
      cat("wle.stepwise: the equal parameter must be greater than tol, using default value: tol+10^(-3)\n")
    equal <- tol+10^(-3)
  }

  if (min.weight<0) {
    if (verbose)
      cat("wle.stepwise: the minimum sum of the weights can not be negative, using default value \n")
    min.weight <- 0.5
  }

  nrep <- 2^nvar-1

  z <- .Fortran("wstep",
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
	as.integer(ntype),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	as.integer(num.sol),
	as.double(min.weight),
	as.double(f.in),
	as.double(f.out),
	as.integer(nmethod),
	wstep=mat.or.vec(nrep,nvar+1),
	param=mat.or.vec(num.sol,nvar),
	var=double(num.sol),
	resid=mat.or.vec(num.sol,size),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
	same=integer(num.sol),
	indice=integer(1),
	info=integer(1),
	imodel=integer(1),
	nsol=integer(1),
	PACKAGE="wle")

  result$wstep <- z$wstep[1:z$imodel,]
  result$coefficients <- z$param[1:z$nsol,]
  result$scale <- sqrt(z$var[1:z$nsol])
  result$residuals <- z$resid[1:z$nsol,]
  result$tot.weights <- z$totweight[1:z$nsol]
  result$weights <- z$weight[1:z$nsol,]
  result$freq <- z$same[1:z$nsol]
  result$index <- z$indice
  result$info <- z$info
  result$call <- cl
  result$contrasts <- attr(xdata, "contrasts")
  result$xlevels <- xlev
  result$terms <- mt
  result$type <- type
  result$method <- method
  result$f.in <- f.in
  result$f.out <- f.out

  if (model)
    result$model <- mf
  if (ret.x)
    result$x <- xdata
  if (ret.y)
    result$y <- ydata

  dn <- colnames(xdata)

  if (is.null(nrow(result$coefficients)))
    names(result$coefficients) <- dn
  else
    dimnames(result$coefficients) <- list(NULL,dn)
  if (z$imodel<=1)
    names(result$wstep) <- c(dn," ")
  else
    dimnames(result$wstep) <- list(NULL,c(dn," "))

  class(result) <- "wle.stepwise"
  return(result)
}

#############################################################
#                                                           #
#	summary.wle.stepwise function                       #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 3, 2001                             #
#	Version: 0.4-1                                      #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

summary.wle.stepwise <- function (object, num.max=20, verbose=FALSE, ...) {

  if (is.null(object$terms))
    stop("invalid \'wle.stepwise\' object")
  if (num.max<1) {
    if (verbose)
      cat("summary.wle.stepwise: num.max can not less than 1, num.max=1 \n")
    num.max <- 1
  }

  ans <- list()
  wstep <- object$wstep
  if(is.null(nmodel <- nrow(wstep)))
    nmodel <- 1
  num.max <- min(nmodel,num.max)
  if (nmodel!=1)
    wstep <- wstep[(nmodel-num.max+1):nmodel,]

  ans$wstep <- wstep
  ans$num.max <- num.max
  ans$type <- object$type
  ans$f.in <- object$f.in
  ans$f.out <- object$f.out
  ans$call <- object$call

  class(ans) <- "summary.wle.stepwise"
  return(ans)
}

#############################################################
#                                                           #
#	print.wle.stepwise function                             #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: October, 27, 2003                                 #
#	Version: 0.4-1                                          #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.wle.stepwise <- function (x, digits = max(3, getOption("digits") - 3), num.max=max(1, nrow(x$wstep)), ...) {
  res <- summary.wle.stepwise(object=x, num.max=num.max, ...)
  print.summary.wle.stepwise(res, digits=digits, ...)
}

#############################################################
#                                                           #
#	print.summary.wle.stepwise function                 #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 3, 2001                             #
#	Version: 0.4                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

print.summary.wle.stepwise <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
  cat("\n",x$type," selection procedure\n")
  if (x$type=="Forward" | x$type=="Stepwise")
    cat("\nF.in: ",x$f.in)
  if (x$type=="Backward" | x$type=="Stepwise")
    cat("\nF.out: ",x$f.out)
  cat(" \n")
  cat("\nLast ",x$num.max," iterations:\n")

  if(x$num.max>1) {
    nvar <- ncol(x$wstep)-1
    x$wstep[,(nvar+1)] <- signif(x$wstep[,(nvar+1)],digits)
  } else {
    nvar <- length(x$wstep)-1
    x$wstep[(nvar+1)] <- signif(x$wstep[(nvar+1)],digits)
  }
  print(x$wstep)
  cat("\n")
  invisible(x)
}
