#############################################################
#                                                           #
#	mle.stepwise function                               #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: February,  10, 2010                           #
#	Version: 0.5                                        #
#                                                           #
#	Copyright (C) 2010 Claudio Agostinelli              #
#                                                           #
#############################################################

mle.stepwise <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, type="Forward", f.in=4.0, f.out=4.0, contransts=NULL, verbose=FALSE) {

  ntype <- switch(type,
    Forward = 1,
    Backward = 2,
    Stepwise = 3,
    stop("The type must be Forward, Backward or Stepwise"))

  ret.x <- x
  ret.y <- y
  result <- list()	
  mt <- terms(formula, data = data)
  mf <- cl <- match.call()
  mf$type <- mf$f.in <- mf$f.out <- NULL
  mf$max.iter <- mf$contrasts <- NULL
  mf$model <- mf$x <- mf$y <- NULL
  mf$verbose <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent()))
  xvars <- as.character(attr(mt, "variables"))[-1]
  inter <- attr(mt, "intercept")
  if((yvar <- attr(mt, "response")) > 0)
    xvars <- xvars[-yvar]
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
    stop("Number of observation must be at least equal to the number of predictors (including intercept) + 1")
  if (f.in<0 | f.out<0)
    stop("f.in and f.out can not be negative")

  nrep <- 2^nvar-1

  z <- .Fortran("step",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(0), 
	as.integer(size),
	as.integer(nvar),
	as.integer(nrep),
	as.integer(ntype),
	as.double(f.in),
	as.double(f.out),
	step=mat.or.vec(nrep,nvar+1),
	info=integer(1),
	imodel=integer(1),
	PACKAGE = "wle")

  result$step <- z$step[1:z$imodel,]
  result$info <- z$info
  result$call <- match.call()
  result$contrasts <- attr(xdata, "contrasts")
  result$xlevels <- xlev
  result$terms <- mt
  result$type <- type
  result$f.in <- f.in
  result$f.out <- f.out

  if (model)
    result$model <- mf
  if (ret.x)
    result$x <- xdata
  if (ret.y)
    result$y <- ydata

  dn <- colnames(xdata)

  if (z$imodel>0) {
    if (z$imodel==1) {
      names(result$step) <- c(dn," ")
    } else {
      dimnames(result$step) <- list(NULL,c(dn," "))
    }
  }

  class(result) <- "mle.stepwise"
  return(result)
}


summary.mle.stepwise <- function (object, num.max=20, verbose=FALSE, ...) {

  if (is.null(object$terms))
    stop("invalid \'mle.stepwise\' object")

  if (num.max<1) {
    if (verbose) cat("summary.mle.stepwise: num.max can not less than 1, num.max=1 \n")
    num.max <- 1
  }

  ans <- list()
  step <- object$step
  if(is.null(nmodel <- nrow(step)))
    nmodel <- 1
  num.max <- min(nmodel,num.max)
  if (nmodel!=1)
    step <- step[(nmodel-num.max+1):nmodel,]

  ans$step <- step
  ans$num.max <- num.max
  ans$type <- object$type
  ans$f.in <- object$f.in
  ans$f.out <- object$f.out
  ans$call <- object$call

  class(ans) <- "summary.mle.stepwise"
  return(ans)
}

print.mle.stepwise <- function (x, digits = max(3, getOption("digits") - 3), num.max=max(1,nrow(x$step)), ...) {
  res <- summary.mle.stepwise(object=x, num.max=num.max, ...)
  print.summary.mle.stepwise(res, digits=digits, ...)
}

print.summary.mle.stepwise <- function (x, digits = max(3, getOption("digits") - 3), ...) {
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
    nvar <- ncol(x$step)-1
    x$step[,(nvar+1)] <- signif(x$step[,(nvar+1)],digits)
  } else {
    nvar <- length(x$step)-1
    x$step[(nvar+1)] <- signif(x$step[(nvar+1)],digits)
  }
  print(x$step)
  cat("\n")
  invisible(x)
}
