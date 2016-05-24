#############################################################
#                                                           #
#	wle.cv.one function                                 #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: October, 13, 2012                             #
#	Version: 0.1-3                                      #
#                                                           #
#	Copyright (C) 2012 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.cv.one <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, monte.carlo=500, split, boot=30, group, num.sol=1, raf="HD", smooth=0.031, tol=10^(-6), equal=10^(-3), max.iter=500, min.weight=0.5, contrasts=NULL, num.max, nmodel, verbose=FALSE, file, append=FALSE, ...) {

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
if (!missing(file)) {
    if (!append) {
        if (file.exists(file)) file.remove(file) 
    }
}

    info <- wcv <- vector(length=0)
    ret.x <- x
    ret.y <- y
    result <- list()	
    mt <- terms(formula, data = data)
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$monte.carlo <- mf$split <- NULL
    mf$boot <- mf$group <- mf$smooth <- NULL
    mf$tol <- mf$equal <- mf$num.sol <- NULL
    mf$min.weight <- mf$max.iter <- mf$raf <- NULL
    mf$contrasts <- NULL
    mf$model <- mf$x <- mf$y <- NULL
    mf$verbose <- NULL
    mf$num.max <- mf$nmodel <- mf$file <- mf$... <- NULL
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
if (missing(nmodel)) {
    nmodel <- 1:nrep
}

if (missing(num.max)) {
    num.max <- length(nmodel)
}

if (size<nvar) {
stop("Number of observations must be at least equal to the number of predictors (including intercept)")
}

if (group<nvar) {
    group <- max(round(size/4),(nvar+1))
    group <- min(size, group)
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

    z <- .Fortran("wleregfix",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(0), 
	as.integer(size),
	as.integer(nvar),
	as.integer(nvar),
	as.integer(boot),
	as.integer(group),
	as.integer(num.sol),
	as.integer(raf),
	as.double(smooth),
	as.double(tol),
	as.double(equal),
	as.integer(max.iter),
	param=mat.or.vec(num.sol,nvar),
	var=double(num.sol),
	resid=mat.or.vec(num.sol,size),
	totweight=double(num.sol),
	weight=mat.or.vec(num.sol,size),
        density=mat.or.vec(num.sol,size),
        model=mat.or.vec(num.sol,size),
        delta=mat.or.vec(num.sol,size),
	same=integer(num.sol),
	nsol=integer(1),
	nconv=integer(1),
	PACKAGE = "wle")

   if (z$nconv==boot) stop("wle.cv.one: No solutions in the full model")
##          info=1

   indice <- 0
   wvaria <- z$var
   wtotpesi <- z$totweight
   dvar <- wvaria[1]+1
   for (i in 1:z$nsol) {
        if (dvar>wvaria[i] & min.weight<wtotpesi[i]) {
            dvar <- wvaria[i]
            indice <- i
        }
   }

   if (indice==0) stop("wle.cv.one: No solutions in the full model satisfied the criterion")
###       info=2
  dpesi <- z$weight[indice,]
  z$index <- indice

  if (verbose) {
      cat("Found ", z$nsol, "solution/s \n")
  }
  
for (imodel in nmodel) {

  if (verbose) {
      cat("Running model with number ", imodel, " on ", nrep, "\n")
  }
    
  zcv <- .Fortran("wlecvone",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(0), 
	as.integer(size),
	as.integer(nvar),
	as.integer(monte.carlo),
        as.integer(imodel),
	as.integer(split),
        as.double(dpesi),
	wcv=double(nvar+1),
	info=integer(1),
	PACKAGE="wle")
   info <- c(info, zcv$info)
   if (missing(file)) {
       if (NROW(wcv)<num.max) {
           wcv <- rbind(wcv, zcv$wcv)
       } else {
           pos <- which.max(wcv[,nvar+1])
           if (wcv[pos,nvar+1] > zcv$wcv[nvar+1]) wcv[pos,] <- zcv$wcv
       }
   } else { 
       write.table(x=matrix(c(imodel, zcv$wcv), nrow=1), file=file, append=TRUE, ...)
   }
  
}

if (missing(file)) {

    delnull <- z$same==0

    result$wcv <- wcv
    result$coefficients <- z$param[!delnull,]
    result$scale <- sqrt(z$var[!delnull])
    result$residuals <- z$resid[!delnull]
    result$weights <- z$weight[!delnull,]
    result$tot.weights <- z$totweight[!delnull]
    result$freq <- z$same[!delnull]
    result$call <- cl
    result$info <- info
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
} else {
    result <- cat("Results saved to file ", file, "\n")
}

    
return(result)
}

