#############################################################
#                                                           #
#	mle.cv.one function                                 #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: October, 13, 2012                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 2012 Claudio Agostinelli              #
#                                                           #
#############################################################

mle.cv.one <- function(formula, data=list(), model=TRUE, x=FALSE, y=FALSE, monte.carlo=500, split, contrasts=NULL, num.max, nmodel, verbose=FALSE, file, append=FALSE, ...) {

if (missing(split)) {
    split <- 0
}
if (!missing(file)) {
    if (!append) {
        if (file.exists(file)) file.remove(file) 
    }
}

    info <- cv <- vector(length=0)
    ret.x <- x
    ret.y <- y
    result <- list()	
    mt <- terms(formula, data = data)
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$monte.carlo <- mf$split <- mf$contrasts <- NULL
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

if (size<nvar+1) {stop("Number of observation must be at least equal to the number of predictors (including intercept) + 1")}

if (split<nvar+2 | split>(size-2)) {
    split <- max(round(size^(3/4)),nvar+2)
    if (verbose) cat("mle.cv: dimension of the split subsample set to default value = ",split,"\n")
}
maxcarlo <- sum(log(1:size))-(sum(log(1:split))+sum(log(1:(size-split))))
if (monte.carlo<1 | log(monte.carlo) > maxcarlo){
    stop("MonteCarlo replication not in the range")
}

for (imodel in nmodel) {

  if (verbose) {
      cat("Running model with number ", imodel, " on ", nrep, "\n")
  }
  
  z <- .Fortran("mlecvone",
	as.double(ydata),
	as.matrix(xdata),
	as.integer(0), 
	as.integer(size),
	as.integer(nvar),
	as.integer(monte.carlo),
        as.integer(imodel),
	as.integer(split),
        cv=double(nvar+1),
	info=integer(1),
	PACKAGE="wle")
   info <- c(info, z$info)

   if (missing(file)) { 
       if (NROW(cv)<num.max) {
           cv <- rbind(cv, z$cv)
       } else {
           pos <- which.max(cv[,nvar+1])
           if (cv[pos,nvar+1] > z$cv[nvar+1]) cv[pos,] <- z$cv
       }
   } else {
       write.table(x=matrix(c(imodel, z$cv), nrow=1), file=file, append=TRUE, ...)
   }
}

if (missing(file)) {

    result$cv <- cv
    result$call <- cl
    result$info <- info
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
} else {
    result <- cat("Results saved to file ", file, "\n")  
}
 
return(result)
}

