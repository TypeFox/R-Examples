tuneSnpRF <- function(x.autosome=NULL, x.xchrom=NULL, xchrom.names=NULL, x.covar=NULL, y, 
       	  	   mtryStart=floor(sqrt(sum(c(ncol(x.autosome),ncol(x.xchrom)/2,ncol(x.covar))))), 
		   ntreeTry=50, stepFactor=2, improve=0.05, trace=TRUE, plot=TRUE, doBest=FALSE, ...) {
  if (improve < 0) stop ("improve must be non-negative.")
  classRF <- is.factor(y)
  errorOld <- if (classRF) {
    snpRF(x.autosome=x.autosome, x.xchrom=x.xchrom, xchrom.names=xchrom.names, x.covar=x.covar, y=y, 
    	  mtry=mtryStart, ntree=ntreeTry, keep.forest=FALSE, ...)$err.rate[ntreeTry,1]
  } else {
     stop("Regression trees not implemented")
  }
  if (errorOld < 0) stop("Initial setting gave 0 error and no room for improvement.")
  if (trace) {
    cat("mtry =", mtryStart, " OOB error =",
        if (classRF) paste(100*round(errorOld, 4), "%", sep="") else
        errorOld, "\n")
  }

  oobError <- list()
  oobError[[1]] <- errorOld
  names(oobError)[1] <- mtryStart  
  
  for (direction in c("left", "right")) {
    if (trace) cat("Searching", direction, "...\n")
    Improve <- 1.1*improve
    mtryBest <- mtryStart
    mtryCur <- mtryStart
    while (Improve >= improve) {
      mtryOld <- mtryCur
      mtryCur <- if (direction == "left") {
        max(1, ceiling(mtryCur / stepFactor))
      } else {
        min(sum(dim(x.autosome)[2],dim(x.xchrom)[2]/2,dim(x.covar)[2],na.rm=T), floor(mtryCur * stepFactor))
      }
      if (mtryCur == mtryOld) break
      errorCur <- if (classRF) {
        snpRF(x.autosome=x.autosome, x.xchrom=x.xchrom, xchrom.names=xchrom.names, x.covar=x.covar, y=y, mtry=mtryCur,
	      ntree=ntreeTry, keep.forest=FALSE, ...)$err.rate[ntreeTry,"OOB"]
      } else {
        snpRF(x.autosome=x.autosome, x.xchrom=x.xchrom, xchrom.names=xchrom.names, x.covar=x.covar, y=y, mtry=mtryCur,
	      ntree=ntreeTry, keep.forest=FALSE, ...)$mse[ntreeTry]
      }
      if (trace) {
        cat("mtry =",mtryCur, "\tOOB error =",
            if (classRF) paste(100*round(errorCur, 4), "%", sep="") else
            errorCur, "\n")
      }
      oobError[[as.character(mtryCur)]] <- errorCur
      Improve <- 1 - errorCur/errorOld
      cat(Improve, improve, "\n")
      if (Improve > improve) {
        errorOld <- errorCur
        mtryBest <- mtryCur
      }
    }
  }
  mtry <- sort(as.numeric(names(oobError)))
  res <- unlist(oobError[as.character(mtry)])
  res <- cbind(mtry=mtry, OOBError=res)

  if (plot) {
    plot(res, xlab=expression(m[try]), ylab="OOB Error", type="o", log="x",
         xaxt="n")
    axis(1, at=res[,"mtry"])
  }

  if (doBest) 
    res <- snpRF(x.autosome=x.autosome, x.xchrom=x.xchrom, xchrom.names=xchrom.names, x.covar=x.covar, y=y,
    	         mtry=res[which.min(res[,2]), 1], ...)
  
  res
}
