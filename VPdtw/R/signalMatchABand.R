signalMatchABand <- function(reference,query, lambda=rep(0.0,length(reference)), maxshift=50) {

### Change on August 2nd 2007.  Created vectors of penalties above

### Performs nonlinear monotone transformation of the time
### axis to align to signals
### s and r are assumed to be two signals of the same length

### Uses implementation of the DTW algorithm to find the shortest path
### through a square matrix from top left to bottom right
### (path length is sum of the elements)
###
### s and r are the signals
### lambda is a non-diagonal step penalty
###
### s is assumed to be reference and r is aligned to it in this version.

  nr <- length(reference)
  nq <- length(query)
  ## lambda <- rep(0,length(reference));maxshift=20
  if(nq + maxshift < nr) stop("length(query) + maxshift should be greater than or equal to length(reference)\n")
  ## if(length(query) != n) stop("Signals need to same length in this implementation")

  pp <- .C("signalMatchWrapABand",
           reference=as.double(reference),
           query=as.double(query),
           nr=as.integer(nr),
           nq=as.integer(nq),
           path=integer(nr),
           lambda=as.double(lambda),
           maxs=as.integer(maxshift),
           dup=FALSE,
           PACKAGE="VPdtw")

  path <- pp$path
  path[path==0] <- NA
  minp <- min(path,na.rm=TRUE)
  maxp <- max(path,na.rm=TRUE)

  xIndices <- path; xVals <- 1:length(path)##pp$reference)
  if(minp>1) {
    xIndices <- c(1:(minp-1),xIndices)
    xVals <- c(seq(to=0,len=minp-1,by=1),xVals)
  }
  if(maxp<length(query)) {
    xIndices <- c(xIndices,(maxp+1):length(query))
    xVals <- c(xVals,seq(from=max(xVals)+1,by=1,len=length(query)-(maxp)))
  }

  if(FALSE) {
    plot(reference,type="l",lwd=2,xlim=c(1-maxshift,nr+maxshift))
    lines(which(!is.na(path)),query[na.omit(path)],col=2)
    lines(xVals,query[xIndices],col=3,lty=2)
  }

  ## Should we warn when we get close to maxshift?
  shift <- xVals - xIndices
  if(max(abs(shift),na.rm=TRUE)>(3*maxshift/4)) cat("Warning: Observed shift more than three quarters of maxshift\n")

  ## Come up with a nice summary
  output <- matrix(NA,length(xVals),4)
  colnames(output) <- c("xVals","reference","warped query","shift")
  output[,"xVals"] <- xVals
  str <- which(xVals==1)
  end <- which(xVals==nr)
  output[seq(str,end,by=1),"reference"] <- reference

  output[,"warped query"] <- query[xIndices]
  output[,"shift"] <- shift


  ##cat(range(pp$path - 1:length(pp$path)),"\n")
  zz <- output
  return(invisible(zz))
}

