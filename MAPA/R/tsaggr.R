tsaggr <- function(y,fout,fmean=c(FALSE,TRUE),outplot=c(FALSE,TRUE)){
# Non-overlapping temporal aggregation
# 
# Inputs:
#   y             Time series vector (can be ts object).
#   fout          Vector containing desirable aggregation levels. Must be positive
#                 and integer. If larger than length(y) then it is ignored.
#   fmean         If TRUE the aggregated is done using mean, otherwise sum is used. 
#   outplot       If TRUE a plot of the original series and the aggregated ones is 
#                 produced.
#
# Outputs:
#   out           List of temporally aggregated series. If y was a ts object, then 
#                 'out' has ts objects with appropriate frequencies. Any non-integer
#                 frequency is set equal to 1. Series are named ALx, where x is the 
#                 aggregation level.
#   all           An array containing all aggregated series in the original frequency.
#                 Series are named ALx, where x is the aggregation level.
#   idx           List of indices used to produce 'out' from 'all': 
#                 y.out[[i]] <- y.all[y.idx[[i]],i]
#                 Series are named ALx, where x is the aggregation level.
#
# Example:
#   out <- tsaggr(admissions,fout=2:12,fmean=TRUE,outplot=TRUE)
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>  

  # Defaults
  fmean <- fmean[1]
  outplot <- outplot[1]
  
  n <- length(y)
  
  # Get original frequency if y is ts
  if (class(y)=="ts"){
    f <- frequency(y)
  } else {
    f <- NULL
  }
  
  # Check fout
  fout <- fout[fout>0]
  fout <- fout[fout<=n]
  if (length(fout)==0){
    stop('fout must be positive integer(s), smaller than length(y).')
  }
  k <- length(fout)
  
  # Aggregation scale
  if (fmean == TRUE){
    scl <- fout
  } else {
    scl <- rep(1,k)
  }
  
  # Initialise array with aggregate series
  y.all <- array(NA,c(n,k))
  y.idx <- vector("list",k)
  y.out <- vector("list",k)
  slost <- vector("numeric",k)
  
  for (i in 1:k){
    # Find number of observations to drop from beginning
    slost[i] <- n %% fout[i]
    y.idx[[i]] <- seq((slost[i]+1),(n-fout[i]+1),fout[i])
    nk <- length(y.idx[[i]])
    # Aggregate
    for (j in 1:nk){
      y.all[y.idx[[i]][j]:(y.idx[[i]][j]+fout[i]-1),i] <- sum(y[y.idx[[i]][j]:(y.idx[[i]][j]+fout[i]-1)])/scl[i]
    }
    # Produce aggregate series
    y.out[[i]] <- y.all[y.idx[[i]],i]
    if (!is.null(f)){
      if (f %% fout[i] == 0){
        ftemp <- f/fout[i]
      } else {
        ftemp <- 1
      }
      y.out[[i]] <- ts(y.out[[i]],frequency=ftemp)
    }
  }
  
  # Name outputs
  names(y.out) <- paste0("AL",fout)
  names(y.idx) <- paste0("AL",fout)
  colnames(y.all) <- paste0("AL",fout)
  
  # Produce plot
  if (outplot==TRUE){
    ymin <- min(cbind(y,y.all),na.rm=TRUE)
    ymax <- max(cbind(y,y.all),na.rm=TRUE)
    ymin <- ymin - 0.1*(ymax-ymin)
    ymax <- ymax + 0.1*(ymax-ymin)
    cmp <- rainbow(k,start=0/6,end=4/6)
    plot(as.numeric(y),col="black",type="l",ylab="",xlab="Period",lwd=2,
         xaxs="i",yaxs="i",ylim=c(ymin,ymax))
    for (i in 1:k){
      lines(y.all[,i],col=cmp[i],type="l")
      points((y.idx[[i]]+fout[i]-1),y.all[y.idx[[i]],i],pch=20,col=cmp[i])
    }
    legend("topleft",c("Original",paste0("AL",fout[c(1,k)])),col=c("black",cmp[1],cmp[k]),lty=1,lwd=c(2,1,1),cex=0.7)
  }

  return(list(out=y.out,all=y.all,idx=y.idx))
  
}