assign("est.variograms",
function(point.obj,pair.obj,a1,a2,trim) {

#  est.variogram takes a "point" object, point.obj, and a "pair" object, pair.obj,
#  calculates empirical variogram estimates.


#  The result is an object of type "variogram" with 4 components: $lags, 
#  $classic, $robust, and $n.
#
#  $lags - lag category
#  $bins - distance bins for plotting
#  $classic - the classic variogram estimator
#  $robust - the robust variogram estimator
#  $med - the median estimator
#  $trimmed.mean - the median trimmed estimator
#  $n - the number of pairs in the lag
#  library(sgeostat, pos=which(search()=="package:gstat")+1)		

  if (!inherits(point.obj,"point")) stop('Point.obj must be of class, "point".\n')

  if (!inherits(pair.obj,"pair")) stop('Pair.obj must be of class, "pair".\n')

  if(missing(a1)) stop('Must enter at least one attribute.\n')
  if(missing(a2)) a2 <- a1

  a1 <- point.obj[[match(a1,names(point.obj))]]
  a2 <- point.obj[[match(a2,names(point.obj))]]

# Allocate some space...
  lags         <- sort(unique(pair.obj$lags))
  classic      <- rep(0,length(lags))
  robust       <- rep(0,length(lags))
  med          <- rep(0,length(lags))
  trimmed.mean <- rep(0,length(lags))  
  n            <- rep(0,length(lags))

  diff <- a1[pair.obj$from]-a2[pair.obj$to]
  bo   <- split(diff,pair.obj$lags)

# this fails sometimes:
#  tmp<-unique(pair.obj$lags[-which.na(pair.obj$lags)])
# so do this:
  if (any(is.na(pair.obj$lags))) 
     tmp <- unique(pair.obj$lags[-which.na(pair.obj$lags)])
  else 
     tmp <- unique(pair.obj$lags)

  for (i in c(1:length(tmp))) {
#  for (i in unique(pair.obj$lags)) {
    n[i] <- length(bo[[i]][!is.na(bo[[i]])])

#   classic, see Matheron
    classic[i] <- sum((bo[[i]])^2,na.rm=TRUE) / n[i]

#   robust , med & trimmed mean, see Cressie, 1990
    robust[i]       <- (sum(abs(bo[[i]])^.5,na.rm=TRUE) / n[i] )^4 / (0.457 + (0.494/n[i]))
    med[i]          <- (median(abs(bo[[i]])^.5,na.rm=TRUE))^4 / (0.457 + (0.494/n[i]))
    trimmed.mean[i] <- (mean(abs(bo[[i]])^.5, trim=trim, na.rm=TRUE))^4 / (0.457 + (0.494/n[i]))
#    med[i] <- median(abs(bo[[i]]),na.rm=TRUE)
  }    
  o.variogram <- data.frame(lags,bins=c(pair.obj$bins,recursive=TRUE),classic,
                      robust,med,trimmed.mean,n=n)
  class(o.variogram) <- c("variogram","data.frame")
  return(o.variogram)

})

