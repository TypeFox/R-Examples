eqwma <-
function(x, length=5, lag=1, start=1, p=1,
  log=FALSE, abs=FALSE, as.vector=TRUE)
{
#zoo related:
if(is.zoo(x)){
  xindex <- index(x)
  x <- coredata(as.vector(x))
  zoochk <- TRUE
}else{
  x <- coredata(as.vector(x))
  zoochk <- FALSE
}

#compute eqwma:
xn <- length(x)
x <- na.trim(x, sides="left")
xnadj <- length(x)
if(xn > xnadj){nachk <- TRUE}else{nachk <- FALSE}
if(abs){xabs <- abs(x)}else{xabs <- x}
if(p!=1){xabsp <- xabs^p}else{xabsp <- xabs}

out <- NULL
for(j in 1:length(length)){
  movavg <- rollmean(xabsp, length[j], na.pad=TRUE, align="r")
  if(start <= I(length[j]-1)){
    for(i in start:I(length[j]-1)){
      movavg[i] <- sum(xabsp[1:i])/i
    }
  }
  movavg <- c(rep(NA,lag), as.vector(movavg[1:I(length(movavg)-lag)]))
  out <- cbind(out, movavg)
}

#apply log?:
if(log){
  out <- log(out)
  colnames(out) <- paste("logEqWMA(", length, ")", sep="")
}else{
  colnames(out) <- paste("EqWMA(", length, ")", sep="")
}
if(nachk) out <- rbind( matrix(NA, I(xn-xnadj), dim(out)[2]), out)
if(as.vector){
  if( dim(out)[2]==1 ){ out <- as.vector(out) }
}
if(zoochk){ out <- zoo(out, order.by=xindex) }

return(out)
}
