dsmooth.variogram <-
function(vout)
{
  ## vout : output from dvariogram
  ## value : interpolated variogram at the first lag
  chks <- with(vout$gamma.s0, mean(n==0))
  ds <- with(vout$gamma.s0,mean(diff(x))) 
  chkt <- with(vout$gamma.0t, mean(n==0))
  dt <- with(vout$gamma.0t,mean(diff(x))) 
  dd <- c(ds,dt)
  if( chks >0.5){
	 ii <- with(vout$gamma.s0, which(n==0))
	 coords <- with(vout, cbind(gamma.s0$x[ii],tdist[1]))
     ghat <- working.smoothvariogram(coords,vout$gamma.long,dd)
	 vout$gamma.s0[ii,c(2,3)] <- cbind(1,ghat) ## assume interpolated n =1
	 ## vout$gamma.long <- rbind(vout$gamma.long, cbind(coords,ghat,1))
     tmp <- cbind(coords,ghat,1)
	 colnames(tmp) <- colnames(vout$gamma.long)
	 vout$gamma.long <- rbind(vout$gamma.long,tmp)
	 vout$gamma[1,ii] <- ghat
  }
  if ( chkt >0.5){
	 ii <- with(vout$gamma.0t, which(n==0))
	 coords <- with(vout, cbind(tdist[1],gamma.0t$x[ii]))
	 ghat <- working.smoothvariogram(coords,vout$gamma.long,dd)
	 vout$gamma.0t[ii,c(2,3)] <- cbind(1,ghat)
	 tmp <- cbind(coords,ghat,1)
	 colnames(tmp) <- colnames(vout$gamma.long)
	 vout$gamma.long <- rbind(vout$gamma.long,tmp)
	 vout$gamma[ii,1] <- ghat
  }
  vout
}
