###############################################################################################
#####calculate the cdf using gaussian kernel with sd=bw for a vector of xx PHI((c0-xx)/bw)
#####################################################################

KernelSmoothing.cdf <-
function(xx,c0,bw)
  {
     
    new.xx <- (c0-xx)/bw
    ks.prob <- pnorm(new.xx,mean=0,sd=1,lower.tail=TRUE)##### #Note: normal CDF PHI((c0-xx)/bw) can be calculated by pnorm(c0,mean=xx,sd=bw) or by first transform yy=(c0-xx)/bw and then pnorm(yy,lower.tail=T)
    mean(ks.prob,na.rm=TRUE)##NOTE: use return(mean(ks.prob.na.rm=T)) takes more time
  }

