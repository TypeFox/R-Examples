"ci.examp" <-
function(mean.sim=100, sd=10, n=25, reps=50,
                     conf.level=0.95, method="z",
                     lower.conf=(1-conf.level)/2,
                     upper.conf=1-(1-conf.level)/2 ) {

# This function demonstrates confidence intervals.  It will simulate
# data from a normal distribution and create multiple confidence
# intervals and plot all intervals of the mean along with a reference
# line indicating the mean.  lower.conf and upper.conf allow you to
# create unbalanced intervals.

  data <- matrix( rnorm( n*reps, mean.sim, sd), ncol=n)
  rmeans <- rowMeans(data)
  switch(method, Z=,z={
    lower <- qnorm( lower.conf, rmeans, sd/sqrt(n))
    upper <- qnorm( upper.conf, rmeans, sd/sqrt(n))
  },
         T=,t= {
           cv.l <- qt(lower.conf, n-1)
           cv.u <- qt(upper.conf, n-1)
           rsds <- sqrt( apply(data,1,var) )/sqrt(n)
           
           lower <- rmeans+cv.l*rsds
           upper <- rmeans+cv.u*rsds
         },
         BOTH=, Both=, both={
           lz <- qnorm( lower.conf, rmeans, sd/sqrt(n))
           uz <- qnorm( upper.conf, rmeans, sd/sqrt(n))
           
           cv.l <- qt(lower.conf, n-1)
           cv.u <- qt(upper.conf, n-1)
           rsds <- sqrt( apply(data,1,var) )/sqrt(n)
           
           lt <- rmeans+cv.l*rsds
           ut <- rmeans+cv.u*rsds
           
           lower <- c(rbind(lt,lz,mean.sim))
           upper <- c(rbind(ut,uz,mean.sim))
           
           reps <- reps*3
           rmeans <- rep(rmeans, each=3)
           rmeans[c(F,F,T)] <- NA
           
         },
         stop("method must be z, t, or both") )
  
  if( any( upper==Inf ) ) upper <- rep( 2*mean.sim-min(lower), reps )
  if( any( lower==-Inf ) ) lower <- rep( 2*mean.sim-max(upper), reps )
  
  xr <- range( upper, lower )
  
  plot(lower,seq(1,reps), type="n", xlim=xr, xlab="Confidence Interval",
       ylab="Index")
  
  abline( v= qnorm(c(1-upper.conf,1-lower.conf), mean.sim, sd/sqrt(n)), col=10)
  
  if( method=="both" || method=="Both" || method=="BOTH"){
    title( main="Confidence intervals based on both distributions",
          sub="Upper interval is Z in each pair")
  } else {
    title( main=paste("Confidence intervals based on",method,"distribution"))
  }
  
  colr <- ifelse( lower > mean.sim, 5, ifelse( upper < mean.sim, 6, 1) )
  
  abline(v=mean.sim)
  
  for( i in seq(1,reps) ){
    
    segments(lower[i], i, upper[i], i, col=colr[i])
    
  }
  
  points( rmeans, seq(along=rmeans), pch="|" )
  invisible(NULL)
}

