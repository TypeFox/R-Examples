`HWUppercl` <- function(r,verbose=FALSE,cex=1,curvecol="black") {
   # draw upper confidence limit for ordinary chisquare test
   ll <- r/(1+r)
   ul <- 1/(1+r)
   
   if(verbose) cat("Roots HW curve",ll,ul,"\n")
         
   if ((ll<=0.5) & (ul>=0.5)) {
      p <- seq(ll,ul,by=0.005)
      q <- 1 - p
      pt <- 2*(p-0.5)/sqrt(3)
      fpup <- 2*p*q + 2*p*q*r
      points(pt,fpup,type="l",col=curvecol,cex=cex)
   }
 }

