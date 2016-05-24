`HWChisqLowerl` <-
function(r,verbose=FALSE,cex=1,curvecol="black") {
   # draw lower confidence limit for ordinary chisquare test
   if (r <=1) {         
      p <- seq(0,1,by=0.005)
      q <- 1 - p
      pt <- 2*(p-0.5)/sqrt(3)
      fplow <- 2*p*q - 2*p*q*r
      points(pt,fplow,type="l",col=curvecol,cex=cex)
    }
 }

