`DposUp` <-
function(r,k,n,cc,chiquant,verbose=FALSE,cex=1,curcol="black",curtyp="solid") {
   # draw the upper cl for D>0, chisquare with cc
   ka <- -k^2
   kb <- 1.5*k^2
   kc <- chiquant/n - 1.5*k^2 + 2*k
   kd <- - 2*chiquant/n - 2*k
   ke <- chiquant/n - 1 + 2*k

   coef <- c(ka,kb,kc,kd,ke)

   out <- polyroot(coef)
   Imout <- round(abs(Im(out)),digits=6)
   Reout <- round(Re(out),digits=8)
   nroot <- length(out)

   uroots <- data.frame(1:nroot,out,Imout,Reout,Imout<1e-10,HWChisqccCurve(1-Reout,Reout,chiquant,n,
   curvetype="DposUp"),round(2*Reout,digits=8))

   colnames(uroots) <- c("Root","Value","Im","Re","Im<1e-10","rAB","check")

   if(verbose) cat("Roots upper curve\n")

   ql <- CheckRoots(uroots,verbose=verbose)

   DposUL <- ql
         
   if (!is.null(ql)) {
      q <- seq(ql,1-ql,by=0.005)
      p <- 1-q
      qt <-  2*(q-0.5)/sqrt(3)
      ul <- 2*p*q+2*cc*(1-p*q)/n+2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
      points(qt,ul,type="l",cex=cex,col=curcol,lty=curtyp)
   }
   return(DposUL=DposUL)
 }

