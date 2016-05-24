`DposLow` <-
function(r,k,n,cc,chiquant,verbose=FALSE,cex=1,curcol="black",curtyp="dotted") {
   # draw the lower cl for D>0, chisquare with cc
   ka <- k^2
   kb <- -1.5*k^2
   kc <- 1.5*k^2 - 2*k - chiquant/n
   kd <- 2*chiquant/n + 2*k
   ke <- 1 - 2*k - chiquant/n

   coef <- c(ka,kb,kc,kd,ke)

   out <- polyroot(coef)
   Imout <- round(abs(Im(out)),digits=6)
   Reout <- round(Re(out),digits=8)
   nroot <- length(out)

   umidroots <- data.frame(1:nroot,out,Imout,Reout,Imout<1e-10,HWChisqccCurve(1-Reout,Reout,chiquant,n,curvetype="DposLow"),round(2*Reout,digits=8))



   colnames(umidroots) <- c("Root","Value","Im","Re","Im<1e-10","rAB","check")

   if(verbose) cat("upper (D>0) mid curve\n")

   ql <- CheckRoots(umidroots,verbose=verbose)

   DposLL <- ql
         
   if (!is.null(ql)) {
      q <- seq(ql,1-ql,by=0.005)
      p <- 1-q
      qt <-  2*(q-0.5)/sqrt(3)
      ul <- 2*p*q+2*cc*(1-p*q)/n-2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
      points(qt,ul,col=curcol,type="l",cex=cex,lty=curtyp)
   }
   return(DposLL=DposLL)
 }

