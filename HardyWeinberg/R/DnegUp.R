`DnegUp` <- function(r,k,n,cc,chiquant,verbose=FALSE,cex=1,curcol="black",curtyp="dotted")
  {
   # draw the upper cl for D<0, chisquare with cc
 
   kra <- k^2
   krb <- -1.5*k^2
   krc <- 1.5*k^2 + 2*k - chiquant/n
   krd <- 2*chiquant/n - 2*k
   kre <- 1 + 2*k - chiquant/n

   coefr <- c(kra,krb,krc,krd,kre)

   out <- polyroot(coefr)
   Imout <- round(abs(Im(out)),digits=6)
   Reout <- round(Re(out),digits=8)
   nroot <- length(out)

   lmidroots <- data.frame(1:nroot,out,Imout,Reout,Imout<1e-10,HWChisqccCurve(1-Reout,Reout,chiquant,n,curvetype="DnegUp"),round(2*Reout,digits=8))

   colnames(lmidroots) <- c("Root","Value","Im","Re","Im<1e-10","rAB","check")

   if(verbose) cat("D<0; Upper Curve\n")

   # paint the upper part of the curve
   ql <- CheckRoots(lmidroots,verbose=verbose,mini=FALSE)

   DnegUL <- ql
           
   if (!is.null(ql)) {
      q <- seq(ql,1-ql,by=0.005)
      p <- 1-q
      qt <-  2*(q-0.5)/sqrt(3)
      ul <- 2*p*q-2*cc*(1-p*q)/n+2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
      points(qt,ul,col=curcol,type="l",cex=cex,lty=curtyp)
   }

   # paint a second lower part of the curve.

   d1 <- (-3*cc^2 - 4*n*cc)/(2*n^2)
   d2 <- (2*cc-chiquant)/n + 1
     
   kla <- (cc/n)^2
   klb <- d1
   klc <- d2-d1
   kld <- -2*d2
   kle <- d2

   lcoef <- c(kla,klb,klc,kld,kle)
         
   lout <- polyroot(lcoef)
   Imlout <- round(abs(Im(lout)),digits=6)
   Relout <- Re(lout)
   nlroot <- length(lout)

   lroots <- data.frame(1:nlroot,lout,Imlout,Relout,Imlout<1e-10,round(HWChisqccCurve(1-Relout,Relout,chiquant,n,curvetype="DnegUp"),digits=10),0)

   colnames(lroots) <- c("Root","Value","Im","Re","Im<1e-10","rAB","check")
   

   if(verbose) {
      cat("Zeros Upper curve::\n")
      print(lroots)
   }
   
   ql <- CheckRoots(lroots,verbose=verbose)

   qlu <- CheckRoots(lmidroots,verbose=verbose,mini=TRUE)

   ql <- min(ql,qlu)
   qlu <- max(ql,qlu)
   if(verbose) print(ql)
         
   DnegUL <- ql
           
   if (!is.null(ql)) {
      q <- seq(ql,qlu,by=0.005)
      p <- 1-q
      qt <-  2*(q-0.5)/sqrt(3)
      ul <- 2*p*q-2*cc*(1-p*q)/n+2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
      points(qt,ul,col=curcol,type="l",cex=cex,lty=curtyp)

      q <- seq(1-qlu,1-ql,by=0.005)
      p <- 1-q
      qt <-  2*(q-0.5)/sqrt(3)
      ul <- 2*p*q-2*cc*(1-p*q)/n+2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
      points(qt,ul,col=curcol,type="l",cex=cex,lty=curtyp)
      
   }
   return(DnegUL=DnegUL)     
}

