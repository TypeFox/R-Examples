`DnegLow` <- function(r,k,n,cc,chiquant,verbose=FALSE,cex=1,curcol="black",curtyp="solid") {
   # draw the lower cl for D<0, chisquare with cc
   ll <- r/(1+r)
   ul <- 1/(1+r)
 
   kla <- k^2
   klb <- -2*k - 1.5*k^2
   klc <- 1 + 4*k + 1.5*k^2 - chiquant/n
   kld <- -2 - 4*k + 2*chiquant/n
   kle <- 1 + 2*k - chiquant/n

   lcoef <- c(kla,klb,klc,kld,kle)
         
   lout <- polyroot(lcoef)
   Imlout <- round(abs(Im(lout)),digits=6)
   Relout <- round(Re(lout),digits=6)
   nlroot <- length(lout)

   lroots <- data.frame(1:nlroot,lout,Imlout,Relout,Imlout<1e-10,
   round(HWChisqccCurve(1-Relout,Relout,chiquant,n,curvetype="DnegLow"),digits=6),
   round(alowcurve(1-Relout,Relout,chiquant,n),digits=6))

   colnames(lroots) <- c("Root","Value","Im","Re","Im<1e-10","rAB","check")

   if(verbose) cat("Roots lower curve\n")

   ql <- CheckRoots(lroots,verbose=verbose)
         
   DnegLL <- ql
         
   if (!is.null(ql)) {
       q <- seq(ql,1-ql,by=0.005)
       p <- 1-q
       pt <-  2*(p-0.5)/sqrt(3)
       llc <- 2*p*q-2*cc*(1-p*q)/n-2*sqrt(cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
       points(pt,llc,type="l",col=curcol,cex=cex,lty=curtyp)
   }
   return(DnegLL=DnegLL)
 }

