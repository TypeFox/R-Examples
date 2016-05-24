"plotlmrdia" <-
function(lmr,
      nopoints=FALSE,
      nolines=FALSE,
      nolimits=FALSE,
      noaep4=FALSE,
      nogev=FALSE,
      noglo=FALSE,
      nogpa=FALSE,
      nope3=FALSE,
      nogno=FALSE,
      nogov=FALSE,
      nocau=TRUE,
      noexp=FALSE,
      nonor=FALSE,
      nogum=FALSE,
      noray=FALSE,
      nosla=TRUE,
      nouni=FALSE,
         xlab="L-SKEW",
         ylab="L-KURTOSIS",
         autolegend=FALSE,xleg=NULL,yleg=NULL,
         ...) {

   entries <- vector(mode = "character")
   Elwd    <- vector(mode = "numeric")
   Epch    <- vector(mode = "numeric")
   Ecol    <- vector(mode = "numeric")
   Elty    <- vector(mode = "numeric")
   Ecex    <- vector(mode = "numeric")
   entryi  <- 0

   plot(lmr$limits, xlab = xlab, ylab = ylab, type = "n",
        ...)

   if(nolimits == FALSE) {
     lines(lmr$limits,lwd=2,col=8)
     entryi <- entryi + 1
     entries[entryi] <- "Theoretical limits"
     Elwd[entryi] <- 2
     Ecol[entryi] <- 8
     Epch[entryi] <- NA
     Elty[entryi] <- 1
     Ecex[entryi] <- 1
   }
   if(nolines == FALSE) {
     if(noaep4 == FALSE) {
        lines(lmr$aep4, col=2, lty=3)
        entryi <- entryi + 1
        entries[entryi] <- "AEP4"
        Elwd[entryi] <- 1
        Ecol[entryi] <- 2
        Epch[entryi] <- NA
        Elty[entryi] <- 4
        Ecex[entryi] <- 1
     }
     if(nogev == FALSE) {
        lines(lmr$gev, col=2,lty=2)
        entryi <- entryi + 1
        entries[entryi] <- "GEV"
        Elwd[entryi] <- 1
        Ecol[entryi] <- 2
        Epch[entryi] <- NA
        Elty[entryi] <- 2
        Ecex[entryi] <- 1
     }
     if(noglo == FALSE) {
        lines(lmr$glo, col=3)
        entryi <- entryi + 1
        entries[entryi] <- "GLO"
        Elwd[entryi] <- 1
        Ecol[entryi] <- 3
        Epch[entryi] <- NA
        Elty[entryi] <- 1
        Ecex[entryi] <- 1
     }
     if(nogno == FALSE) {
        lines(lmr$gno, col=4, lty=2)
        entryi <- entryi + 1
        entries[entryi] <- "GNO"
        Elwd[entryi] <- 1
        Ecol[entryi] <- 4
        Epch[entryi] <- NA
        Elty[entryi] <- 2
        Ecex[entryi] <- 1
     }
     if(nogov == FALSE) {
        lines(lmr$gov, col=6, lty=2)
        entryi <- entryi + 1
        entries[entryi] <- "GOV"
        Elwd[entryi] <- 1
        Ecol[entryi] <- 6
        Epch[entryi] <- NA
        Elty[entryi] <- 2
        Ecex[entryi] <- 1
     }
     if(nogpa == FALSE) {
        lines(lmr$gpa, col=4)
        entryi <- entryi + 1
        entries[entryi] <- "GPA"
        Elwd[entryi] <- 1
        Ecol[entryi] <- 4
        Epch[entryi] <- NA
        Elty[entryi] <- 1
        Ecex[entryi] <- 1
     }
     if(nope3 == FALSE) {
        lines(lmr$pe3, col=6)
        entryi <- entryi + 1
        entries[entryi] <- "PE3"
        Elwd[entryi] <- 1
        Ecol[entryi] <- 6
        Epch[entryi] <- NA
        Elty[entryi] <- 1
        Ecex[entryi] <- 1
     }
   }
   if(nopoints == FALSE) {
     if(nocau == FALSE) {
        points(lmr$cau,pch=13,col=3,cex=1.25)
        entryi <- entryi + 1
        entries[entryi] <- "CAU (limiting TL1)"
        Elwd[entryi] <- NA
        Ecol[entryi] <- 3
        Epch[entryi] <- 13
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.25
     }
     if(noexp == FALSE) {
        points(lmr$exp,pch=16,col=2,cex=1.5)
        entryi <- entryi + 1
        entries[entryi] <- "EXP"
        Elwd[entryi] <- NA
        Ecol[entryi] <- 2
        Epch[entryi] <- 16
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.5
     }
     if(nonor == FALSE) {
        points(lmr$nor,pch=15,col=2,cex=1.5)
        entryi <- entryi + 1
        entries[entryi] <- "NOR"
        Elwd[entryi] <- NA
        Ecol[entryi] <- 2
        Epch[entryi] <- 15
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.5
     }
     if(nogum == FALSE) {
        points(lmr$gum,pch=17,col=2,cex=1.5)
        entryi <- entryi + 1
        entries[entryi] <- "GUM"
        Elwd[entryi] <- NA
        Ecol[entryi] <- 2
        Epch[entryi] <- 17
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.5
     }
     if(noray == FALSE) {
        points(lmr$ray,pch=18,col=2,cex=1.5)
        entryi <- entryi + 1
        entries[entryi] <- "RAY"
        Elwd[entryi] <- NA
        Ecol[entryi] <- 2
        Epch[entryi] <- 18
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.5
     }
     if(nosla == FALSE) {
        points(lmr$sla,pch=10,cex=1.25,col=3)
        entryi <- entryi + 1
        entries[entryi] <- "SLA (TL1)"
        Elwd[entryi] <- NA
        Ecol[entryi] <- 3
        Epch[entryi] <- 10
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.25
     }
     if(nouni == FALSE) {
        points(lmr$uniform,pch=12,cex=1.25,col=2)
        entryi <- entryi + 1
        entries[entryi] <- "UNI"
        Elwd[entryi] <- NA
        Ecol[entryi] <- 2
        Epch[entryi] <- 12
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.25
     }
   }

   if(autolegend == TRUE & length(entries) > 0) {
     if(is.null(xleg)) warning("xleg is NULL, but needed")
     if(is.null(yleg)) warning("yleg is NULL, but needed")

     legend(xleg,yleg,entries,
            lwd=Elwd,
            col=Ecol,
            pch=Epch,
            lty=Elty,
            pt.cex=Ecex,
            xjust=0.5,bty="n",cex=0.9)
   }
}
