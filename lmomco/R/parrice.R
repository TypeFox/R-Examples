"parrice" <-
function(lmom, checklmom=TRUE) {
   para <- vector(mode="numeric", length=2)
   names(para) <- c("nu", "alpha")
   if(length(lmom$L1) == 0) { # convert to named L-moments
     lmom <- lmorph(lmom)     # nondestructive conversion!
   }
   if(checklmom & ! are.lmom.valid(lmom)) {
     warning("L-moments are invalid")
     return()
   }

   L1  <- lmom$L1
   LCV <- lmom$LCV
   if(is.null(L1) | is.null(lmom$LCV)) {
     warning("NULL L-moments")
     return()
   }
   IFAIL <- 0
   IFAILTEXT <- "Successful parameter estimation"
   minLCV <- min(.lmomcohash$RiceTable$LCV)
   maxLCV <- max(.lmomcohash$RiceTable$LCV)
   if(LCV >  maxLCV) {
      LCV <- maxLCV
      IFAIL <- 1
      IFAILTEXT <- "LCV too large for Rice (greater than Rayleigh), fitting Rayleigh instead"
   } else if(LCV < minLCV) {
      warning("LCV too small (<",minLCV,") for Rice as implemented")
      return()
   }
   SNR  <- approx(.lmomcohash$RiceTable$LCV, .lmomcohash$RiceTable$SNR, xout=LCV)$y
   G    <- approx(.lmomcohash$RiceTable$LCV, .lmomcohash$RiceTable$G,   xout=LCV, rule=1:2)$y
   A    <- L1/G
   V    <- A*SNR
   para[1] <- V
   para[2] <- A
   return(list(type='rice', para=para, source="parrice",
               ifail=IFAIL,
               ifailtext=IFAILTEXT))
}
