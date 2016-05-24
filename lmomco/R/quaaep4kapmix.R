"quaaep4kapmix" <-
function(f, lmom, checklmom=TRUE) {
    if(! check.fs(f)) return()
    docheck <- TRUE
    if(checklmom) {
      if(! are.lmom.valid(lmom)) {
        warning("L-moments are invalid.")
        return()
      }
      docheck <- FALSE
    }
    # nondestructive conversion!
    if(length(lmom$L1) != 0) lmom <- lmorph(lmom)
    T3 <- lmom$ratios[3]; T4 <- lmom$ratios[4]
    A <- abs(T3)
    co <- c(0.7755464,  -3.3354852,  14.1955782, -29.9090294,
            37.2141451, -24.7411869,   6.7997646)
    T4.aep4low <-   co[1]*A   + co[2]*A^2 + co[3]*A^3 +
                    co[4]*A^4 + co[5]*A^5 + co[6]*A^6 + co[7]*A^7
    T4.kapup <- (5*A^2+1)/6
    if(T4.aep4low < T4 && T4 < T4.kapup) {
       W <- (T4.kapup - T4) / (T4.kapup - T4.aep4low)
       aep4 <- paraep4(lmom, method="A",
                       kapapproved=FALSE, checklmom=docheck)
       if(aep4$ifail > 0) {
          warning("It seems a mixture should be ok, but failure by the aep4 algorithms ",
                  "reporting invalid parameters, reverting to kappa only if possible")
          kap  <- parkap(lmom, checklmom=docheck)
          if(kap$ifail > 0) {
             return(rep(NA, length(f)))
          } else {
             return(par2qua(f, kap))
          }
       } else {
          kap  <- parkap(lmom, checklmom=docheck)
          return(par2qua2(f, aep4, kap, weight=W))
       }
    } else {
      aep4.or.kap <- paraep4(lmom, method="A",
                             kapapproved=TRUE, checklmom=docheck)
      return(par2qua(f, aep4.or.kap))
    }
}

