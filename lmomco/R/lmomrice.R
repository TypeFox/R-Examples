"lmomrice" <-
function(para, ...) {
    V   <- para$para[1]
    A   <- para$para[2]
    SNR <- V/A
    if(V == 0) {
      ray <- vec2par(c(0,A), type="ray")
      lmr <- lmomray(para=ray)
      lmr$source <- "lmomrice"
      return(lmr)
    }
    if(SNR > 52) {
      # Bessels are limited out, must compute remainder at
      # apparent numerical limits
      xbar <- A * SNR # just V no noise
      xvar <- A^2 # as SNR --> infinity: 2*A^2 + V^2 - A^2 * SNR^2
      nor  <- vec2par(c(xbar,sqrt(xvar)), type="nor")
      lmr  <- lmomnor(para=nor)
      lmr$source <- "lmomrice via lmomnor (super high SNR)"
      return(lmr)
    } else if(SNR > 24) {
      # pdfrice() can no longer integrate correctly
      # have numerically approached the Normal, but still can
      # compute the Bessels for Laguerre Polynomial
      L05  <- LaguerreHalf(-V^2/(2*A^2))
      xbar <- A * sqrt(pi/2) * L05
      xvar <- 2*A^2 + V^2 - A^2 * (pi/2) * L05^2
      nor  <- vec2par(c(xbar,sqrt(xvar)), type="nor")
      lmr  <- lmomnor(para=nor)
      lmr$source <- "lmomrice via lmomnor (high SNR)"
      return(lmr)
    }
    if(SNR < 0.08) {
      # theoLmoms.max.ostat() can no longer integrate correctly
      # have numerically approached the Rayleigh
      ray <- vec2par(c(0,A), type="ray")
      lmr <- lmomray(para=ray)
      lmr$source <- "lmomrice via lmomray (very low SNR)"
      return(lmr)
    }
    lmr <- theoLmoms.max.ostat(para=para,
                cdf=cdfrice, pdf=pdfrice,
                lower=0, upper=.Machine$double.max, ...)
    lmr$source <- "lmomrice"
    if(! are.lmom.valid(lmr)) {
      warning("The Rician parameters are producing invalid L-moments or L-moments outside ",
              "of implementation of Rice distribution in lmomco")
      warning("Printing the parameters")
      print(para)
      warning("Printing the Lmoments")
      print(lmr)
      return()
    }
    return(lmr)
}

