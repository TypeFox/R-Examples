"pdfrice" <-
function(x,para) {
   if(! are.parrice.valid(para)) return()
   V   <- para$para[1]
   A   <- para$para[2]

   if(V == 0) {
      ray <- vec2par(c(0,A), type="ray")
      return(pdfray(x,para=ray))
   }
   SNR <- V/A
   if(SNR > 52) {
      xbar <- A * SNR
      xvar <- A^2; # as SNR --> infinity: 2*A^2 + V^2 - A^2 * SNR^2
      nor  <- vec2par(c(xbar,sqrt(xvar)), type="nor")
      return(pdfnor(x,para=nor))
   } else if(SNR > 24) {
      L05  <- LaguerreHalf(-V^2/(2*A^2))
      xbar <- A * sqrt(pi/2) * L05
      xvar <- 2*A^2 + V^2 - A^2 * (pi/2) * L05^2
      nor  <- vec2par(c(xbar,sqrt(xvar)), type="nor")
      return(pdfnor(x,para=nor))
   }
   B  <- V/A^2

   TMP  <- x/(A^2) * exp( -(x^2 + V^2) / (2*A^2) )
   toIo <- x*B
   ops  <- options(warn=-1)
   Bo   <- besselI(toIo, nu=0, expon.scaled = TRUE)
   f    <- exp(log(TMP) + (log(Bo) + toIo))
   options(ops)

   names(f) <- NULL
   f[! is.finite(f)] <- NA
   f[is.na(f)] <- 0 # decision Dec. 2015
   return(f)
}

