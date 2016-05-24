"cdfrice" <- function(x, para=NULL) {
   if(! are.parrice.valid(para)) return()
   V   <- para$para[1]
   A   <- para$para[2]

   if(V == 0) {
     ray <- vec2par(c(0,A), type="ray")
     return(cdfray(x,para=ray))
   }
   SNR <- V/A
   if(SNR > 52) {
      xbar <- A * SNR
      xvar <- A^2; # as SNR --> infinity: 2*A^2 + V^2 - A^2 * SNR^2
      nor  <- vec2par(c(xbar,sqrt(xvar)), type="nor")
      return(cdfnor(x,para=nor))
   } else if(SNR > 24) {
      L05  <- LaguerreHalf(-V^2/(2*A^2))
      xbar <- A * sqrt(pi/2) * L05
      xvar <- 2*A^2 + V^2 - A^2 * (pi/2) * L05^2
      nor  <- vec2par(c(xbar,sqrt(xvar)), type="nor")
      return(cdfnor(x,para=nor))
   }
   f <- sapply(1:length(x), function(i) {
                       if(x[i] < 0)    return(0)
                       if(x[i] == Inf) return(1)
                       return(integrate(pdfrice, 0, x[i], para=para)$value)  })
   names(f) <- NULL
   return(f)
}

