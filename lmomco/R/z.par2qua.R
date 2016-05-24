"z.par2qua" <-
function(f,p,para,z=0,...) {
   if(is.null(p)) {
      warning("p is NULL, this function will not assume p=0, returning NULL")
      return(NULL)
   }
   if(length(p) != 1) {
      warning("only the first element of scalar argument p will be used")
      p <- p[1]
   }
   if(! check.fs(p)) return()

   if(length(z) != 1) {
      warning("only the first element of scalar argument z will be used")
      z <- z[1]
   }
   # assume f and para are valid and qlmomco() will check that anyway
   z.of.fit <- par2qua(0,para,...)
   if(z.of.fit <= z) {
      warning("evidently inconsistent z argument relative to that of the ",
              "fitted distribution, returning NULL")
   }

   x <- vector(mode="numeric", length(f))
   nF <- (f - p)/(1 - p)
   x[f >  p] <- par2qua(nF[f > p], para, ...)
   x[f <= p] <- z
   names(x) <- NULL
   return(x)
}
