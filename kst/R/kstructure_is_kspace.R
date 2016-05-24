### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### kstructure_is_kspace.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets)
###
### 2008-04-24: created
###

kstructure_is_kspace <- function(x) {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   ### check for space
   y <- inherits(x, "kspace") || length(x)==length(kspace(x))

   ### return results
   y
}
