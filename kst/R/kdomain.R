### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### kdomain.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets)
###
### 2008-06-18: created
###

kdomain <- function(x) {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   ### compute domain
   domain <- as.set(unique(unlist(as.list(x))))

   ### return domain
   domain
}
