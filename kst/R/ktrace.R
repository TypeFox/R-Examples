### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### ktrace.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets)
###
### 2008-04-24: created
###

ktrace <- function(x, items) {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   ### check items
   if (!inherits(items, "set")) {
      stop(sprintf("%s must be of class %s.", dQuote("items"), dQuote("set")))
   }

   ### compute trace
   items <- as.set(lapply(items, as.character))
   trace <- as.set(lapply(x, set_intersection, items))
   class(trace) <- class(x)

   ### return results
   trace
}
