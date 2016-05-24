### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### plot.kstructure.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(relations), library(sets), library(Rgraphviz)
###
### 2008-04-17: created
###

plot.kstructure <- function(x, ...) {

   ### check if package Rgraphviz is available
   if (!require("Rgraphviz")) {
     stop(sprintf("Plotting requires package 'Rgraphviz'."))
   }

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   ### compute structure matrix
   relmat <- set_outer(x, set_is_subset)

   ### plot results
   plot(relation(incidence=relmat, domain=x), ...)

}
