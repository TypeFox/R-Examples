### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### kspace.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets)
###
### 2008-04-17: created
###

kspace <- function(x) {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   ### compute knowledge space
   dom <- kdomain(x)
   space <- c(x, set(dom), set(set()))
   class(space) <- class(x)
   space <- closure(space, operation="union")
   class(space) <- unique(c("kspace", class(space)))

   ### return space
   space
}
