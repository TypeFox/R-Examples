### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### kbase.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets)
###
### 2008-04-24: created
###

kbase <- function(x) {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }
   if (!kstructure_is_kspace(x)) {
      stop("'x' must be a knowledge space.")
   }

   ### compute base
   atoms <- katoms(x, items=as.set(unique(unlist(x))))
   base <- do.call(c,atoms)
   names(base) <- NULL
   class(base) <- unique(c("kstructure", class(base)))

   ### return base
   base
}
