### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### knotions.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets)
###
### 2008-04-24: created
###

knotions <- function(x) {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   ### calculate notions
   dom <- kdomain(x)
   kstates <- sapply(x, function(z) dom %in% z)
   rownames(kstates) <- dom
   notions <- set()
   for (i in seq_len(nrow(kstates))) {
      n <- which(apply(kstates,1,function(z)all(z==kstates[i,])))
      notions <- c(notions, set(as.set(names(n))))
   }

   ### return notions
   notions
}
