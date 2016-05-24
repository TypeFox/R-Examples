### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### reduction.kstructure.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets)
###
### 2008-04-24: created
###

reduction.kstructure <- function(x, operation=c("discrimination", "union", "intersection"),...) {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   if(operation == "discrimination") {
   ### compute discriminative reduction

      dom <- kdomain(x)
      notions <- knotions(x)
      kstates <- sapply(x, function(z) dom %in% z)
      rownames(kstates) <- dom
      redu <- set()
      for (i in x) {
         for (j in notions) {
            if (set_is_subset(j,i)) {
               i <- set_union(set_symdiff(i,j),set(paste(j, collapse="")))
            }
         }
         redu <- c(redu, set(as.set(i)))
      }

   } else
      redu <- NextMethod()

   class(redu) <- class(x)

   ### return results
   redu

}
