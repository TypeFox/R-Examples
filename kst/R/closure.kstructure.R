### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### closure.kstructure.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets)
###
### 2008-04-17: created
###

closure.kstructure <- function(x, operation=c("union", "intersection"), ...) {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   clos <- if(operation == "union") {
   ### compute closure under union
      dom <- kdomain(x)
      relmat <- t(sapply(x, function(z) dom %in% z))
      relmat <- binary_closure(relmat, operation)
      relmat <- relmat[order(rowSums(relmat)),]
      colnames(relmat) <- dom
      as.set(apply(relmat,1,function(z)as.set(names(which(z)))))      

   } else
      NextMethod()

   class(clos) <- class(x)
      
   ### return closure
   clos

}
