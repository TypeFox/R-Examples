### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### as.relation.kstructure.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(relations), library(sets)
###
### 2008-05-29: created
###

as.relation.kstructure <- function(x,...) {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   ### compute relation
   dom <- kdomain(x)
   atoms <- katoms(x, items=dom)
   relmat <- mat.or.vec(length(dom), length(dom))
   colnames(relmat) <- dom
   rownames(relmat) <- dom
   for (i in 1:length(dom)) {
      items <- unique(unlist(atoms[[i]]))
      relmat[items,i] <- 1
   }
   rel <- as.relation(relmat)

   ### return results
   rel
}
