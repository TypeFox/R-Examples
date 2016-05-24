### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### lpath.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets), library(relations)
###
### 2009-03-11: created
###

lpath <- function(x) {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   ### convert sets to relations
   relmat <- set_outer(x, set_is_proper_subset)
   relmat <- relation(domain=x, incidence=relmat)
   reldom <- as.list(relation_domain(relmat)[[1]])
   relmat <- relation_incidence(transitive_reduction(relmat)) == 1L

   ### compute learning paths
   lpath <- list()
   ts <- which(rowSums(relmat)==0)

   FUN <- function(j, path=NULL) {
       ns <- which(relmat[,j])
       if (length(ns)<1) {
           path <- as.tuple(reldom[c(path, ts[i])]) ## ts[i]??
           class(path) <- unique(c("lpath", class(path)))
           lpath <<- c(lpath, list(path))
       } else {
           for (k in 1:length(ns)) {
               FUN(ns[k], c(ns[k], path))
           }
       }
   }

   for (i in seq_along(ts)) {
      FUN(j = ts[i])
   }

   ### return learning paths
   lpath
}



