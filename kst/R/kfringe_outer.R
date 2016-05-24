### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### kfringe_outer.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(relations), library(sets)
###
### 2009-07-02: created
###

kfringe_outer <- function(x, state=NULL) {

   ### convert sets to relations
   relmat <- set_outer(x, set_is_proper_subset)
   relmat <- relation(domain=x, incidence=relmat)
   reldom <- as.list(relation_domain(relmat)[[1]])
   relmat <- relation_incidence(transitive_reduction(relmat))

   ### compute inner fringes
   ofringe <- list()
   for (i in 1:length(reldom)) {
      cs <- reldom[[i]]
      if (sum(relmat[i,])==0) {
         ofringe[[i]] <- set()
      } else {
         ns <- reldom[which(relmat[i,]==1)]
         ofr <- set()
         for (j in 1:length(ns)) {
            nsp <- ns[[j]]
            ofr <- c(ofr, set_symdiff(cs, nsp))
         }
         ofringe[[i]] <- ofr
      }
      names(ofringe)[[i]] <- colnames(relmat)[i]
   }

   ### select specified fringe
   if (!is.null(state)) {
      state <- LABEL(state)
      ofringe <- ofringe[[state]]
   }

   ### return results
   ofringe

}
