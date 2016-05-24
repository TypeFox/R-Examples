### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### kfringe_inner.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(relations), library(sets)
###
### 2009-06-03: created
###

kfringe_inner <- function(x, state=NULL) {

   ### convert sets to relations
   relmat <- set_outer(x, set_is_proper_subset)
   relmat <- relation(domain=x, incidence=relmat)
   reldom <- as.list(relation_domain(relmat)[[1]])
   relmat <- relation_incidence(transitive_reduction(relmat))

   ### compute inner fringes
   ifringe <- list()
   for (i in 1:length(reldom)) {
      cs <- reldom[[i]]
      if (sum(relmat[,i])==0) {
         ifringe[[i]] <- cs
      } else {
         ps <- reldom[which(relmat[,i]==1)]
         ifr <- set()
         for (j in 1:length(ps)) {
            psp <- ps[[j]]
            ifr <- c(ifr, set_symdiff(cs, psp))
         }
         ifringe[[i]] <- ifr
      }
      names(ifringe)[[i]] <- colnames(relmat)[i]
   }

   ### select specified fringe
   if (!is.null(state)) {
      state <- LABEL(state)
      ifringe <- ifringe[[state]]
   }

   ### return results
   ifringe

}
