### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### katoms.R
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###
### dependencies: library(sets)
###
### 2008-04-24: created
###

katoms <- function(x, items) {

   ### check x
   if (!inherits(x, "kstructure")) {
      stop(sprintf("%s must be of class %s.", dQuote("x"), dQuote("kstructure")))
   }

   ### check items
   if (!inherits(items, "set")) {
      stop(sprintf("%s must be of class %s.", dQuote("items"), dQuote("set")))
   }

   ### compute atoms
   x <- as.list(x)
   atoms <- list()
   items <- as.set(lapply(items, as.character))
   for (i in items) {
      states <- x[grep(i,x)]
      atom <- set()
      for (j in seq_along(states)) {
         subsets <- lapply(states[-j],set_is_subset, states[[j]])
         if (!any(unlist(subsets))) {
            atom <- c(atom, set(as.set(states[[j]])))
         }
      }
      atoms[[i]] <- atom
   }
   names(atoms) <- unlist(items)

   ### return atoms
   atoms
}
