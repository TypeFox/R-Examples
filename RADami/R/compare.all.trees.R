compare.all.trees <-
function(treeset, ...) {
  ntrees = length(treeset)
  outmat = matrix(NA, ntrees, ntrees)
  for(i in 1:ntrees) {
   for(j in 1:i) {
    outmat[i, j] <- all.equal(treeset[[i]], treeset[[j]], ...)
	} # close j
  } # close i
  outmat
  }
