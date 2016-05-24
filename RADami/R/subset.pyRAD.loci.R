subset.pyRAD.loci <-
function(x, loci, taxa, format = 'DNAStringSet', reportInterval = 500, nucVarType = c("verystrict", "strict", "relaxed"), ...) {
## Arguments:
##  x = pyRAD.loci object
##  loci = loci to use (by name)
##  taxa = taxa to use (by name)
##  format = only DNAStringSet export supported now
##  reportInterval = interval to use for reporting subsetting progress
##  nucVarType = either strict for requiring variability to be due to only unambiguous nucleotides (strict, verystrict) or allowing ambiguities to encode variability (relaxed), and insisting on parsimony informativeness (verystrict) or just variability (strict, relaxed) at at least one site
  excludedNucs <- switch(nucVarType[1], verystrict = 5:17, strict = 5:17, relaxed = 15:17)
  nucThresh <- switch(nucVarType[1], verystrict = 2, strict = 1, relaxed = 1)
  out <- list(
    DNA = structure(vector('list', length(loci)), names = loci),
    variable = structure(logical(length(loci)), names = loci),
	ntaxa = structure(integer(length(loci)), names = loci)
	)
  inds.vector <- x$tips %in% taxa
  counter = 0
  start.time <- Sys.time()
  for(i in loci) {
	counter <- counter + 1
	if(counter / reportInterval - counter %/% reportInterval == 0) {
  	   message(paste('... subsetting', counter, 'of', length(loci), 
 	   '-- Estimated time remaining =', round(((Sys.time() - start.time) / counter) * (length(loci) - counter), 1), attr(Sys.time() - start.time, 'units')
  	   ))
	}
    seq.index <- x$locus.index == i & inds.vector
	out$DNA[[i]] <- DNAStringSet(x$seqs[seq.index])
	names(out$DNA[[i]]) <- x$tips[seq.index]
	out$variable[i] <- any(apply(consensusMatrix(out$DNA[[i]])[-c(excludedNucs), ], 2, function(x) sum(x >= nucThresh) > 1))
	out$ntaxa[i] <- sum(seq.index)
	}
  return(out)
  }
