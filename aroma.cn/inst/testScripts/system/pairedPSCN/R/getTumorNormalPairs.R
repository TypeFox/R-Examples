getTumorNormalPairs <- function(ds, pairs, ..., verbose=TRUE) {
  ## Split data set in (tumor, normal) pairs
  sets <- list(tumor=list(), normal=list());
  for (type in colnames(pairs)) {
    idxs <- match(pairs[,type], getNames(ds));
    sets[[type]] <- extract(ds, idxs);
  }

  verbose && print(verbose, sets);
  sets;
}

############################################################################
## HISTORY
## 2010-08-12
## o Created.
############################################################################
