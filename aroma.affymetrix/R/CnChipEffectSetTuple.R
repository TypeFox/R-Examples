setConstructorS3("CnChipEffectSetTuple", function(dsList=list(), ..., .setClass="CnChipEffectSet") {
  # Nothing do to?
  if (inherits(dsList, "CnChipEffectSetTuple")) {
    return(dsList);
  }

  extend(ChipEffectSetTuple(dsList, ..., .setClass=.setClass), c("CnChipEffectSetTuple", uses("CopyNumberDataSetTuple")));
})



setMethodS3("as.CnChipEffectSetTuple", "CnChipEffectSetTuple", function(this, ...) {
  # Nothing do to
  this;
})


setMethodS3("as.CnChipEffectSetTuple", "default", function(this, ...) {
  CnChipEffectSetTuple(this, ...);
})


setMethodS3("hasAlleleBFractions", "CnChipEffectSetTuple", function(this, ...) {
  cesList <- getSets(this);
  res <- sapply(cesList, FUN=hasAlleleBFractions);

  # Sanity check
  if (length(unique(res)) != 1) {
    throw("Inconsistent data sets: The ", class(this)[1], " contains data sets where some have combined the alleles and some have not.");
  }

  res <- res[1];
  res;
})

setMethodS3("hasStrandiness", "CnChipEffectSetTuple", function(this, ...) {
  cesList <- getSets(this);
  res <- sapply(cesList, FUN=hasStrandiness);

  # Sanity check
  if (length(unique(res)) != 1) {
    throw("Inconsistent data sets: The ", class(this)[1], " contains data sets where some have data by strand and some have not.");
  }

  res <- res[1];
  res;
})



##############################################################################
# HISTORY:
# 2009-11-16
# o Created.
##############################################################################
