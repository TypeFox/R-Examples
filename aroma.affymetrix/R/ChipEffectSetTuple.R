setConstructorS3("ChipEffectSetTuple", function(dsList=list(), ..., .setClass="ChipEffectSet") {
  # Nothing do to?
  if (inherits(dsList, "ChipEffectSetTuple")) {
    return(dsList);
  }

  extend(AffymetrixCelSetTuple(dsList, ..., .setClass=.setClass), "ChipEffectSetTuple");
})


setMethodS3("getFullNames", "ChipEffectSetTuple", function(this, ..., exclude="chipEffects") {
  NextMethod("getFullNames", exclude=exclude);
})


setMethodS3("as.ChipEffectSetTuple", "ChipEffectSetTuple", function(this, ...) {
  # Nothing do to
  this;
})


setMethodS3("as.ChipEffectSetTuple", "default", function(this, ...) {
  ChipEffectSetTuple(this, ...);
})


##############################################################################
# HISTORY:
# 2009-11-18
# o Added as.ChipEffectSetTuple().
# 2007-03-19
# o Created.
##############################################################################
