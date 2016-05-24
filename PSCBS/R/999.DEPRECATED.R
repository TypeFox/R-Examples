## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## DEFUNCT
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Defunct since v0.41.0 (2015-03-30)
setMethodS3("bootstrapDHByRegion", "PairedPSCBS", function(fit, B=100, statsFcn=function(x) quantile(x, probs=c(0.025, 0.050, 0.95, 0.975)), by=c("betaTN", "betaT"), ..., force=FALSE, verbose=FALSE) {
  .Defunct("bootstrapTCNandDHByRegion");
}, deprecated=TRUE) # bootstrapDHByRegion()



##############################################################################
# HISTORY
# 2015-02-22
# o CLEANUP: bootstrapDHByRegion() is defunct (was deprecated since 2013).
# 2013-04-20
# o CLEANUP: Formally deprecated bootstrapDHByRegion().
# 2013-01-15
# o Now bootstrapDHByRegion() uses the params$avgDH estimator, iff given.
# 2012-02-24
# o Added argument 'force' to bootstrapDHByRegion().
# 2011-06-14
# o Updated code to recognize new column names.
# 2010-11-22
# o DEPRECATED: bootstrapDHByRegion() should no longer be used.
# 2010-11-03 [HB]
# o ROBUSTNESS: Now bootstrapDHByRegion() uses resample() of R.utils.
# 2010-11-01 [HB]
# o Now bootstrapDHByRegion() estimates more quantiles.
# o BUG FIX: bootstrapDHByRegion() would give an error if only a single
#   quantile was requested.
# o BUG FIX: bootstrapDHByRegion() would give "Error in if (nbrOfUnits >
#   segJJ[, "dh.num.mark"]) { : missing value where TRUE/FALSE needed" when
#   'dh.num.mark' was NA.
# 2010-10-25 [HB]
# o BUG FIX: Now bootstrapDHByRegion() for PairedPSCBS handles the rare case
#   when markers with the same positions are split in two different segments.
# o BUG FIX: bootstrapDHByRegion() for PairedPSCBS would bootstrap from the
#   incorrect set of loci when the DH region contained only one locus.
# o BUG FIX: bootstrapDHByRegion() for PairedPSCBS would bootstrap from the
#   incorrect set of loci if more than one chromosome was available.
# 2010-09-16 [HB]
# o Added bootstrapDHByRegion(), which is what is used by paired PSCBS.
# o Created.
##############################################################################
