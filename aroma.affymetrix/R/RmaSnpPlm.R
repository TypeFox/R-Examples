###########################################################################/**
# @RdocClass RmaSnpPlm
#
# @title "The RmaSnpPlm class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RmaPlm".}
#   \item{mergeStrands}{If @TRUE, the sense and the anti-sense strands are
#      fitted together, otherwise separately.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("RmaSnpPlm", function(..., mergeStrands=FALSE) {
  extend(RmaPlm(...), c("RmaSnpPlm", uses(SnpPlm())),
    mergeStrands = mergeStrands
  )
})


setMethodS3("getAsteriskTags", "RmaSnpPlm", function(this, collapse=NULL, ...) {
  # Returns 'RMA[,<flavor>]'
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Add class specific parameter tags
  if (!this$mergeStrands)
    tags <- c(tags, "+-");

  # Collapse
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2007-12-06
# o Added getAsteriskTag() for RmaSnpPlm.
# 2006-09-11
# o Simple benchmarking [Thinkpad A31]: Fitting 1000 units (with merged
#   strands) across 22 arrays (100K Xba) takes in total 114 sec, that is,
#   5.1ms/unit/array.
#   For all 59015 SNPs it takes ~5.0min/array or ~112min/22 arrays.
#   4545 units and 22 arrays: 60s to read all data, 50s to fit the model,
#   30s to store probe affinities, and 120s to store chip-effects.
#   In total 274s, that is, 2.7ms/unit/array.
#   We are still spending [(60+120)/274 =] 65% on I/O.
#   For all 59015 SNPs it takes ~2.7min/array or ~60min/22 arrays.
# o The fit function now returns all required fields.
# 2006-08-25
# o Created from the corresponding Li & Wong model.
############################################################################
