getPlasqTypes <- function(cdf, ...) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  cdfGetFields <- affxparser::cdfGetFields
#  cdfMergeAlleles <- affxparser::cdfMergeAlleles


  cdf <- .applyCdfGroups(cdf, .cdfAddPlasqTypes);
  cdf <- .applyCdfGroups(cdf, cdfGetFields, c("plasqType"));
#  cdf <- .applyCdfGroups(cdf, cdfMergeAlleles);
  cdf <- lapply(cdf, FUN=.subset2, "groups");
  cdf;
}

############################################################################
# HISTORY:
# 2006-12-30
# o Created.
############################################################################
