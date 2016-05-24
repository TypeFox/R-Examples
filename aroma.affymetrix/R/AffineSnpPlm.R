###########################################################################/**
# @RdocClass AffineSnpPlm
#
# @title "The AffineSnpPlm class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffinePlm".}
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
setConstructorS3("AffineSnpPlm", function(..., mergeStrands=FALSE) {
  extend(AffinePlm(...), c("AffineSnpPlm", uses(SnpPlm())),
    mergeStrands = mergeStrands
  )
})



setMethodS3("getAsteriskTags", "AffineSnpPlm", function(this, collapse=NULL, ...) {
  # Returns 'AFF[,<flavor>]'
  tags <- NextMethod("getAsteriskTags", collapse=collapse);

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
# o Added getAsteriskTag() for AffineSnpPlm.
# 2006-09-11
# o Created from the MbeiSnpPlm.
############################################################################
