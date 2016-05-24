###########################################################################/**
# @RdocClass AvgSnpPlm
#
# @title "The AvgSnpPlm class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AvgPlm".}
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
setConstructorS3("AvgSnpPlm", function(..., mergeStrands=FALSE) {
  extend(AvgPlm(...), c("AvgSnpPlm", uses(SnpPlm())),
    mergeStrands = mergeStrands
  )
})


setMethodS3("getAsteriskTags", "AvgSnpPlm", function(this, collapse=NULL, ...) {
  # Returns 'AVG[,<flavor>]'
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
# o Added getAsteriskTag() for AvgSnpPlm.
# 2007-09-08
# o Created from the MbeiSnpPlm.R.
############################################################################
