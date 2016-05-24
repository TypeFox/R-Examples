###########################################################################/**
# @RdocClass MbeiSnpPlm
#
# @title "The MbeiSnpPlm class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "MbeiPlm".}
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
setConstructorS3("MbeiSnpPlm", function(..., mergeStrands=FALSE) {
  extend(MbeiPlm(...), c("MbeiSnpPlm", uses(SnpPlm())),
    mergeStrands = mergeStrands
  )
})


setMethodS3("getAsteriskTags", "MbeiSnpPlm", function(this, collapse=NULL, ...) {
  # Returns 'MBEI[,<flavor>]'
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
# o Added getAsteriskTag() for MbeiSnpPlm.
# 2006-09-11
# o Created from the RmaSnpPlm.
############################################################################
