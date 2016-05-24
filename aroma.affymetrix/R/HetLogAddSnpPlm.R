###########################################################################/**
# @RdocClass HetLogAddSnpPlm
#
# @title "The HetLogAddSnpPlm class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "HetLogAddPlm".}
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
setConstructorS3("HetLogAddSnpPlm", function(..., mergeStrands=FALSE) {
  extend(HetLogAddPlm(...), c("HetLogAddSnpPlm", uses(SnpPlm())),
    mergeStrands = mergeStrands
  )
})


setMethodS3("getAsteriskTags", "HetLogAddSnpPlm", function(this, collapse=NULL, ...) {
  # Returns 'HLA[,<flavor>]'
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
# o Added getAsteriskTag() for HetLogAddSnpPlm.
# 2007-10-06
# o Created.
############################################################################
