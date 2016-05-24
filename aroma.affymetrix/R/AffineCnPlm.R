###########################################################################/**
# @RdocClass AffineCnPlm
#
# @title "The AffineCnPlm class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffineSnpPlm".}
#   \item{combineAlleles}{If @FALSE, allele A and allele B are treated
#      seperately, otherwise together.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("AffineCnPlm", function(..., combineAlleles=FALSE) {
  extend(AffineSnpPlm(...), c("AffineCnPlm", uses(CnPlm())),
    combineAlleles = combineAlleles
  )
})


setMethodS3("getAsteriskTags", "AffineCnPlm", function(this, collapse=NULL, ...) {
  # Returns 'AFF[,<flavor>][,+-]'
  tags <- NextMethod("getAsteriskTags", collapse=collapse);

  # Add class specific parameter tags
  if (this$combineAlleles)
    tags <- c(tags, "A+B");

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2007-12-06
# o Added getAsteriskTag() for AffineCnPlm.
# 2007-01-07
# o Created from MbeiCnPlm.R.
############################################################################
