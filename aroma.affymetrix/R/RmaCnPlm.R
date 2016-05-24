###########################################################################/**
# @RdocClass RmaCnPlm
#
# @title "The RmaCnPlm class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RmaSnpPlm".}
#   \item{combineAlleles}{If @FALSE, allele A and allele B are treated
#      seperately, otherwise together.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#   TO DO.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("RmaCnPlm", function(..., combineAlleles=FALSE) {
  extend(RmaSnpPlm(...), c("RmaCnPlm", uses(CnPlm())),
    combineAlleles = combineAlleles
  )
})


setMethodS3("getAsteriskTags", "RmaCnPlm", function(this, collapse=NULL, ...) {
  # Returns 'RMA[,<flavor>][,+-]'
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Add class specific parameter tags
  if (this$combineAlleles)
    tags <- c(tags, "A+B");

  # Collapse
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2007-12-06
# o Added getAsteriskTag() for RmaSnpPlm.
# 2006-09-10
# o Recreated.
############################################################################
