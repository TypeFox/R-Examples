###########################################################################/**
# @RdocClass AvgCnPlm
#
# @title "The AvgCnPlm class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AvgSnpPlm".}
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
setConstructorS3("AvgCnPlm", function(..., combineAlleles=FALSE) {
  extend(AvgSnpPlm(...), c("AvgCnPlm", uses(CnPlm())),
    combineAlleles = combineAlleles
  )
})


setMethodS3("getAsteriskTags", "AvgCnPlm", function(this, collapse=NULL, ...) {
  # Returns 'AVG[,<flavor>][,+-]'
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
# o Added getAsteriskTag() for AvgSnpPlm.
# 2007-09-08
# o Created from MbeiCnPlm.R.
############################################################################
