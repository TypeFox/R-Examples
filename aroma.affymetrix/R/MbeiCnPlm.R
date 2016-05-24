###########################################################################/**
# @RdocClass MbeiCnPlm
#
# @title "The MbeiCnPlm class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "MbeiSnpPlm".}
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
setConstructorS3("MbeiCnPlm", function(..., combineAlleles=FALSE) {
  extend(MbeiSnpPlm(...), c("MbeiCnPlm", uses(CnPlm())),
    combineAlleles = combineAlleles
  )
})


setMethodS3("getAsteriskTags", "MbeiCnPlm", function(this, collapse=NULL, ...) {
  # Returns 'MBEI[,<flavor>][,+-]'
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
# o Added getAsteriskTag() for MbeiCnPlm.
# 2006-09-12
# o Recreated.
############################################################################
