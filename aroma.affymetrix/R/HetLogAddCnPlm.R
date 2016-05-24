###########################################################################/**
# @RdocClass HetLogAddCnPlm
#
# @title "The HetLogAddCnPlm class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "HetLogAddCnPlm".}
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
setConstructorS3("HetLogAddCnPlm", function(..., combineAlleles=FALSE) {
  extend(HetLogAddSnpPlm(...), c("HetLogAddCnPlm", uses(CnPlm())),
    combineAlleles = combineAlleles
  )
})


setMethodS3("getAsteriskTags", "HetLogAddCnPlm", function(this, collapse=NULL, ...) {
  # Returns 'HLA[,<flavor>][,+-]'
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Add class-specific parameter tags
  if (this$combineAlleles)
    tags <- c(tags, "A+B");

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2007-12-06
# o Added getAsteriskTag() for HetLogAddCnPlm.
# 2007-10-06
# o Created.
############################################################################
