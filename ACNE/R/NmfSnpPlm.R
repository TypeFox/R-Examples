###########################################################################/**
# @RdocClass NmfSnpPlm
#
# @title "The NmfSnpPlm class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "NmfPlm".}
#   \item{mergeStrands}{If @TRUE, the sense and the anti-sense strands are
#      fitted together, otherwise separately.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/########################################################################### 
setConstructorS3("NmfSnpPlm", function(..., mergeStrands=FALSE) {
  extend(NmfPlm(...), c("NmfSnpPlm", uses(SnpPlm())),
    mergeStrands=mergeStrands
  )
})


setMethodS3("getAsteriskTags", "NmfSnpPlm", function(this, collapse=NULL, ...) {
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
# 2009-03-24 [HB]
# o Added Rdoc comments.
# 2009-01-28 [HB]
# o Dropped argument 'prevpaf' from constructor.  The new ProbeLevelModel
#   class will support this.
# 2008-07-03 [MO]
# o Created.
############################################################################
