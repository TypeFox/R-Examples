###########################################################################/**
# @RdocClass SmoothRmaModel
#
# @title "The SmoothRmaModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Chromosomal Smoothing Robust Multichip Analysis
#  method.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#              @see "SmoothMultiarrayModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("SmoothSaModel", function(...) {
  extend(SmoothMultiarrayModel(...), "SmoothSaModel");
})


setMethodS3("getAsteriskTags", "SmoothSaModel", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Replace first tags
  tags[1] <- "SA";

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }

  tags;
}, protected=TRUE)


setMethodS3("getFitUnitGroupFunction", "SmoothSaModel", function(this, ...) {
  smoothWSA;
}, protected=TRUE)


##############################################################################
# HISTORY:
# 2007-09-20
# o Created.
##############################################################################
