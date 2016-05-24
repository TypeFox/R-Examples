###########################################################################/**
# @RdocClass RawCopyNumberModel
#
# @title "The RawCopyNumberModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents an identity copy-number model which returns
#  the input as is.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Passed to the constructor of the superclass.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("RawCopyNumberModel", function(...) {
  extend(CopyNumberChromosomalModel(...), "RawCopyNumberModel");
})

setMethodS3("getAsteriskTags", "RawCopyNumberModel", function(this, ...) {
  "";
}, protected=TRUE)

setMethodS3("getSetTag", "RawCopyNumberModel", function(this, ...) {
  "raw";
}, protected=TRUE)


##############################################################################
# HISTORY:
# 2007-10-18
# o Created.
##############################################################################
