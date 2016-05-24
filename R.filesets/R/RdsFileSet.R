###########################################################################/**
# @RdocClass RdsFileSet
# @alias byPath.RdsFileSet
#
# @title "The RdsFileSet class"
#
# \description{
#  @classhierarchy
#
#  An RdsFileSet object represents a set of @see "RdsFile":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("RdsFileSet", function(...) {
  extend(GenericDataFileSet(...), "RdsFileSet");
})


setMethodS3("byPath", "RdsFileSet", function(static, ..., pattern="[.]rds$") {
  NextMethod("byPath", pattern=pattern);
})



###########################################################################
# HISTORY:
# 2013-11-20
# o Created.
############################################################################
