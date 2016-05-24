###########################################################################/**
# @RdocClass RDataFileSet
# @alias byPath.RDataFileSet
#
# @title "The RDataFileSet class"
#
# \description{
#  @classhierarchy
#
#  An RDataFileSet object represents a set of @see "RDataFile":s.
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
setConstructorS3("RDataFileSet", function(...) {
  extend(GenericDataFileSet(...), "RDataFileSet")
})


setMethodS3("byPath", "RDataFileSet", function(static, ..., pattern="[.](RData|Rdata|rdata)$") {
  NextMethod("byPath", pattern=pattern)
})


###########################################################################
# HISTORY:
# 2015-01-12
# o Created.
############################################################################
