###########################################################################/**
# @RdocDefault objectSize
#
# @title "Gets the size of the object in bytes"
#
# \description{
#   @get "title".
#   This method is just a wrapper for @see "utils::object.size".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "utils::object.size".}
# }
#
# \value{
#  Returns an @integer.
# }
#
# @author
#
# \seealso{
#   Internally @see "utils::object.size".
# }
#
# \keyword{attribute}
# \keyword{utilities}
#*/###########################################################################
setMethodS3("objectSize", "default", function(...) {
  args <- list(...)
  args$.scannedEnvs <- NULL ## Used by objectSize() for environment
  do.call(object.size, args=args)
})


############################################################################
# HISTORY:
# 2005-03-23
# o Created for completeness. Now objectSize() can be used for all kind
#   of objects.
############################################################################
