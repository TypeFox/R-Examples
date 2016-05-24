###########################################################################/**
# @set "class=error"
# @RdocMethod throw
#
# @title "Throws (rethrows) an object of class 'error'"
#
# \description{
#  Rethrows an 'error' object. The 'error' class was introduced in R v1.8.1
#  with the new error handling mechanisms.
# }
#
# @synopsis
#
# \arguments{
#   \item{error}{An object or class 'error'.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns nothing.
# }
#
# @author
#
# \seealso{
#   See the \code{tryCatch()} method etc.
#   See the @see "Exception" class for more detailed information.
# }
#
# \keyword{error}
#*/###########################################################################
setMethodS3("throw", "error", function(error, ...) {
  base::stop(error);
})



############################################################################
# HISTORY:
# 2012-02-29
# o CLEANUP: throw() for 'error' is now just a wrapper for stop(). 
#   Previously it had to do what stop() now does for 'condition' objects.
# 2005-02-15
# o Added arguments '...' in order to match any generic functions.
# 2004-03-03
# o Added throw() for the error class (new in R v1.8.0).
############################################################################
