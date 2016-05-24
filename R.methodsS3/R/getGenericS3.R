###########################################################################/**
# @RdocDefault getGenericS3
#
# @title "Gets an S3 generic function"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{The name of the generic function.}
#   \item{envir}{The @environment from which the search for the
#     generic @function is done.}
#   \item{inherits}{A @logical specifying whether the enclosing frames
#     should be searched or not.}
#   \item{...}{Not used.}
# }
#
# \seealso{
#   @see "setGenericS3".
#   @see "getMethodS3".
#   @see "isGenericS3".
# }
#
# @author
#
# @keyword programming
# @keyword methods
#*/###########################################################################
setMethodS3("getGenericS3", "default", function(name, envir=parent.frame(), inherits=TRUE, ...) {
  fcn <- .findFunction(name, envir=envir, inherits=inherits)$fcn;
  if (is.null(fcn)) {
    throw("No such function found: ", name);
  } else if (!isGenericS3(fcn)) {
    throw("The function found is not an S3 generic function: ", name);
  }
  fcn;
})



############################################################################
# HISTORY:
# 2013-10-06
# o Now getGenericS3() uses .findFunction().
# 2013-10-05
# o Added argument 'inherits' to getGenericS3().
# 2010-09-18
# o BUG FIX: getGenericS3() failed to locate generic functions created
#   in the global enviroment.
# 2008-05-08
# o Added getGenericS3().
############################################################################
