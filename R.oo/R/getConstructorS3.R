###########################################################################/**
# @RdocDefault getConstructorS3
#
# @title "Get a constructor method"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{The name of the constructor function.}
#   \item{...}{Not used.}
# }
#
# \seealso{
#   @see "setConstructorS3".
#   @see "R.methodsS3::getMethodS3".
#   @see "R.methodsS3::isGenericS3".
# }
#
# @author
#
# @keyword "programming"
# @keyword "methods"
#*/###########################################################################
setMethodS3("getConstructorS3", "default", function(name, ...) {
  # TO DO/FIX ME: This part only works when packages are attached.
  # /HB 2013-10-08
  if (!exists(name, mode="function")) {
    throw("No such function found: ", name);
  }

  fcn <- get(name, mode="function");
  if (isGenericS3(fcn)) {
    throw("The function found is an S3 generic function: ", name);
  }

  fcn;
})


############################################################################
# HISTORY:
# 2008-05-08
# o Added getConstructorS3().
############################################################################
