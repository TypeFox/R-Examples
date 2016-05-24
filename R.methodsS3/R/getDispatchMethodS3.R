###########################################################################/**
# @RdocDefault getDispatchMethodS3
#
# @title "Gets the S3 method that a generic function would call"
#
# \description{
#  @get "title" according to an S3 @see "base::class" @vector.
# }
#
# @synopsis
#
# \arguments{
#   \item{methodName}{A @character string specifying the name of a
#     generic function.}
#   \item{classNames}{A @character @vector of @see "base::class" names.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @function, or throws an exception if not found.
# }
#
# \seealso{
#   @see "findDispatchMethodsS3".
# }
#
# @author
#
# @keyword programming
# @keyword methods
# @keyword internal
#*/###########################################################################
setMethodS3("getDispatchMethodS3", "default", function(methodName, classNames, ...) {
  res <- findDispatchMethodsS3(methodName, classNames, firstOnly=TRUE, ...);
  if (length(res) == 0) {
    throw(sprintf("No method %s() for this class structure: %s", methodName, paste(classNames, collapse=", ")));
  }

  res[[1]]$fcn;
}, private=TRUE)


############################################################################
# HISTORY:
# 2010-12-02
# o Added Rdoc comments.
# o Made getDispatchMethodS3() a default method.
# 2009-11-20
# o Added getDispatchMethodS3().
############################################################################
