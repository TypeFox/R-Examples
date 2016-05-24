###########################################################################/**
# @RdocDefault dimension
#
# @title "Gets the dimension of the object"
#
# \description{
#   Gets the dimension of the object similar to what \code{dim()} does,
#   but instead of @NULL it will return the length of a vector.
#   If a function is passed, @NULL is returned.
# }
#
# @synopsis
#
# \arguments{
#   \item{object}{The object for which the dimension should be obtained.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer @vector or @NULL.
# }
#
# \examples{
#   dimension(matrix(1:100, ncol=10))     # 10 10
#   dimension(1:14)                       # 14
#   dimension(data.frame(a=1:10, b=10:1)) # 10  2
#   dimension(print)                      # NULL
# }
#
# @author
#
# \seealso{
#   @see "ll.default".
#   @see "base::dim" and @see "base::length".
# }
#
# @keyword attribute
# @keyword utilities
# @keyword internal
#*/###########################################################################
setMethodS3("dimension", "default", function(object, ...) {
  if (is.function(object))
    return(NULL);
  size <- dim(object);
  if (is.null(size))
    size <- length(object);
  as.integer(size);
})




############################################################################
# HISTORY:
# 2005-02-15
# o Added arguments '...' in order to match any generic functions.
# 2002-10-17
# o Created to be used by ll().
############################################################################
