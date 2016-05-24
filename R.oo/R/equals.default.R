###########################################################################/**
# @RdocDefault equals
#
# @title "Compares an object with another"
#
# \description{
#  @get "title" and returns @TRUE if they are equal.
#  The equal property must be
#
#  1) \emph{reflexive}, i.e. \code{equals(o1,o1)} should be @TRUE.
#
#  2) \emph{symmetric}, i.e. \code{equals(o1,o2)} is @TRUE if and only
#  if \code{equals(o2,o1)} is @TRUE.
#
#  3) \emph{transitive}, i.e. \code{equals(o1,o2)} is @TRUE and
#  \code{equals(o2,o3)} is @TRUE, then \code{equals(o1,o3)} should
#  be @TRUE.
#
#  5) \emph{consistent}, i.e. \code{equals(o1,o2)} should return the same
#  result on multiple invocations as long as nothing has changed.
#
#  6) \code{equals(o1,}@NULL\code{)} should return @FALSE, unless
#  \code{o1} is also @NULL.
#
#  By default @see "base::identical" is used.
# }
#
# @synopsis
#
# \arguments{
#   \item{object, other}{Objects to be compared.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE if the objects are equal, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @see "base::identical".
# }
#
# @keyword attribute
# @keyword utilities
# @keyword internal
#*/###########################################################################
setMethodS3("equals", "default", function(object, other, ...) {
  (is.null(object) && is.null(other)) || identical(object, other);
})



############################################################################
# HISTORY:
# 2005-02-15
# o Added arguments '...' in order to match any generic functions.
# 2004-10-18
# o Added Rdoc comments.
# 2002-12-08
# o Created because it was needed for convience in the R.util::Tree class.
############################################################################
