###########################################################################/**
# @RdocFunction pseudoinverse
#
# @title "Calculates the pseudo inverse of a matrix"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{A @numeric @matrix.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @matrix.
# }
#
# \seealso{
#   Internally @see "base::svd" is used.
# }
#
# @keyword "internal"
#*/###########################################################################
pseudoinverse <- function(X, ...) {
  svd <- svd(X)
  d <- svd$d

  if(length(d) == 0L) {
    array(0, dim=dim(X)[2:1])
  } else {
    svd$v %*% (1/d * t(svd$u))
  }
} # pseudoinverse()

############################################################################
# HISTORY:
# 2009-03-24 [HB]
# o Added Rdoc comments.
# o Tidied up code. Minor minor speed up.
# 2008-xx-xx
# o Created.
############################################################################
