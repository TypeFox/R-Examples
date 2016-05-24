###########################################################################/**
# @RdocFunction .argAssertRange
#
# @title "A private function"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{x}{A @numeric @vector to be validated.}
#   \item{range}{A @numeric @vector of length two.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns \code{x} if valid, otherwise an error is thrown.
# }
#
# @author
#
# @keyword "programming"
# @keyword "internal"
#*/########################################################################### 
.argAssertRange <- function(x, range, ...) {
  name <- substitute(x);
  r <- range(x);

  if (any(is.na(r))) {
    stop("Argument '", name, "' has missing values.");
  }

  if (r[1] < range[1] || r[2] > range[2]) {
    stop("Argument '", name, "' is out of range [", range[1], ",", range[2], "]: [", r[1], ",", r[2], "]");
  }

  x;
} # .argAssertRange()

##############################################################################
# HISTORY:
# 2008-08-20
# o Created to avoid the dependancy on R.utils::Arguments.
##############################################################################
