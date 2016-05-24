###########################################################################/**
# @RdocDefault throw
#
# @title "Throws an exception"
#
# \description{
#  Throws an exception by calling stop().
#
#  Note that \code{throw()} can be defined for specific classes, which can
#  then be caught (or not) using \code{\link[base:conditions]{tryCatch}}().
#
#  \emph{This default function will be overridden by ditto in the \pkg{R.oo}
#  package, if that is loaded.  The latter @see "R.oo::throw" implementation
#  is fully backward compatible with this one, but the error object thrown
#  is of class @see "R.oo::Exception".}
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{One or several strings that are concatenated and collapsed
#       into on message string.}
# }
#
# \value{
#   Returns nothing.
# }
#
# @examples "../incl/throw.Rex"
#
# @author
#
# \keyword{error}
#*/###########################################################################
setMethodS3("throw", "default", function(...) {
  stop(...);
})



############################################################################
# HISTORY:
# 2005-09-17
# o Added to R.methodsS3 since it is so commonly used by my packages.
# 2005-02-20
# o Updated broken link to tryCatch().
# 2005-02-10
# o Making use of tryCatch() only.
# 2002-10-17
# o Now throw() always throws an Exception.
# 2002-05-25
# * Bug fix in Rd \examples{}. Forgot a comment.
# 2002-04-21
# * Redefined throw.default() so it takes several arguments, which are then
#   pasted together with sep="". In other words, instead of doing
#     stop(paste("bla bla", "value:", x, ".\n", sep=""))
#   one can just do
#     throw("bla bla", "value:", x, ".\n")
#   This is also a step towards the new exception model that supports
#   classes.
# * Extract the throw() functions from trycatch.R, which relies on them, but
#   the throw()'s are stand-alone.
############################################################################
