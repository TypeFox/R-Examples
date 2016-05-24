###########################################################################/**
# @RdocDefault throw
#
# @title "Throws an Exception"
#
# \description{
#  Throws an exception similar to \code{stop()}, but with support for
#  @see "Exception" classes. The first argument (\code{object}) is by 
#  default pasted together with other arguments (@...) and with seperator
#  \code{sep=""}.  For instance, to throw an exception, write
#
#    \code{throw("Value out of range: ", value, ".")}.
#
#  which is short for
#
#    \code{throw(Exception("Value out of range: ", value, "."))}.
#
#  Note that \code{throw()} can be defined for classes inheriting
#  @see "Exception", which can then be caught (or not)
#  using \code{\link[base:conditions]{tryCatch}()}.
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
# \examples{
#   rbern <- function(n=1, prob=1/2) {
#     if (prob < 0 || prob > 1)
#       throw("Argument 'prob' is out of range: ", prob)
#     rbinom(n=n, size=1, prob=prob)
#   }
#
#   rbern(10, 0.4)
#   # [1] 0 1 0 0 0 1 0 0 1 0
#   tryCatch(rbern(10, 10*0.4),
#     error=function(ex) {}
#   )
# }
#
# @author
#
# \seealso{
#   See the @see "Exception" class for more detailed information.
# }
#
# \keyword{error}
#*/###########################################################################
setMethodS3("throw", "default", function(...) {
  throw(Exception(...));
}, overwrite=TRUE, conflict="quiet")


## setGenericS3("throw", force=TRUE);


############################################################################
# HISTORY:
# 2012-06-17
# o Override generic function throw() of R.methodsS3 with one here.
# 2012-03-08
# o Now the default throw() of R.methodsS3 is "quietly" overwritten,
#   i.e. there is no longer a warning about it when R.oo is loaded.
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
