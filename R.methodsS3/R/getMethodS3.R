###########################################################################/**
# @RdocDefault getMethodS3
#
# @title "Gets an S3 method"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{The name of the method.}
#   \item{class}{The class of the method.}
#   \item{envir}{The @environment from which the search for the
#     S3 method is done.}
#   \item{...}{Not used.}
# }
#
# \seealso{
#   This is just a conveniency wrapper around @see "utils::getS3method"
#   that have arguments consistent with @see "setMethodS3".
#   @see "getGenericS3".
# }
#
# @author
#
# @keyword programming
# @keyword methods
#*/###########################################################################
setMethodS3("getMethodS3", "default", function(name, class="default", envir=parent.frame(), ...) {
  args <- list(name, class=class, optional=FALSE);
  do.call(getS3method, args, envir=envir);
})



############################################################################
# HISTORY:
# 2011-11-17
# o CLEANUP: Dropped example(getMethodS3), which was for setMethodS3().
# 2010-09-18
# o BUG FIX: getMethodS3() failed to locate S3 methods created in the
#   global enviroment.
# 2008-05-08
# o Added getMethodS3().
############################################################################
