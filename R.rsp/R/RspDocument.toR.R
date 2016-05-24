#########################################################################/**
# @set "class=RspDocument"
# @RdocMethod toR
#
# @title "Translates the RSP document into R source code"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{factory}{A @see "RspSourceCodeFactory".}
#   \item{...}{Optional arguments passed to \code{toSourceCode()} for
#              the @see "RspSourceCodeFactory".}
# }
#
# \value{
#  Returns the R source code as an @see "RspRSourceCode".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("toR", "RspDocument", function(object, factory=RspRSourceCodeFactory(), ...) {
  # Argument 'factory':
  factory <- Arguments$getInstanceOf(factory, "RspSourceCodeFactory");

  toSourceCode(factory, object, ...);
}) # toR()


#########################################################################/**
# @RdocMethod evaluate
#
# @title "Parses, translates, and evaluates the RSP document"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{envir}{The @environment where the RSP document is evaluated.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns the last evaluated expression, iff any.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("evaluate", "RspDocument", function(object, envir=parent.frame(), ...) {
  code <- toR(object);
  process(code, envir=envir, ...);
}, createGeneric=FALSE)



##############################################################################
# HISTORY:
# 2013-02-09
# o Created.
##############################################################################
