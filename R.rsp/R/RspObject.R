###########################################################################/**
# @RdocClass RspObject
#
# @title "The abstract RspObject class"
#
# \description{
#  @classhierarchy
#
#  An RspObject represents an instance a specific RSP class.
# }
#
# @synopsis
#
# \arguments{
#   \item{value}{An R object.}
#   \item{attrs}{RSP attributes as a named @list, e.g. \code{type},
#      \code{language}, and \code{source}.}
#   \item{...}{Additional named RSP attributes.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setConstructorS3("RspObject", function(value=NA, attrs=list(), ...) {
  # Argument 'attrs':
  if (!is.list(attrs)) {
    throw("Argument 'attrs' is not a list: ", mode(attrs)[1L]);
  }

  # Argument '...':
  userAttrs <- list(...);


  this <- extend(value, "RspObject");
  this <- setAttributes(this, attrs);
  this <- setAttributes(this, userAttrs);
  this;
})



#########################################################################/**
# @RdocMethod print
# @alias print.RspDocument
# @alias print.RspFileProduct
# @alias print.RspProduct
# @alias print.RspSourceCode
# @alias print.RspString
# @alias print.RspStringProduct
#
# @title "Prints a summary of an RSP object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("print", "RspObject", function(x, ...) {
  s <- NextMethod("print");
  s <- c(sprintf("%s:", class(x)[1L]), s);
  s <- paste(s, collapse="\n");
  cat(s, "\n", sep="");
}, protected=TRUE)



##############################################################################
# HISTORY:
# 2015-02-04
# o CLEANUP: Turned all attribute methods into "default" methods.
# 2013-??-??
# o Created from RspString.R.
##############################################################################
