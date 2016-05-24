###########################################################################/**
# @RdocClass RspCutDirective
# @alias RspEndcutDirective
# @alias RspCopyDirective
# @alias RspEndcopyDirective
# @alias RspPasteDirective
#
# @title "The RspCopyDirective class"
#
# \description{
#  @classhierarchy
#
#  An RspCutDirective and RspPasteDirective are
#  @see "RspDirective" that will copy, cut, and paste a set of
#  @see "RspConstruct":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{value}{A @character string.}
#   \item{...}{Additional arguments passed to @see "RspDirective".}
#   \item{.validate}{A @logical specifying whether the attributes should
#    be validated or not.}
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
setConstructorS3("RspCutDirective", function(value="cut", ..., .validate=TRUE) {
  # Argument '.validate':
  if (missing(.validate)) .validate <- !missing(value);

  # Argument 'value':
  value <- match.arg(value, choices=c("cut", "copy"));

  this <- extend(RspDirective(value, ...), "RspCutDirective")
  if (.validate) {
    requireAttributes(this, c("name"));
  }
  this;
})

setConstructorS3("RspEndcutDirective", function(value="endcut", ...) {
  value <- match.arg(value, choices=c("endcut", "endcopy"));
  extend(RspDirective(value, ...), "RspEndcutDirective")
})

setConstructorS3("RspCopyDirective", function(value="copy", ..., .validate=TRUE) {
  # Argument '.validate':
  if (missing(.validate)) .validate <- !missing(value);

  extend(RspCutDirective(value, ..., .validate=.validate), "RspCopyDirective")
})

setConstructorS3("RspEndcopyDirective", function(value="endcopy", ...) {
  extend(RspEndcutDirective(value, ...), "RspEndcopyDirective")
})


setConstructorS3("RspPasteDirective", function(value="paste", ...) {
  this <- extend(RspDirective(value, ...), "RspPasteDirective")
  if (!missing(value)) {
    requireAttributes(this, c("name"));
  }
  this;
})

##############################################################################
# HISTORY:
# 2014-07-02
# o Added RspCopyDirective and RspEndcopyDirective.
# 2014-07-01
# o Created.
##############################################################################
