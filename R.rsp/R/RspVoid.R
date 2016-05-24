###########################################################################/**
# @RdocClass RspVoid
#
# @title "The RspVoid class"
#
# \description{
#  @classhierarchy
#
#  An RspVoid is an @see "RspConstruct" that contains nothing and
#  outputs nothing.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
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
setConstructorS3("RspVoid", function(...) {
  extend(RspConstruct(), "RspVoid")
})

setMethodS3("asRspString", "RspVoid", function(object, ...) {
  RspString()
})

##############################################################################
# HISTORY:
# 2015-02-04
# o Added Rdoc.
# 2014-10-19
# o Created.
##############################################################################
