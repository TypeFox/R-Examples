###########################################################################/**
# @RdocClass RspException
# @alias RspParseException
# @alias RspPreprocessingException
#
# @title "The RspException class"
#
# \description{
#  @classhierarchy
#
#  An RspException is an @see "R.oo::Exception" that is thrown during
#  the processing of an RSP document.
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
setConstructorS3("RspException", function(...) {
  extend(Exception(...), "RspException");
})


setConstructorS3("RspParseException", function(...) {
  extend(RspException(...), "RspParseException");
})


setConstructorS3("RspPreprocessingException", function(..., item=NULL) {
  extend(RspException(...), "RspPreprocessingException",
    item = item
  );
})

setMethodS3("getMessage", "RspPreprocessingException", function(this, ...) {
  ## The following is not possible due to bug in R.oo 1.13.0:
  ##  msg <- NextMethod("getMessage");
  msg <- this$.msg;

  item <- this$item;
  if (!is.null(item)) {
    itemStr <- asRspString(item);
    itemStr <- sprintf(" (%s)", itemStr);
  } else {
    itemStr <- "";
  }
  sprintf("An error occured while preprocessing RSP directive%s: %s", itemStr, msg);
}, createGeneric=FALSE)

setMethodS3("getItem", "RspPreprocessingException", function(this, ...) {
  this$item;
})




##############################################################################
# HISTORY:
# 2013-03-17
# o Created.
##############################################################################
