###########################################################################/**
# @RdocClass HtmlRspLanguage
#
# @title "The HtmlRspLanguage class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of the @see "RspLanguage".}
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
setConstructorS3("HtmlRspLanguage", function(...) {
  extend(RspLanguage(language="html", ...), "HtmlRspLanguage")
})


setMethodS3("getComment", "HtmlRspLanguage", function(object, ...) {
  s <- paste(..., collapse="\n", sep="");
  s <- paste("<!-- ", s, " -->", sep="");
  s;
})


setMethodS3("escape", "HtmlRspLanguage", function(object, ...) {
  s <- paste(..., collapse="\n", sep="");
  s <- gsub("<", "&lt;", s);
  s <- gsub(">", "&gt;", s);
  s;
})


setMethodS3("getVerbatim", "HtmlRspLanguage", function(object, ..., newline=NULL) {
  s <- NextMethod("getVerbatim", newline="");
  if (is.null(newline))
    newline <- getNewline(object);
  if (is.character(newline)) {
    s <- gsub("\n\r|\r\n|\r", "\n", s);
    s <- gsub("\n", newline, s);
  }
  s;
})


##############################################################################
# HISTORY:
# 2005-08-01
# o Added Rdoc comments.
# 2005-07-29
# o Created.
##############################################################################
