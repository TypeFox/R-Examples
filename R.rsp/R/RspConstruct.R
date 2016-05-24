###########################################################################/**
# @RdocClass RspConstruct
#
# @title "The RspConstruct class"
#
# \description{
#  @classhierarchy
#
#  An RspConstruct object represents an RSP construct, which can either be
#  (i) an RSP text (a plain text section), (ii) an RSP comment,
#  (iii) an RSP preprocessing directive, or (iv) an RSP expression.
# }
#
# @synopsis
#
# \arguments{
#   \item{object}{A R object.}
#   \item{...}{Arguments passed to @see "RspObject".}
#   \item{comment}{An optional @character string.}
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
setConstructorS3("RspConstruct", function(object=character(), ..., comment=NULL) {
  this <- extend(RspObject(object, ...), "RspConstruct");
  attr(this, "#comment") <- comment;
  this;
})


#########################################################################/**
# @RdocMethod getInclude
# @alias getInclude.RspText
# @alias getInclude.RspCodeChunk
# @alias getInclude.RspVariableDirective
#
# @title "Checks whether an RSP construct will include text to the output or not"
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
#  Returns @TRUE of @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getInclude", "RspConstruct", function(object, ...) {
  FALSE;
})



#########################################################################/**
# @RdocMethod getComment
#
# @title "Gets the comment of an RSP construct"
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
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getComment", "RspConstruct", function(object, ...) {
  getAttribute(object, "#comment");
})


#########################################################################/**
# @RdocMethod getSuffixSpecs
#
# @title "Gets the suffix specifications"
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
#  Returns a trimmed @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getSuffixSpecs", "RspConstruct", function(object, ...) {
  specs <- attr(object, "suffixSpecs");
  if (is.null(specs)) return(NULL);
  specs <- trim(specs);
##  specs <- gsub("^\\[[ \t\v]*", "", specs);
##  specs <- gsub("[ \t\v]*\\]$", "", specs);
  specs;
})



#########################################################################/**
# @RdocMethod "asRspString"
# @alias asRspString.RspCode
# @alias asRspString.RspCodeChunk
# @alias asRspString.RspComment
# @alias asRspString.RspDirective
# @alias asRspString.RspDocument
# @alias asRspString.RspText
# @alias asRspString.RspUnParsedDirective
# @alias asRspString.RspUnparsedDirective
#
# @title "Recreates an RSP string from an RspConstruct"
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
#  Returns an @see "RspString".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("asRspString", "RspConstruct", function(object, ...) {
  throw(sprintf("Do not know how to construct an RSP string from %s: %s", class(object)[1L], capture.output(str(object))));
})



##############################################################################
# HISTORY:
# 2013-03-26
# o Added getInclude() to RspConstruct that defaults to FALSE.
# 2013-03-15
# o Now RspConstruct extends RspObject.
# 2013-02-22
# o Added RspUnparsedExpression.
# 2013-02-11
# o Added Rdoc help.
# 2013-02-09
# o Created.
##############################################################################
