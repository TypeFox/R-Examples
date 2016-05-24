###########################################################################/**
# @RdocClass RspSourceCode
#
# @title "The RspSourceCode class"
#
# \description{
#  @classhierarchy
#
#  An RspSourceCode object is a @character @vector holding RSP generated
#  source code for a particular programming language.
# }
#
# @synopsis
#
# \arguments{
#   \item{code}{@character @vector.}
#   \item{...}{Additional arguments passed to the @see "RspProduct"
#     constructor.}
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
setConstructorS3("RspSourceCode", function(code=character(), ...) {
  extend(RspProduct(code, ...), "RspSourceCode");
})


setMethodS3("print", "RspSourceCode", function(x, ...) {
  code <- x;
  code <- paste(code, collapse="\n");
  cat(code);
  cat("\n");
})



#########################################################################/**
# @RdocMethod evaluate
#
# @title "Evaluates the source code"
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
#  Returns the last evaluated expression, iff any.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("evaluate", "RspSourceCode", abstract=TRUE, createGeneric=FALSE);



#########################################################################/**
# @RdocMethod tangle
# @alias tangle.RspRSourceCode
#
# @title "Drops all text-outputting calls from the source code"
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
#  Returns a @see "RspSourceCode" objects.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("tangle", "RspSourceCode", abstract=TRUE);


#########################################################################/**
# @RdocMethod tidy
# @alias tidy.RspRSourceCode
#
# @title "Tidy up the RSP source code"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{format}{A @character string specifying how the source code
#     should be tidied.}
#   \item{collapse}{How source code lines should be collapsed.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @RspSourceCode of the same class as the input source code.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("tidy", "RspSourceCode", function(object, format=c("asis"), collapse="\n", ...) {
  # Argument 'format':
  format <- match.arg(format);

  # Record attributes
  attrs <- attributes(object);

  # Collapse?
  if (!is.null(collapse)) {
    object <- paste(object, collapse=collapse);
  }

  # Restore attributes (if lost above)
  attributes(object) <- attrs;

  object;
})


##############################################################################
# HISTORY:
# 2013-02-16
# o Now RspSourceCode extends RspProduct,
# o Renamed SourceCode to RspSourceCode.
# 2013-02-14
# o Added tangle() for SourceCode.
# 2013-02-13
# o Added getType() for SourceCode.
# 2013-02-11
# o Added Rdoc help.
# 2013-02-09
# o Created.
##############################################################################
