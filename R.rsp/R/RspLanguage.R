###########################################################################/**
# @RdocClass RspLanguage
# @alias escape.HtmlRspLanguage
#
# @title "The RspLanguage class"
#
# \description{
#  @classhierarchy
#
#  An RspLanguage object specifies what the markup language of the
#  response/output document is, e.g. plain text and HTML.
#  The RspLanguage class provides methods to obtain language specific
#  strings/output such as how newlines and comments are written.
#  The RspLanguage class describes a plain text languages.  For HTML
#  see the @see "HtmlRspLanguage" subclass.
# }
#
# @synopsis
#
# \arguments{
#   \item{language}{A @character string.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
# @keyword internal
#*/###########################################################################
setConstructorS3("RspLanguage", function(language="plain", ...) {
  language <- Arguments$getCharacter(language, length=1, nchar=c(1,64));

  extend(Object(), "RspLanguage",
    language=language
  );
})



#########################################################################/**
# @RdocMethod getLanguage
#
# @title "Gets the language string"
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
#
# @keyword IO
#*/#########################################################################
setMethodS3("getLanguage", "RspLanguage", function(this, ...) {
  this$language;
})


#########################################################################/**
# @RdocMethod getNewline
#
# @title "Gets the newline string specific for a given RSP response language"
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
#
# @keyword IO
#*/#########################################################################
setMethodS3("getNewline", "RspLanguage", function(this, ...) {
  "\n";
})



#########################################################################/**
# @RdocMethod getComment
# @alias getComment.HtmlRspLanguage
#
# @title "Gets a comment string specific for a given RSP response language"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{R objects to be pasted together.}
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
#
# @keyword IO
#*/#########################################################################
setMethodS3("getComment", "RspLanguage", function(this, ...) {
  s <- paste(..., collapse="\n", sep="");
  # By default, no output!
  "";
})


#########################################################################/**
# @RdocMethod escape
#
# @title "Escapes a string specifically for a given RSP response language"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{R objects to be pasted together.}
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
#
# @keyword IO
#*/#########################################################################
setMethodS3("escape", "RspLanguage", function(this, ...) {
  paste(..., collapse="\n", sep="");
})


#########################################################################/**
# @RdocMethod getVerbatim
# @alias getVerbatim.HtmlRspLanguage
#
# @title "Gets a verbatim string specific for a given RSP response language"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{R objects to be pasted together.}
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
#
# @keyword IO
#*/#########################################################################
setMethodS3("getVerbatim", "RspLanguage", function(this, ...) {
  escape(this, ...);
})



##############################################################################
# HISTORY:
# 2006-07-04
# o Added more Rdoc help.
# 2005-08-01
# o Made into an Object class.
# o Added Rdoc comments.
# 2005-07-29
# o Created.
##############################################################################
