###########################################################################/**
# @RdocClass RspString
#
# @title "The RspString class"
#
# \description{
#  @classhierarchy
#
#  An RspString is a @character @vector with RSP markup.
# }
#
# @synopsis
#
# \arguments{
#   \item{s}{A @character @vector.}
#   \item{...}{Arguments passed to @see "RspObject".}
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
setConstructorS3("RspString", function(s=character(), ...) {
  # Argument 's':
  s <- paste(s, collapse="\n");

  extend(RspObject(s, ...), "RspString");
})


setMethodS3("print", "RspString", function(x, ...) {
  s <- sprintf("%s:", class(x)[1L]);
  s <- c(s, sprintf("Content type: %s", getAttribute(x, "type", NA)));
  s <- c(s, sprintf("Language: %s", getAttribute(x, "language", NA)));
  metadata <- getMetadata(x, local=FALSE);
  if (length(metadata) > 0L) {
    metadata <- unlist(metadata, use.names=TRUE);
    s <- c(s, sprintf("Metadata '%s': %s", names(metadata), metadata));
  } else {
    s <- c(s, "Metadata to available.");
  }
  s <- c(s, sprintf("Number of characters: %s", nchar(x)));
  s <- c(s, sprintf("Number of lines: %s", nbrOfLines(x)));
  ruler <- paste(rep("#", times=getOption("width")-2L), collapse="");
  s <- c(s, ruler, x);
  s <- paste(s, collapse="\n");
  cat(s, "\n", sep="");
}, protected=TRUE)


#########################################################################/**
# @RdocMethod nbrOfLines
#
# @title "Gets the number of lines in an RSP string"
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
#  Returns a non-negative @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("nbrOfLines", "RspString", function(object, ...) {
  length(unlist(strsplit(object, split="\n", fixed=TRUE), use.names=FALSE));
})



#########################################################################/**
# @RdocMethod getType
#
# @title "Gets the type of an RSP string"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{default}{If unknown/not set, the default content type to return.}
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
setMethodS3("getType", "RspString", function(object, default=NA, as=c("text", "IMT"), ...) {
  as <- match.arg(as);
  res <- getAttribute(object, "type", default=as.character(default));
  res <- tolower(res);
  if (as == "IMT" && !is.na(res)) {
    res <- parseInternetMediaType(res);
  }
  res;
}, protected=TRUE)



#########################################################################/**
# @RdocMethod getSource
#
# @title "Gets the source reference of an RSP string"
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
setMethodS3("getSource", "RspString", function(object, ...) {
  getAttribute(object, "source", default=as.character(NA));
}, protected=TRUE, createGeneric=FALSE)




#########################################################################/**
# @RdocMethod parse
#
# @title "Parses the RSP string"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments passed to the RSP parser.}
#   \item{envir}{The @environment where the RSP document is parsed.}
#   \item{parser}{An @see "RspParser".}
# }
#
# \value{
#  Returns a @see "RspDocument" (unless \code{until != "*"} in case it
#  returns a deparsed @see "RspString".)
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("parse", "RspString", function(object, ..., envir=parent.frame(), parser=RspParser()) {
  # Argument 'parser':
  parser <- Arguments$getInstanceOf(parser, "RspParser");

  parse(parser, object, ..., envir=envir);
}, createGeneric=FALSE, protected=TRUE) # parse()



##############################################################################
# HISTORY:
# 2013-03-14
# o Added a print() method for RspStrings.
# 2013-03-09
# o Moved all parsing code to the new RspParser class.
# 2013-03-07
# o Added annotation attributes to RspString and RspDocument.
# 2013-02-13
# o Added getType() for RspString.
# 2013-02-11
# o Added Rdoc help.
# 2013-02-09
# o Created.
##############################################################################
