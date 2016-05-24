###########################################################################/**
# @RdocClass RspSourceCodeFactory
#
# @title "The RspSourceCodeFactory class"
#
# \description{
#  @classhierarchy
#
#  An RspSourceCodeFactory is language-specific engine that knows how to translate
#  individual @see "RspExpression":s into source code of a specific
#  programming language.
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
#
# @keyword internal
#*/###########################################################################
setConstructorS3("RspSourceCodeFactory", function(language=NA, ...) {
  language <- Arguments$getCharacter(language);
  extend(language, "RspSourceCodeFactory");
})


#########################################################################/**
# @RdocMethod getLanguage
#
# @title "Gets the language"
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
#  Returns an @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getLanguage", "RspSourceCodeFactory", function(this, ...) {
  as.character(this);
})


#########################################################################/**
# @RdocMethod makeSourceCode
#
# @title "Makes a RspSourceCode object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{code}{A @character @vector of code strings.}
#   \item{...}{Arguments passed to the language-specific
#      @see "RspSourceCode" constructor, e.g.
#      \code{type} and \code{metadata}.}
# }
#
# \value{
#  Returns a @see "RspSourceCode" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("makeSourceCode", "RspSourceCodeFactory", function(this, code, ...) {
  lang <- getLanguage(this);
  className <- sprintf("Rsp%sSourceCode", capitalize(lang));
  ns <- getNamespace("R.rsp");
  clazz <- Class$forName(className, envir=ns);
  code <- clazz(code, ...);

  # Get source code header, body, and footer.
  code <- getCompleteCode(this, code, ...);
  code <- c(code$header, code$body, code$footer);

  # Made code object
  code <- clazz(code, ...);

  code;
}, protected=TRUE)




#########################################################################/**
# @RdocMethod exprToCode
# @alias exprToCode.RspRSourceCodeFactory
# @alias exprToCode.RspShSourceCodeFactory
#
# @title "Translates an RspExpression into source code"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{expr}{An @see "RspExpression".}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("exprToCode", "RspSourceCodeFactory", abstract=TRUE);




#########################################################################/**
# @RdocMethod getCompleteCode
# @alias getCompleteCode.RspRSourceCodeFactory
#
# @title "Gets the source code header, body, and footer"
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
setMethodS3("getCompleteCode", "RspSourceCodeFactory", function(this, object, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'object':
  object <- Arguments$getInstanceOf(object, "RspSourceCode");
  lang <- getLanguage(this);
  className <- sprintf("Rsp%sSourceCode", capitalize(lang));
  object <- Arguments$getInstanceOf(object, className);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create header and footer code
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Default header and footer
  header <- '';
  footer <- '';


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Merge all code
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  list(header=header, body=object, footer=footer);
}, protected=TRUE) # getCompleteCode()



#########################################################################/**
# @RdocMethod toSourceCode
#
# @title "Translates an RSP document to source code"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{expr}{An @see "RspDocument" that has been preprocessed
#               and flattened.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns the generated source code as a @see "RspSourceCode" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("toSourceCode", "RspSourceCodeFactory", function(object, doc, ...) {
  # Argument 'doc':
  doc <- Arguments$getInstanceOf(doc, "RspDocument");

  if (length(doc) == 0L) {
    code <- makeSourceCode(object, "", ..., type=getType(doc), metadata=getMetadata(doc, local=TRUE));
    return(code);
  }

  # Assert that the RspDocument 'doc' contains no RspDocument:s
  if (any(sapply(doc, FUN=inherits, "RspDocument"))) {
    throw(sprintf("%s argument 'doc' contains other RspDocuments, which indicates that it has not been flattened.", class(doc)[1L]));
  }

  # Assert that the RspDocument 'doc' contains no RspDirective:s
  if (any(sapply(doc, FUN=inherits, "RspDirective"))) {
    throw(sprintf("%s argument 'doc' contains RSP preprocessing directives, which indicates that it has not been preprocessed.", class(doc)[1L]));
  }

  # Assert that 'doc' contains only RspText:s and RspExpression:s
  nok <- sapply(doc, FUN=function(expr) {
    if (inherits(expr, "RspText") || inherits(expr, "RspExpression")) {
      NA;
    } else {
      class(expr);
    }
  });
  nok <- nok[!is.na(nok)];
  nok <- unique(nok);
  if (length(nok) > 0L) {
    throw(sprintf("%s argument 'doc' contains RSP preprocessing directives, which indicates that it has not been preprocessed: %s", class(doc)[1L], hpaste(nok)));
  }

  # Unescape RspText
  isText <- sapply(doc, FUN=inherits, "RspText");
  doc[isText] <- lapply(doc[isText], FUN=function(expr) {
    RspText(getContent(expr, unescape=TRUE));
  });

  # Coerce all RspConstruct:s to source code
  code <- vector("list", length=length(doc));
  for (kk in seq_along(doc)) {
    code[[kk]] <- exprToCode(object, doc[[kk]], index=kk);
  }
  code <- unlist(code, use.names=FALSE);

  code <- makeSourceCode(object, code, ..., type=getType(doc), metadata=getMetadata(doc, local=TRUE));

  code;
}) # toSourceCode()


##############################################################################
# HISTORY:
# 2013-02-14
# o Added a default getCompleteCode() for RspSourceCodeFactory.
# 2013-02-13
# o ROBUSTNESS: Now toSourceCode() for RspDocument asserts that the document
#   has been flattened and preprocessed.
# 2013-02-11
# o Added Rdoc help.
# 2013-02-10
# o Created.
##############################################################################
