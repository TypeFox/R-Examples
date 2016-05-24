###########################################################################/**
# @RdocClass RspDirective
#
# @title "The abstract RspDirective class"
#
# \description{
#  @classhierarchy
#
#  An RspDirective is an @see "RspConstruct" that represents an
#  RSP preprocesing directive of format \code{<\%@ ... \%>}.
#  The directive is independent of the underlying programming language.
# }
#
# @synopsis
#
# \arguments{
#   \item{value}{A @character string.}
#   \item{...}{Arguments passed to the constructor of @see "RspConstruct".}
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
setConstructorS3("RspDirective", function(value=character(), ...) {
  extend(RspConstruct(value, ...), "RspDirective");
})


#########################################################################/**
# @RdocMethod "requireAttributes"
#
# @title "Asserts that certain attributes exist"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{condition}{A @character specifying the condition to be tested.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns itself (invisibly).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("requireAttributes", "RspDirective", function(this, names, condition=c("all", "any"), ...) {
  # Argument 'condition':
  condition <- match.arg(condition);

  attrs <- getAttributes(this);
  ok <- is.element(names, names(attrs));

  if (condition == "all") {
    if (!all(ok)) {
      throw(RspPreprocessingException(sprintf("One or more required attributes (%s) are missing", paste(sQuote(names[!ok]), collapse=", ")), item=this));
    }
  } else if (condition == "any") {
    if (!any(ok)) {
      throw(RspPreprocessingException(sprintf("At least one of the required attributes (%s) must be given",  paste(sQuote(names[!ok]), collapse=", ")), item=this));
    }
  }

  invisible(this);
}, protected=TRUE)


setMethodS3("getNameContentDefaultAttributes", "RspDirective", function(item, known=NULL, doc=NULL, ...) {
  name <- getAttribute(item, "name");
  content <- getAttribute(item, "content");
  default <- getAttribute(item, "default");
  file <- getAttribute(item, "file");

  # Was directive given in short format <@<directive> file="<content>">?
  if (is.null(name) && is.null(content) && !is.null(file)) {
    name <- "file";
    content <- file;
    file <- NULL;
  }

  # Was directive given in short format <@<directive> <name>="<content>">?
  if (is.null(name) && is.null(content)) {
    attrs <- getAttributes(item);
    names <- setdiff(names(attrs), c("file", "default", known));
    if (length(names) == 0L) {
      throw(RspPreprocessingException("At least one of attributes 'name' and 'content' must be given", item=item));
    }
    name <- names[1L];
    content <- attrs[[name]];
  }

  # Was directive given with 'file' attribute?
  if (!is.null(file) && !is.null(doc)) {
    path <- getPath(doc);
    if (!is.null(path)) {
      pathname <- file.path(getPath(doc), file);
    } else {
      pathname <- file;
    }
    # Sanity check
    stopifnot(!is.null(pathname));
    content <- .readText(pathname);
  }


  # Use default?
  if (!is.null(content) && (is.na(content) || content == "NA")) {
    value <- default;
  } else {
    value <- content;
  }

  list(name=name, value=value, content=content, file=file, default=default);
}, protected=TRUE) # getNameContentDefaultAttributes()


setMethodS3("asRspString", "RspDirective", function(object, ...) {
  body <- unclass(object);
  attrs <- getAttributes(object);
  if (length(attrs) == 0L) {
    attrs <- "";
  } else {
    attrs <- sprintf('%s="%s"', names(attrs), attrs);
    attrs <- paste(c("", attrs), collapse=" ");
  }

  comment <- getComment(object);
  if (length(comment) == 0L) {
    comment <- "";
  } else {
    comment <- sprintf(" #%s", comment);
  }
  suffixSpecs <- attr(object, "suffixSpecs");
  if (length(suffixSpecs) == 0L) {
    suffixSpecs <- "";
  }
  fmtstr <- "<%%@%s%s%s%s%%>";
  s <- sprintf(fmtstr, body, attrs, comment, suffixSpecs);
  RspString(s);
})



###########################################################################/**
# @RdocClass RspUnparsedDirective
#
# @title "The RspUnparsedDirective class"
#
# \description{
#  @classhierarchy
#
#  An RspUnparsedDirective is an @see RspDirective that still has not
#  been parsed for its class and content.  After @see "parse":ing such
#  an object, the class of this RSP directive will be known.
# }
#
# @synopsis
#
# \arguments{
#   \item{value}{A @character string.}
#   \item{...}{Arguments passed to @see "RspDirective".}
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
setConstructorS3("RspUnparsedDirective", function(value="unparsed", ...) {
  extend(RspDirective(value, ...), "RspUnparsedDirective");
})



#########################################################################/**
# @RdocMethod parse
#
# @title "Parses the unknown RSP directive for its class"
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
#  Returns an @see "RspDirective" of known class.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("parse", "RspUnparsedDirective", function(expr, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  parseAttributes <- function(rspCode, known=mandatory, mandatory=NULL, ...) {
    bfr <- rspCode;

    # Argument 'known':
    known <- unique(union(known, mandatory));

    # Remove all leading white spaces
    pos <- regexpr("^[ \t\n\r]+", bfr);
    len <- attr(pos, "match.length");
    bfr <- substring(bfr, first=len+1L);

    attrs <- list();
    if (nchar(bfr) > 0L) {
      # Add a white space
      bfr <- paste(" ", bfr, sep="");
      while (nchar(bfr) > 0L) {
        # Read all (mandatory) white spaces
        pos <- regexpr("^[ \t\n\r]+", bfr);
        if (pos == -1L) {
          throw(Exception("Error when parsing attributes of RSP preprocessing directive. Expected white space: ", code=sQuote(rspCode)));
        }
        len <- attr(pos, "match.length");
        bfr <- substring(bfr, first=len+1L);

        # Nothing left?
        if (nchar(bfr) == 0L) {
          break;
        }

        # Is the remaining part a comment?
        if (regexpr("^#", bfr) != -1L) {
          # ...then add it as an (R) attribute to 'attrs'.
          comment <- gsub("^#", "", bfr);
          attr(attrs, "comment") <- comment;
          # ...and finish.
          break;
        }

        # Read the attribute name
        pos <- regexpr("^[abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_][abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0-9_]*", bfr);
        if (pos == -1L) {
          throw(Exception("Error when parsing attributes of RSP preprocessing directive. Expected an attribute name: ", code=sQuote(rspCode)));
        }
        len <- attr(pos, "match.length");
        name <- substring(bfr, first=1L, last=len);
        bfr <- substring(bfr, first=len+1L);

        # Read the '=' with optional white spaces around it
        pos <- regexpr("^[ \t\n\r]*=[ \t\n\r]*", bfr);
        if (pos == -1L) {
          throw(Exception("Error when parsing attributes of RSP preprocessing directive. Expected an equal sign: ", code=sQuote(rspCode)));
        }
        len <- attr(pos, "match.length");
        bfr <- substring(bfr, first=len+1L);

        # Work with a raw buffer
        bfrR <- charToRaw(bfr)

        # Read the value with mandatory brackets around it
        # (a) Identify the bracket symbols
        lbracketR <- bfrR[1L]
        lbracket <- rawToChar(lbracketR)
        rbracket <- c("{"="}", "("=")", "["="]", "<"=">")[lbracket];

        # (b) Single brackets or paired ones?
        if (is.na(rbracket)) {
          # (i) Single, e.g. '...', "...", @...@ etc.
          bfrR <- bfrR[-1L];
          wbracket <- 1L;

          # Find first non-escape symbol
          pos <- which(bfrR == lbracketR)

          # Failed to locate a string enclosed in quotation marks
          if (length(pos) == 0L) {
            throw(Exception("Error when parsing attributes of RSP preprocessing directive. Expected an attribute value within quotation marks: ", code=sQuote(rspCode)));
          }

          # An empty value?
          if (pos[1L] == 1L) {
            value <- "";
          } else {
            # Drop escaped brackets
            keep <- (bfrR[pos-1L] != charToRaw("\\"))
            pos <- pos[keep]
            # Failed to locate a string enclosed in quotation marks
            if (length(pos) == 0L) {
              throw(Exception("Error when parsing attributes of RSP preprocessing directive. Expected an attribute value within quotation marks: ", code=sQuote(rspCode)));
            }
            pos <- pos[1L];
            bfrR <- bfrR[1:(pos-1)];
            value <- rawToChar(bfrR);
          }

          # Record brackets
          brackets <- c(lbracket, lbracket);

          # Update buffer
          bfr <- substring(bfr, first=pos+2L);
        } else {
          # (ii) Paired brackets, e.g. {...}, [...], <<...>>

          # Width of left bracket, i.e. how many symbols?
          for (wbracket in seq_len(nchar(bfr))) {
            ch <- substring(bfr, first=wbracket, last=wbracket);
            if (ch != lbracket) {
              wbracket <- wbracket - 1L;
              break;
            }
          }
          bfr <- substring(bfr, first=wbracket+1L);

          # (c) Identify right bracket symbol (escaped for regexpr)
          rbracket <- c("{"="\\}", "("="\\)", "["="\\]", "<"=">",
                        "+"="\\+", "."="\\.", "?"="\\?", "|"="\\|")[lbracket];
          if (is.na(rbracket)) rbracket <- lbracket;

          # Right bracket sequence
          rbrackets <- paste(rep(rbracket, times=wbracket), collapse="");
          # .*? is a non-greedy .* expression
          pattern <- sprintf("^(.*?)([^\\]?)%s", rbrackets);
          pos <- regexpr(pattern, bfr);

          # Failed to locate a string enclosed in brackets
          if (pos == -1L) {
            throw(Exception("Error when parsing attributes of RSP preprocessing directive. Expected a attribute value within brackets: ", code=sQuote(rspCode)));
          }

          # Extract value
          len <- attr(pos, "match.length");
          value <- substring(bfr, first=1L, last=len-wbracket);

          # Record brackets
          lbrackets <- paste(rep(lbracket, times=wbracket), collapse="");
          rbrackets <- gsub("\\\\", "\\", rbrackets);
          brackets <- c(lbrackets, rbrackets);

          # Consume buffer
          bfr <- substring(bfr, first=len+wbracket);
        } # if (is.na(rbracket))

        # Set the name of the value
        names(value) <- name;

        # TODO: Record brackets used
        # ...

        attrs <- c(attrs, value);
      }
    } # if (nchar(bfr) > 0L)

    # Check for duplicated attributes
    if (length(names(attrs)) != length(unique(names(attrs))))
        throw(Exception("Duplicated attributes in RSP preprocessing directive.", code=sQuote(rspCode)));

    # Check for unknown attributes
    if (!is.null(known)) {
      nok <- which(is.na(match(names(attrs), known)));
      if (length(nok) > 0L) {
        nok <- paste("'", names(attrs)[nok], "'", collapse=", ", sep="");
        throw(Exception("Unknown attribute(s) in RSP preprocessing directive: ", nok, code=sQuote(rspCode)));
      }
    }

    # Check for missing mandatory attributes
    if (!is.null(mandatory)) {
      nok <- which(is.na(match(mandatory, names(attrs))));
      if (length(nok) > 0L) {
        nok <- paste("'", mandatory[nok], "'", collapse=", ", sep="");
        throw(Exception("Missing attribute(s) in RSP preprocessing directive: ", nok, code=sQuote(rspCode)));
      }
    }

    # Return parsed attributes.
    attrs;
  } # parseAttributes()


  body <- expr;

  pattern <- "^[ ]*([abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ][abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0-9]*)([ \t\n\r]+(.*))*";

  # Sanity check
  if (regexpr(pattern, body) == -1L) {
    throw("Not an RSP preprocessing directive: ", body);
  }

  # <%@foo attr1="bar" attr2="geek"%> => ...
  directive <- gsub(pattern, "\\1", body);
  directive <- tolower(directive);

  # Parse the attributes
  attrs <- gsub(pattern, "\\2", body);
  attrs <- parseAttributes(attrs, known=NULL);
  comment <- attr(attrs, "comment");

  # Infer the class name
  className <- sprintf("Rsp%sDirective", capitalize(directive));

  # Get constructor
  clazz <- tryCatch({
    ns <- getNamespace("R.rsp");
    Class$forName(className, envir=ns);
  }, error = function(ex) {
    NULL;
  })

  # Instantiate object
  if (!is.null(clazz)) {
    res <- newInstance(clazz, attrs=attrs, comment=comment);
  } else {
    res <- RspUnknownDirective(directive, attrs=attrs);
  }

  # Preserve attributes
  attr(res, "suffixSpecs") <- attr(expr, "suffixSpecs");

  res;
}, createGeneric=FALSE) # parse()


setMethodS3("asRspString", "RspUnparsedDirective", function(object, ...) {
  body <- unclass(object);
  suffixSpecs <- attr(object, "suffixSpecs");
  fmtstr <- "<%%@%s%s%%>";
  s <- sprintf(fmtstr, body, suffixSpecs);
  RspString(s);
})



###########################################################################/**
# @RdocClass RspIncludeDirective
#
# @title "The RspIncludeDirective class"
#
# \description{
#  @classhierarchy
#
#  An RspIncludeDirective is an @see "RspDirective" that causes the
#  RSP parser to include (and parse) an external RSP file.
# }
#
# @synopsis
#
# \arguments{
#   \item{value}{A @character string.}
#   \item{attributes}{A named @list, which must contain either
#      a 'file' or a 'text' element.}
#   \item{...}{Optional arguments passed to the constructor
#              of @see "RspDirective".}
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
setConstructorS3("RspIncludeDirective", function(value="include", ...) {
  this <- extend(RspDirective(value, ...), "RspIncludeDirective");
  if (!missing(value)) {
    requireAttributes(this, names=c("file", "text"), condition="any");
  }
  this;
})



#########################################################################/**
# @RdocMethod getFile
#
# @title "Gets the file attribute"
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
setMethodS3("getFile", "RspIncludeDirective", function(directive, ...) {
  getAttribute(directive, "file");
})

#########################################################################/**
# @RdocMethod getContent
#
# @title "Gets the content of the RSP include directive"
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
setMethodS3("getContent", "RspIncludeDirective", function(directive, ...) {
  getAttribute(directive, "content");
})


#########################################################################/**
# @RdocMethod getVerbatim
#
# @title "Checks if verbatim include should be used or not"
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
#  Returns a @logical.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getVerbatim", "RspIncludeDirective", function(directive, ...) {
  res <- getAttribute(directive, "verbatim", default=FALSE);
  res <- as.logical(res);
  res <- isTRUE(res);
  res;
})


#########################################################################/**
# @RdocMethod getWrap
#
# @title "Get the wrap length"
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
#  Returns an @integer, or @NULL.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getWrap", "RspIncludeDirective", function(directive, ...) {
  res <- getAttribute(directive, "wrap");
  if (!is.null(res)) {
    res <- as.integer(res);
  }
  res;
})




###########################################################################/**
# @RdocClass RspEvalDirective
#
# @title "The RspEvalDirective class"
#
# \description{
#  @classhierarchy
#
#  An RspEvalDirective is an @see "RspDirective" that causes the
#  RSP parser to evaluate a piece of R code (either in a text string
#  or in a file) as it is being parsed.
# }
#
# @synopsis
#
# \arguments{
#   \item{value}{A @character string.}
#   \item{attributes}{A named @list, which must contain a 'file'
#      or a 'text' element.}
#   \item{...}{Optional arguments passed to the constructor
#              of @see "RspDirective".}
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
setConstructorS3("RspEvalDirective", function(value="eval", ...) {
  this <- extend(RspDirective(value, ...), "RspEvalDirective");
  if (!missing(value)) {
    requireAttributes(this, names=c("file", "text"), condition="any");
    lang <- getAttribute(this, default="R");
    this <- setAttribute(this, "language", lang);
  }
  this;
})


#########################################################################/**
# @RdocMethod getFile
#
# @title "Gets the file attribute"
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
setMethodS3("getFile", "RspEvalDirective", function(directive, ...) {
  getAttribute(directive, "file");
})


#########################################################################/**
# @RdocMethod getContent
#
# @title "Gets the content of the RSP eval directive"
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
setMethodS3("getContent", "RspEvalDirective", function(directive, ...) {
  getAttribute(directive, "content");
})


###########################################################################/**
# @RdocClass RspPageDirective
#
# @title "The RspPageDirective class"
#
# \description{
#  @classhierarchy
#
#  An RspPageDirective is an @see "RspDirective" that annotates the
#  content of the RSP document, e.g. the content type.
# }
#
# @synopsis
#
# \arguments{
#   \item{value}{A @character string.}
#   \item{...}{Arguments passed to the constructor of @see "RspDirective".}
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
setConstructorS3("RspPageDirective", function(value="page", ...) {
  extend(RspDirective(value, ...), "RspPageDirective")
})


#########################################################################/**
# @RdocMethod getType
#
# @title "Gets the content type"
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
setMethodS3("getType", "RspPageDirective", function(directive, default=NA, as=c("text", "IMT"), ...) {
  as <- match.arg(as);
  res <- getAttribute(directive, "type", default=as.character(default));
  res <- tolower(res);
  if (as == "IMT" && !is.na(res)) {
    res <- parseInternetMediaType(res);
  }
  res;
})




###########################################################################/**
# @RdocClass RspUnknownDirective
#
# @title "The RspUnknownDirective class"
#
# \description{
#  @classhierarchy
#
#  An RspUnknownDirective is an @see "RspDirective" that is unknown.
# }
#
# @synopsis
#
# \arguments{
#   \item{value}{A @character string.}
#   \item{...}{Arguments passed to the constructor of @see "RspDirective".}
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
setConstructorS3("RspUnknownDirective", function(value="unknown", ...) {
  extend(RspDirective(value, ...), "RspUnknownDirective")
})



###########################################################################/**
# @RdocClass RspErrorDirective
#
# @title "The RspErrorDirective class"
#
# \description{
#  @classhierarchy
#
#  An RspErrorDirective is an @see "RspDirective" that generates an
#  RSP preprocessing error (if processed).
# }
#
# @synopsis
#
# \arguments{
#   \item{value}{A @character string.}
#   \item{...}{Arguments passed to the constructor of @see "RspDirective".}
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
setConstructorS3("RspErrorDirective", function(value="error", ...) {
  extend(RspDirective(value, ...), "RspErrorDirective");
})



##############################################################################
# HISTORY:
# 2014-06-28
# o GENERALIZATION: Now it is possible to use any symbol for enclosing
#   attribute values in RSP directives in addition to the current x='y'
#   and x="y" ones, e.g. x=.y., x=|y|, and so on.  Furthermore, paired
#   brackets may also be used, e.g. x={y}, x=[y], x=<y> and x=(y), and
#   then also in matched replicated, e.g. x={{{y}}}.
# 2014-06-02
# o BUG FIX: getNameContentDefaultAttributes() for RspDirective would
#   return value=NULL if 'content' was an empty string, which was why
#   <%@string empty=''%> would not set 'empty' but instead look it up.
# 2014-05-30
# o Now getNameContentDefaultAttributes() only sets the variable value
#   by the 'default' attribute, iff 'content' is specified.
# 2013-03-26
# o Added getNameContentDefaultAttributes() - used to be a local function
#   of preprocess() for RspDocument.
# 2013-03-25
# o BUG FIX: Forgot to add 'suffixSpecs' to the asRspString() string.
# 2013-03-24
# o BUG FIX: RspEvalDirective() would add erroneous attributes.
# 2013-03-15
# o Added requireAttributes() to RspDirective.
# o Added RSP meta directive.
# 2013-02-23
# o Added asRspString() for RspDirective and RspUnparsedDirective.
# 2013-02-22
# o Added RspUnparsedDirective.
# 2013-02-19
# o Added support for attribute 'text' of RspIncludeDirective:s.
# 2013-02-18
# o Added RspIfeqDirective, RspElseDirective, and RspEndifDirective.
# 2013-02-13
# o Added RspPageDirective.
# o Added 'language' attribute to RspEvalDirective.
# 2013-02-11
# o Added Rdoc help.
# 2013-02-09
# o Created.
##############################################################################
