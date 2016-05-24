###########################################################################/**
# @RdocClass RspShSourceCodeFactory
#
# @title "The RspShSourceCodeFactory class"
#
# \description{
#  @classhierarchy
#
#  An RspShSourceCodeFactory is an @see "RspSourceCodeFactory" for
#  the shell ('sh') script language.
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
setConstructorS3("RspShSourceCodeFactory", function(...) {
  extend(RspSourceCodeFactory("sh"), "RspShSourceCodeFactory");
})



setMethodS3("exprToCode", "RspShSourceCodeFactory", function(object, expr, ..., index=NA) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  escapeRspText <- function(text) {
    text <- deparse(text);
    text <- substring(text, first=2L, last=nchar(text)-1L);
    text;
  } # escapeRspText()

  makeCode <- function(code, echo=FALSE, include=FALSE, ...) {
    code <- unlist(strsplit(code, split="\n", fixed=TRUE), use.names=FALSE);
    codeT <- trim(code);

    n <- length(code);
    codeE <- sapply(code, FUN=escapeRspText);
    codeE <- sprintf("printf \"%s\"", codeE);
    suffixR <- rep(" > /dev/null", times=n);
    codeR <- sprintf("%s%s", codeT, suffixR);
    if (include) {
      # Output the last out
      codeR[n] <- sprintf("printf \"%s\"", code[n]);
    }

    codeS <- matrix(c(codeE, codeR), nrow=2L, byrow=TRUE);
    rownames(codeS) <- c("echo", "include");

    if (echo && !include) {
      code <- codeS[1L,,drop=TRUE];
    } else if (echo && include) {
      code <- codeS;
    } else if (!echo && include) {
      code <- codeS[2L,,drop=TRUE];
    } else if (!echo && !include) {
      code <- codeS[2L,,drop=TRUE];
    }

    code;
  } # makeCode()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'expr':
  reqClasses <- c("RspText", "RspExpression");
  if (!inherits(expr, reqClasses)) {
    throw("Argument 'expr' must be of class RspText or RspExpression: ", class(expr)[1L]);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # RspText => echo "<text>"
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (inherits(expr, "RspText")) {
    text <- getContent(expr);

    code <- NULL;
    while (nchar(text) > 0L) {
      textT <- substring(text, first=1L, last=1024L);
      textT <- escapeRspText(textT);
      codeT <- sprintf("printf \"%s\"", textT);
      code <- c(code, codeT);
      text <- substring(text, first=1025L);
    }
    if (is.null(code)) {
      code <- "printf \"\\n\"";
    }

    return(code);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # RspCodeChunk => ...
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (inherits(expr, "RspCodeChunk")) {
    code <- makeCode(getCode(expr), echo=getEcho(expr), include=getInclude(expr));
    return(code);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # RspCode => <code>
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (inherits(expr, "RspCode")) {
    code <- makeCode(getCode(expr), echo=getEcho(expr), include=FALSE);
    return(code);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # RspComment => [void]
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (inherits(expr, "RspComment")) {
    return("");
  }


  throw(sprintf("Unknown class of RSP expression (#%d): %s", index, class(expr)[1L]));
}, protected=TRUE) # exprToCode()



##############################################################################
# HISTORY:
# 2013-03-14
# o Created.
##############################################################################
