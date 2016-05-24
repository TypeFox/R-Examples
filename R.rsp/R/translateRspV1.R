###########################################################################/**
# @RdocDefault translateRspV1
#
# @title "Translates an RSP file to an R servlet"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{file}{A filename, a URL, or a @connection to be read.
#               Ignored if \code{text} is not @NULL.}
#   \item{text}{If specified, a @character @vector of RSP code to be
#               translated.}
#   \item{path}{A pathname setting the current include path.
#               If \code{file} is a filename and its parent directory
#               is different from this one, \code{path} is added
#               to the beginning of \code{file} before the file is read.}
#   \item{rspLanguage}{An @see "RspLanguage" object.}
#   \item{trimRsp}{If @TRUE, white space is trimmed from RSP blocks.}
#   \item{verbose}{Either a @logical, a @numeric, or a @see "R.utils::Verbose"
#     object specifying how much verbose/debug information is written to
#     standard output. If a Verbose object, how detailed the information is
#     is specified by the threshold level of the object. If a numeric, the
#     value is used to set the threshold of a new Verbose object. If @TRUE,
#     the threshold is set to -1 (minimal). If @FALSE, no output is written.
#     [Currently not used.]
#   }
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string of \R source code.
# }
#
# @author
#
# \seealso{
#   @see "sourceRsp".
# }
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("translateRspV1", "default", function(file="", text=NULL, path=getParent(file), rspLanguage=getOption("rspLanguage"), trimRsp=TRUE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  trimTextParts <- function(parts, ...) {
    ## cat("TRIMMING...\n");
    # Identify RSP-only lines by looking at the preceeding
    # and succeeding text parts of each RSP part

    # This code assumes that the first and the last part in 'parts'
    # is always a "text" part.
    stopifnot(names(parts)[1] == "text");
    stopifnot(names(parts)[length(parts)] == "text");

    # Identify all text parts
    idxs <- which(names(parts) == "text");
    partsT <- unlist(parts[idxs], use.names=FALSE);

    # Find text parts that ends with a new line
    endsWithNewline <- (regexpr("(\n|\r|\r\n)[ \t\v]*$", partsT[-length(partsT)]) != -1);
    endsWithNewline <- which(endsWithNewline);

    # Any candidates?
    if (length(endsWithNewline) > 0) {
      # Check the following text part
      nextT <- endsWithNewline + 1L;
      partsTT <- partsT[nextT];

      # Among those, which starts with a new line?
      startsWithNewline <- (regexpr("^[ \t\v]*(\n|\r|\r\n)", partsTT) != -1);
      startsWithNewline <- nextT[startsWithNewline];

      # Any remaining candidates?
      if (length(startsWithNewline) > 0) {
        # Trim matching text blocks
        endsWithNewline <- startsWithNewline - 1L;

        # Trim to the right (excluding new line because it belongs to text)
        partsT[endsWithNewline] <- sub("[ \t\v]*$", "", partsT[endsWithNewline]);

        # Trim to the left (including new line because it belongs to RSP)
        partsT[startsWithNewline] <- sub("^[ \t\v]*(\n|\r|\r\n)", "", partsT[startsWithNewline]);

        parts[idxs] <- partsT;
      }
    }
    ## cat("TRIMMING...done\n");

    parts;
  } # trimTextParts()


  splitRspTags <- function(..., trimRsp=FALSE) {
    bfr <- paste(..., collapse="\n", sep="");

    START <- 0;
    STOP <- 1;

    parts <- list();
    state <- START;
    while(TRUE) {
      if (state == START) {
        # The start tag may exists *anywhere* in static code
        pos <- regexpr("<%", bfr);
        if (pos == -1)
          break;

        part <- list(text=substring(bfr, first=1, last=pos-1));
        bfr <- substring(bfr, first=pos+2);
        state <- STOP;
      } else if (state == STOP) {
        pos <- indexOfNonQuoted(bfr, "%>");
        if (pos == -1)
          break;

        part <- list(rsp=substring(bfr, first=1, last=pos-1));
        bfr <- substring(bfr, first=pos+2);
        state <- START;
      }

      parts <- c(parts, part);
    } # while(TRUE);

    # Add the rest of the buffer as text
    parts <- c(parts, list(text=bfr));

    if (trimRsp) {
      parts <- trimTextParts(parts);
    }

    parts;
  } # splitRspTags()



  parseAttributes <- function(rspCode, known=mandatory, mandatory=NULL, ...) {
    bfr <- rspCode;

    # Argument 'known':
    known <- unique(union(known, mandatory));

    # Remove all leading white spaces
    pos <- regexpr("^[ \t]+", bfr);
    len <- attr(pos, "match.length");
    bfr <- substring(bfr, len+1);

    attrs <- list();
    if (nchar(bfr) >= 0) {
      # Add a white space
      bfr <- paste(" ", bfr, sep="");
      while (nchar(bfr) > 0) {
        # Read all (mandatory) white spaces
        pos <- regexpr("^[ \t]+", bfr);
        if (pos == -1)
          throw(Exception("Error when parsing attributes. Expected a white space.", code=rspCode));
        len <- attr(pos, "match.length");
        bfr <- substring(bfr, len+1);
        # Read the attribute name
        pos <- regexpr("^[abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ][abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0-9]*", bfr);
        if (pos == -1)
          throw(Exception("Error when parsing attributes. Expected an attribute name.", code=rspCode));
        len <- attr(pos, "match.length");
        name <- substring(bfr, 1, len);
        bfr <- substring(bfr, len+1);

        # Read the '=' with optional white spaces around it
        pos <- regexpr("^[ ]*=[ ]*", bfr);
        if (pos == -1)
          throw(Exception("Error when parsing attributes. Expected an equal sign.", code=rspCode));
        len <- attr(pos, "match.length");
        bfr <- substring(bfr, len+1);

        # Read the value with mandatory quotation marks around it
        pos <- regexpr("^\"[^\"]*\"", bfr);
        if (pos == -1)
          throw(Exception("Error when parsing attributes. Expected a quoted attribute value string.", code=rspCode));
        len <- attr(pos, "match.length");
        value <- substring(bfr, 2, len-1);
        bfr <- substring(bfr, len+1);
        names(value) <- name;
        attrs <- c(attrs, value);
      }
    } # if (nchar(bfr) > 0)

    # Check for duplicated attributes
    if (length(names(attrs)) != length(unique(names(attrs))))
        throw(Exception("Duplicated attributes.", code=rspCode));

    # Check for unknown attributes
    if (!is.null(known)) {
      nok <- which(is.na(match(names(attrs), known)));
      if (length(nok) > 0) {
        nok <- paste("'", names(attrs)[nok], "'", collapse=", ", sep="");
        throw(Exception("Unknown attribute(s): ", nok, code=rspCode));
      }
    }

    # Check for missing mandatory attributes
    if (!is.null(mandatory)) {
      nok <- which(is.na(match(mandatory, names(attrs))));
      if (length(nok) > 0) {
        nok <- paste("'", mandatory[nok], "'", collapse=", ", sep="");
        throw(Exception("Missing attribute(s): ", nok, code=rspCode));
      }
    }

    # Return parsed attributes.
    attrs;
  } # parseAttributes()

  # 2005-08-12, Ana-Catarina 2.8kg, kl. 17.17 lokal tid, 48.5cm

  # Function to escape characters so that they can be included within an
  # R character string, e.g. to put 'size="-1"' becomes "size=\"-1\"".
  ASCII.ESCAPED <- ASCII;
  ASCII.ESCAPED[0] <- "\\\\x";
  ASCII.ESCAPED[1:31] <- sprintf("\\\\%03d", as.integer(intToOct(1:31)));
  ASCII.ESCAPED[1+7:13] <- c("\\\\a", "\\\\b", "\\\\t", "\\\\n",
                                              "\\\\v", "\\\\f", "\\\\r");
  # Using non-standard character turns out to be non-supported in
  # some locales. See HISTORY.
  # MAGIC.STRING <- "\255\001\255\002\255\003";
  MAGIC.STRING <- "THISISMYMAGICSTRINGIGUESSNOBODYELSEWOULDPUTTHESAMEINARSPPAGE";

  escapeRspText <- function(text) {
    # Substitute all '\' with '\\'.
    # Comment: A '\' has to be escaped in the C regexpr() function, that
    #          is '\\', which each in turn is written as "\\" in R.
    text <- gsub("\\\\([^\"])", "\\\\\\\\\\1", text);

    # Substitute all '"' with tempory ASCII sequence 'MAGIC.STRING'
    text <- gsub("([^\\]|^)\"", paste("\\1", MAGIC.STRING, sep=""), text);

    # Substitute all '\"' with '\\"'
    text <- gsub("\\\"", "\\\\\\\\\"", text);

    # Substitute all ASCII sequences 'MAGIC.STRING' with '\\"'
    text <- gsub(MAGIC.STRING, "\\\\\"", text);

    # Escape all non-printable characters
    for (kk in 1+1:31) {
      text <- gsub(ASCII[kk], ASCII.ESCAPED[kk], text);
    }

    # Substitute all '<\%' with '<%'.
    text <- gsub("<\\\\%", "<%", text);

    text;
  } # escapeRspText()

  escapeRspText <- function(text) {
    text <- deparse(text);
    text <- substring(text, 2, nchar(text)-1);
    text;
  } # escapeRspText()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # MAIN
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'text'
  if (!is.null(text) && length(as.character(text)) == 0)
    return("");

  # Argument 'rspLanguage'
  if (is.null(rspLanguage)) {
    rspLanguage <- RspLanguage();
  } else if (is.character(rspLanguage)) {
    rspLanguage <- Arguments$getCharacter(rspLanguage);
    clazz <- paste(capitalize(rspLanguage), "RspLanguage", sep="");
    tryCatch({
      clazz <- Class$forName(clazz);
      rspLanguage <- newInstance(clazz);
    }, error=function(ex) {
      throw("No such 'rspLanguage' (\"", rspLanguage, "\"): ", clazz);
    })
  } else if (!inherits(rspLanguage, "RspLanguage")) {
    throw("Argument 'rspLanguage' is not a RspLanguage object: ",
                                                     class(rspLanguage)[1]);
  }

  # Argument 'path'
  if (is.null(path)) {
    path <- ".";
  } else {
    path <- Arguments$getReadablePathname(path, mustExist=FALSE);
  }

  # Argument 'file'
  pathname <- "";
  if (is.null(text) && is.character(file)) {
    if (file == "") {
      text <- readLines(warn=FALSE);
    } else {
      if (isUrl(file)) {
        pathname <- file;
      } else {
        if (!identical(getParent(file), path)) {
          pathname <- filePath(path, file);
        } else {
          pathname <- file;
        }

        if (!isFile(pathname))
          throw("Cannot translate RSP file. File not found: ", pathname);
      }
      text <- readLines(pathname, warn=FALSE);
    }
  } else {
    # When does this happen? /HB 2006-07-04
  }

  # Argument 'trimRsp'
  trimRsp <- Arguments$getLogical(trimRsp);


  text <- paste(paste(text, collapse="\n"), "\n", sep="");

  # Split in non-RSP and RSP parts, e.g splitting by '<%...%>'.
  parts <- splitRspTags(text, trimRsp=trimRsp);
  text <- NULL; # Not needed anymore

  error <- NULL;

  # Translate RSP document to R code
  rCode <- paste(
  "#######################################################################\n",
  "# DO NOT EDIT!  DO NOT EDIT!  DO NOT EDIT!  DO NOT EDIT!  DO NOT EDIT! \n",
  "#                                                                      \n",
  "# This R code was translated from RSP by the R.rsp package.            \n",
  "#                                                                      \n",
  "# Details:                                                             \n",
  "# File: ", file, "\n",
  "# Path: ", path, "\n",
  "#######################################################################\n",
  "\n", sep="");

  code <- "# Assert that write() of R.rsp is used below\n";
  rCode <- c(rCode, code);
  code <- "write <- R.rsp:::write;\n";
  rCode <- c(rCode, code);

  code <- "# Sets the public RspPage 'page' object\n";
  rCode <- c(rCode, code);
  pageCode <- paste("page <- RspPage(pathname=\"", pathname, "\");\n", sep="");
  rCode <- c(rCode, pageCode);
  code <- "# Gets the output connection (or filename) for the response [OBSOLETE]\n";
  rCode <- c(rCode, code);
  code <- "out <- getOutput(response);\n";
  rCode <- c(rCode, code);

  types <- names(parts);
  for (kk in seq(length=length(parts))) {
    part <- parts[[kk]];
    type <- types[[kk]];

    if (type == "text") {
      # [text] => write(response, "[escaped text]");
      if (nchar(part) > 0) {
        while (nchar(part) > 0) {
          currPart <- substring(part, 1, 1024);
          value <- escapeRspText(currPart);
          code <- c("write(response, \"", value, "\");\n");
          rCode <- c(rCode, code);
          part <- substring(part, 1025);
        }
      } else {
        code <- part;
      }
      next;
    }

    if (type == "rsp") {
      rspTag <- paste("<%", part, "%>", sep="");
      rspTagE <- gsub("\n", "\\n", rspTag, fixed=TRUE)
      rspTagE <- gsub("\r", "\\r", rspTagE, fixed=TRUE)
      codeComment <- paste("# ", rspTagE, "\n", sep="");
      rspCode <- trim(part);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # RSP Scripting Elements and Variables
      #
      # <%--[comment]--%>
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pattern <- "^--(.*)--$";
      if (regexpr(pattern, part) != -1) {
        # <%--[comment]--%>  => # [comment]
        comment <- gsub(pattern, "\\1", part);
#        rCode <- c(rCode, comment);
        next;
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # RSP Scripting Elements and Variables
      #
      # <%=[expression]%>
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pattern <- "^=(.*)$";
      if (regexpr(pattern, part) != -1) {
        # <%=[expression]%> => write(response, [expression]\n);
        value <- gsub(pattern, "\\1", part);
        value <- trim(value);
        # TODO: Try to parse here to catch invalid code as soon as possible?
        code <- c(codeComment, "write(response, {", value, "});\n");

        rCode <- c(rCode, code);
        next;
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # RSP Directives
      #
      # <%@ directive attr1="foo" attr2="bar" ...%>
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pattern <- "^@[ ]*([abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ][abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0-9]*)[ ]+(.*)$";
      if (regexpr(pattern, part) != -1) {
        # <%@foo attr1="bar" attr2="geek"%> => ...
        directive <- gsub(pattern, "\\1", part);
        attrs <- gsub(pattern, "\\2", part);
        attrs <- parseAttributes(attrs, known=NULL);


        # <%@include file="url"%> => add what translateRsp(url) returns.
        if (directive == "include") {
          file <- attrs[["file"]];
          if (is.null(file))
            throw("Attribute 'file' is missing: ", rspTag);

          verbatim <- identical(as.logical(attrs[["verbatim"]]), TRUE);
          wrap <- attrs[["wrap"]];
          newline <- attrs[["newline"]];

          if (isUrl(file)) {
            fh <- url(file);
            lines <- readLines(fh, warn=FALSE);
          } else {
            if (!isAbsolutePath(file)) {
              file <- filePath(path, file);
              file <- getAbsolutePath(file);
            }

            if (!isFile(file)) {
              throw("Cannot include file. File not found: ", file);
            }
            lines <- readLines(file, warn=FALSE);
          }

          if (verbatim) {
            if (!is.null(wrap)) {
              wrap <- as.integer(wrap);
              lines <- unlist(sapply(lines, FUN=function(line) {
                first <- seq(from=1, to=nchar(line), by=wrap);
                last <- first + wrap - 1;
                substring(line, first, last);
              }), use.names=FALSE)
            }
            value <- getVerbatim(rspLanguage, lines, newline=newline);
            value <- paste("write(response, \"",
                                  escapeRspText(value), "\");\n", sep="");
          } else {
            # Process and include file.
            value <- translateRspV1(text=lines, path=getParent(file));
            value <- c(value, "\n# Resets the 'page' object\n");
            value <- c(value, pageCode);
          }

          code <- c(codeComment, value);

          rCode <- c(rCode, code);
          next;
        }


        # <%@import file="url"%> => import(response, url)
        if (directive == "import") {
          file <- attrs[["file"]];
          if (is.null(file))
            throw("Attribute 'file' is missing: ", rspTag);

          code <- c(codeComment,
                    "import(response, \"", file, "\", path=\"", path, "\");\n");

          rCode <- c(rCode, code);
          next;
        }

        # <%@page ...%> => ...
        if (directive == "page") {
          code <- c();

          import <- attrs[["import"]];
          if (!is.null(import)) {
            packages <- strsplit(import, split=";|,|:")[[1]];
            code <- paste("library(", packages, ");\n", sep="");
            rCode <- c(rCode, code);
          }

          language <- attrs[["language"]];
          if (!is.null(language)) {
            rCode <- c(rCode, code);
          }

          contentType <- attrs[["contentType"]];
          if (!is.null(contentType)) {
            tmp <- strsplit(contentType, split=";", fixed=TRUE)[[1]];
            mime <- tmp[1];
            args <- tmp[-1];
            if (mime == "text/html") {
              rspLanguage <- HtmlRspLanguage();
##            } else if (mime == "text/latex") {
##              rspLanguage <- LaTeXRspLanguage();
##            } else if (mime == "text/plain") {
##              rspLanguage <- TextRspLanguage();
            }
          }

          info <- attrs[["info"]];
          if (!is.null(info)) {
            comment <- paste("# ", info, "\n", sep="");
            value <- getComment(rspLanguage, info);
            code <- c(comment, "write(response, \"", escapeRspText(value), "\");\n");
            rCode <- c(rCode, code);
          }

          next;
        }

        # <%@directive attr1="foo" attr2="bar"%>
        #     => write(response, directive(attr1="foo", attr2="bar"))
        # TODO: Try to parse here to catch invalid code as soon as possible?
        rspDirective <- directive;
        args <-paste(names(attrs), attrs, sep="=");
        args <- paste(args, collapse=", ");
        value <- paste(rspDirective, "(response, ", args, ")", sep="");
        code <- c(codeComment, "write(response, ", value, ");\n");

        rCode <- c(rCode, code);
        next;
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # RSP Scripting Elements and Variables
      #
      # <%: [expressions] %>  - Output the code and evaluate it
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pattern <- "^:(.*)";
      if (regexpr(pattern, part) != -1) {
        expressions <- gsub(pattern, "\\1", part);

        expressions <- paste(trim(expressions), "\n", sep="");
        code <- c("write(response, \"", escapeRspText(expressions), "\");\n",
                                                    trim(expressions), "\n");

        rCode <- c(rCode, code);
        next;
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # RSP Scripting Elements and Variables
      #
      # <% [expressions] %>  - Include [expression]\n
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      code <- paste(rspCode, "\n", sep="");
      rCode <- c(rCode, code);
    } # if (type == "rsp")
  } # for (kk in ...)

  # Paste all code snippets together
  rCode <- paste(rCode, collapse="", sep="");

  attr(rCode, "pathname") <- pathname;

  rCode;
})

##############################################################################
# HISTORY:
# 2014-10-18
# o CLEANUP/ROBUSTNESS: translateRsp() and translateRspV1(), which are
#   both deprecated, no longer assume that write() is exported from R.rsp.
# 2011-11-17
# o Now the generated R script adds 'write <- R.rsp::write' at the
#   beginning, to assure that it is used instead of base::write().
# 2011-03-12
# o Now the trimming of RSP handles all newline types, i.e. LF, CR+LF, CR.
# 2011-03-08
# o Added argument 'trimRsp' to translateRspV1() for trimming white spaces
#   surrounding RSP blocks that have preceeding and succeeding white space
#   and that are followed by a newline.
# 2009-02-23
# o There is a new translateRsp(). The old version is keep for backward
#   compatibility as translateRspV1().
# 2007-04-07
# o Replace regexpr pattern "^[ \]*=[ \]*" with "^[ \]*=[ \]*".
# 2006-07-20
# o BUG FIX: An RSP comment tag would also replicate last text or R code.
# 2006-07-17
# o BUG FIX: translateRsp("\\\n") would convert to "\\n".  Update internal
#   escapeRspText().  Thanks Peter Dalgaard for the suggestions.
# 2006-07-05
# o BUG FIX: If argument 'path' was NULL, translateRsp() gave an error.
# 2006-07-04
# o Now translateRsp() returns attribute 'pathname' too.  Used by sourceRsp().
# o The assigned 'out' object is obsolete.  Instead there should be an
#   RspResponse object.
# 2006-01-14 (Julien Gagneur)
# o BUG FIX: Changed value of variable MAGIC.STRING, the former was not
#   compatible with gsub under some locales. /JG
# 2005-08-15
# o Now all output is written as GString:s; updated the RspResponse class.
# o Now static text '<\%' is outputted as '<%'.
# o Now the 'out' (connection or filename) is available in the servlet code.
# o Replaced tag '<%#' with '<%:'.
# o Added support for page directive 'import'.
# 2005-08-13
# o BUG FIX: Forgot to add newline after translating a scripting element.
# o Updated escapeRspText() too escape all ASCII characters from 1 to 31 (not
#   zero though).  It also escapes the double quote character where needed.
# 2005-08-01
# o Replace importRsp() with import(response, ...).
# o Added Rdoc comments.
# 2005-07-31
# o Recreated again. Before the RSP code was translated to an output document
#   immediately, but now an intermediate R code is created.  That is, when
#   before the process was RSP -> HTML, it is now RSP -> R -> HTML.  This
#   makes it possible to create much richer RSP documents. Specifically, it
#   is possible to write R code statements spanning more than one RSP tag.
# 2005-07-29
# o Recreated from previous RspEngine class in the R.io package, which was
#   first written in May 2002. See source of old R.io package for details.
##############################################################################
