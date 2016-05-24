###########################################################################/**
# @RdocDefault parseRsp
#
# @title "Parse an RSP code string to an R RSP code string"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{rspCode}{A @character @vector of RSP code to parsed.}
#   \item{rspLanguage}{An @see "RspLanguage" object.}
#   \item{trimRsp}{If @TRUE, white space is trimmed from RSP blocks.}
#   \item{validate}{If @TRUE, the parsed RSP code is validated through the
#     \R parser.}
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
setMethodS3("parseRsp", "default", function(rspCode, rspLanguage=getOption("rspLanguage"), trimRsp=TRUE, validate=TRUE, verbose=FALSE, ...) {
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


  dropRspComments <- function(rspCode, trimRsp=FALSE, ...) {
    pattern <- "<%--.*?--%>";
    if (trimRsp) {
      pattern <- sprintf("%s(|[ \t\v]*(\n|\r|\r\n))", pattern);
    }
    gsub(pattern, "", rspCode, ...);
  } # dropRspComments()



  preprocessRspDirectives <- function(rspCode, ...) {
    rspPattern <- "<%[#@][ ]*([abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ][abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0-9]*)[ ]+(.*?)%>";

    bfr <- rspCode;
    rspCode <- c();
    while (nchar(bfr) > 0) {
      pos <- regexpr(rspPattern, bfr);
      if (pos == -1) {
        break;
      }

      len <- attr(pos, "match.length");
      head <- substring(bfr, first=1, last=pos-1L);
      rspCode <- c(rspCode, head);

      part <- substring(bfr, first=pos, last=pos+len-1L);
      bfr <- substring(bfr, first=pos+len);


      # <%@foo attr1="bar" attr2="geek"%> => ...
      directive <- gsub(rspPattern, "\\1", part);
#      printf("directive: %s\n", directive);

      attrs <- gsub(rspPattern, "\\2", part);
#      printf("attrs: %s\n", attrs);

      if (directive == "insert") {
        attrList <- parseAttributes(attrs, known=NULL);
#        str(attrList);

        path <- attrList$path;
        path <- Arguments$getReadablePath(path, mustExist=TRUE);

        file <- attrList$file;
        pattern <- attrList$pattern;
        if (!is.null(file) && !is.null(pattern)) {
          throw(sprintf("Incorrect RSP directive <%%@insert ...%%>: Only one of attributes 'file' and 'pattern' may be given: <%%@%s %s%%>", directive, attrs));
        }

#str(pattern);
#str(file);
        pathnames <- NULL;
        if (!is.null(pattern)) {
          if (is.null(path))
            path <- ".";
          pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE);
          # Keep only files
          pathnames <- pathnames[sapply(pathnames, FUN=isFile)];
          # Guarantee lexicographic ordering
          pathnames <- sort(pathnames);
        } else if (!is.null(file)) {
          pathname <- Arguments$getReadablePathname(file, path=path, mustExist=TRUE);
#          printf("Pathname: %s\n", pathname);
          pathnames <- pathname;
        } else {
          throw(sprintf("Incomplete RSP directive <%%@insert ...%%>: Either attribute 'file' or 'pattern' must be given: <%%@%s %s%%>", directive, attrs));
        }

        part <- NULL;
        for (pathname in pathnames) {
          bfrT <- readLines(pathname, warn=FALSE);
##          bfrT <- gsub("\\", "\\\\", bfrT, fixed=TRUE);
          bfrT <- paste(bfrT, collapse="\n");
          part <- c(part, bfrT);
        }

        collapse <- attrList$collapse;
        if (is.null(collapse)) collapse <- "\n";
        part <- paste(part, collapse=collapse);
      } # if (directive == ...)

      rspCode <- c(rspCode, part);
    } # while(...)

    rspCode <- c(rspCode, bfr);
    rspCode <- paste(rspCode, collapse="");

    rspCode;
  } # preprocessRspDirectives()



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
##        # Trim trailing white space from RSP tag, if there is
##        # nothing else on the rest of the line?
##        if (trimRsp) {
##          bfr <- sub("^[ \t\v]*(\n|\r|\r\n)", "", bfr);
##        }
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
          throw(Exception("Error when parsing attributes for RSP preprocessing directive. Expected a white space.", code=rspCode));
        len <- attr(pos, "match.length");
        bfr <- substring(bfr, len+1);
        # Read the attribute name
        pos <- regexpr("^[abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ][abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0-9]*", bfr);
        if (pos == -1)
          throw(Exception("Error when parsing attributes for RSP preprocessing directive. Expected an attribute name.", code=rspCode));
        len <- attr(pos, "match.length");
        name <- substring(bfr, 1, len);
        bfr <- substring(bfr, len+1);

        # Read the '=' with optional white spaces around it
        pos <- regexpr("^[ ]*=[ ]*", bfr);
        if (pos == -1)
          throw(Exception("Error when parsing attributes for RSP preprocessing directive. Expected an equal sign.", code=rspCode));
        len <- attr(pos, "match.length");
        bfr <- substring(bfr, len+1);

        # Read the value with mandatory quotation marks around it
        pos <- regexpr("^\"[^\"]*\"", bfr);
        if (pos == -1)
          throw(Exception("Error when parsing attributes for RSP preprocessing directive. Expected a quoted attribute value string.", code=rspCode));
        len <- attr(pos, "match.length");
        value <- substring(bfr, 2, len-1);
        bfr <- substring(bfr, len+1);
        names(value) <- name;
        attrs <- c(attrs, value);
      }
    } # if (nchar(bfr) > 0)

    # Check for duplicated attributes
    if (length(names(attrs)) != length(unique(names(attrs))))
        throw(Exception("Duplicated attributes in RSP preprocessing directive.", code=rspCode));

    # Check for unknown attributes
    if (!is.null(known)) {
      nok <- which(is.na(match(names(attrs), known)));
      if (length(nok) > 0) {
        nok <- paste("'", names(attrs)[nok], "'", collapse=", ", sep="");
        throw(Exception("Unknown attribute(s) in RSP preprocessing directive: ", nok, code=rspCode));
      }
    }

    # Check for missing mandatory attributes
    if (!is.null(mandatory)) {
      nok <- which(is.na(match(mandatory, names(attrs))));
      if (length(nok) > 0) {
        nok <- paste("'", mandatory[nok], "'", collapse=", ", sep="");
        throw(Exception("Missing attribute(s) in RSP preprocessing directive: ", nok, code=rspCode));
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
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'rspCode'
  rspCode <- Arguments$getCharacters(rspCode);

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

  # Argument 'validate'
  validate <- Arguments$getLogical(validate);




  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tidy up RSP code
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Concatenate rsp code strings
  rspCode <- paste(rspCode, collapse="\n");
  rspCode <- paste(rspCode, "\n", sep="");



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (a) Preprocess RSP code
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  lastRspCode <- "";
  while (rspCode != lastRspCode) {
    lastRspCode <- rspCode;

    # Drop RSP comments
    # [1] Example Depot, Greedy and Nongreedy Matching in a Regular
    #     Expression, 2009.
    #     http://www.exampledepot.com/egs/java.util.regex/Greedy.html
    rspCode <- dropRspComments(rspCode, trimRsp=trimRsp);

    # Preprocessing RSP directives should go here, e.g. @insert.
    rspCode <- preprocessRspDirectives(rspCode);
  } # while (...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (b) Parse RSP code
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Split in non-RSP and RSP parts, e.g splitting by '<%...%>'.
  parts <- splitRspTags(rspCode, trimRsp=trimRsp);
  rspCode <- NULL; # Not needed anymore



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (c) Translate RSP code (to R code)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  error <- NULL;
  rCode <- NULL;

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
          code <- sprintf("write(response, \"%s\");\n", value);
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
      codeComment <- unlist(strsplit(rspTag, split="\n", fixed=TRUE), use.names=FALSE);
      codeComment <- paste("# ", codeComment, sep="");
      codeComment <- c(codeComment, "");
      codeComment <- paste(codeComment, collapse="\n");
      rspCode <- trim(part);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # RSP Scripting Elements and Variables
      #
      # <%--[comment]--%>
      #
      # NOTE: With dropRspComments() above, this will never occur here.
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
      # <%# [comment] %>  [MAY BE AMBIGOUS! /HB 2011-03-16]
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##       pattern <- "^#(.*)$";
##       if (regexpr(pattern, part) != -1) {
##         # <%# [comment] %>  => # [comment]
##         comment <- gsub(pattern, "\\1", part);
##         rCode <- c(rCode, comment);
##         next;
##       }

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
        code <- c(codeComment, sprintf("write(response, {%s});\n", value));

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


        # <%@include file="url"%> => add what parseRsp(url) returns.
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
              path <- getwd();
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
            value <- parseRsp(text=lines, path=getParent(file));
#            value <- c(value, "\n# Resets the 'page' object\n");
#            value <- c(value, pageCode);
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
        } # if (directive == ...)

        # <%@directive attr1="foo" attr2="bar"%>
        #     => write(response, directive(attr1="foo", attr2="bar"))
        # TODO: Try to parse here to catch invalid code as soon as possible?
        rspDirective <- directive;
        args <- sprintf("%s=\"%s\"", names(attrs), attrs);
        args <- paste(args, collapse=", ");
        cmd <- sprintf("%s(%s)", rspDirective, args);
        code <- c(codeComment, "write(response, ", cmd, ");\n");

        rCode <- c(rCode, code);
        next;
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # RSP Scripting Elements and Variables
      #
      # <%:: [expressions] %>  - Output the code and evaluate it
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pattern <- "^::(.*)";
      if (regexpr(pattern, part) != -1) {
        expressions <- gsub(pattern, "\\1", part);
        expressions <- trim(expressions);
        expressions <- paste(expressions, "\n", sep="");

        flavor <- c("echo", "chunk")[2];
        if (flavor == "echo") {
          expressions <- strsplit(expressions, split="\n", fixed=TRUE);
          expressions <- unlist(expressions, use.names=FALSE);
          code <- sprintf("write(response, R.utils::withCapture({%s}), collapse='\\n');\n", expressions);
        }

        if (flavor == "chunk") {
          # Sanity check (this code chunk needs to be complete!)
          tryCatch({
            parse(text=expressions);
          }, error = function(ex) {
            throw(ex);
          });
          expressions <- paste(expressions, collapse="");
          # FIXME: The following does not handled "nested" code
          # inside strings, e.g. cat("cat('cat(\''hello\\')')").
          # /HB 2011-03-16
          expressions <- gsub("'", "\\'", expressions, fixed=TRUE);
          code <- sprintf("write(response, sourceWithTrim('%s', echo=TRUE), collapse='\\n');\n", expressions);
        }
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



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (d) Validate R code (via parsing)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (validate) {
    # Parse the R RSP code so that error messages contains line numbers.
    con <- textConnection(rCode);
    on.exit(close(con));

    tryCatch({
      # Parse parsed R code
      rExpr <- parse(con);
    }, error = function(ex) {
      # If an parse error occurs, show tranlated code.
      msg <- ex$message;
      line <- gsub(".*line *([0-9]+).*", "\\1", msg);
      code <- displayCode(code=rCode, highlight=line, pager="none");
      code <- unlist(strsplit(code, split="\n", fixed=TRUE), use.names=FALSE);
      ex$code <- code;
      stop(ex);
    });
  }
  attr(rCode, "validate") <- validate;

  rCode;
})

##############################################################################
# HISTORY:
# 2011-11-14
# o ROBUSTNESS: Now <%=[expr]%> is translated with curly brackets around
#   the expression, i.e. write(response, {[expr]}).  This allows for
#   writing <%= x <- 1; x^2 %> instead of <%={ x <- 1; x^2 }%>.
# 2011-11-07
# o BUG FIX: <@fcn foo="bar"> for "fallback" directives would try to call
#   fcn(foo=bar) instead of fcn(foo="bar").
# o CLEANUP: Now translating to evalCapture() instead of evalWithEcho().
# o BUG FIX: One of the tags would generate invalid evalWithEcho() code.
# 2011-04-12
# o Change the preprocess directives to have format <%#insert ...%>.
# 2011-04-01
# o Added support for <%insert path="<path>" pattern="<pattern>"%>.
# o Added support for <%insert file="<filename>" path="<path>"%>.
# o Added support for <%insert file="<pathname>"%>.
# o Added internal processRspInserts().
# 2011-03-30
# o Now parseRsp() drops RSP comments, i.e. '<%-- {anything} --%>'.
# o Added internal dropRspComments() for dropping '<%-- {anything} --%>'.
# 2011-03-15
# o Added RSP markup <%;{R code}%> for evaluating and echoing code.
# 2011-03-12
# o Now the trimming of RSP handles all newline types, i.e. LF, CR+LF, CR.
# 2011-03-08
# o Added argument 'trimRsp' to parseRsp() for trimming white space
#   surrounding RSP blocks that have preceeding and succeeding white space
#   and that are followed by a newline.
# 2011-02-13
# o BUG FIX: parseRsp() would generate invalid R code/R comments for
#   multiline <%=...%> statements.
# 2009-02-23
# o Now parseRsp() does the validation of the R code.
# o Renamed translateRsp() to parseRsp().
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
