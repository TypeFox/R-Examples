###########################################################################/**
# @RdocClass RspDocument
#
# @title "The RspDocument class"
#
# \description{
#  @classhierarchy
#
#  An RspDocument represents a @list of @see "RspConstruct":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{expressions}{A @list of @see "RspConstruct":s and
#      @see "RspDocument":s.}
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
setConstructorS3("RspDocument", function(expressions=list(), ...) {
  # Argument 'expressions':
  if (!is.list(expressions)) {
    throw("Argument 'expressions' is not a list: ", mode(expressions)[1L]);
  }

  extend(RspObject(expressions, ...), "RspDocument");
})


setMethodS3("print", "RspDocument", function(x, ...) {
  s <- sprintf("%s:", class(x)[1L]);
  s <- c(s, sprintf("Source: %s", getSource(x)));
  s <- c(s, sprintf("Total number of RSP constructs: %d", length(x)));
  if (length(x) > 0L) {
    types <- sapply(x, FUN=function(x) class(x)[1L]);
    if (length(types) > 0L) {
      tbl <- table(types);
      for (kk in seq_along(tbl)) {
        s <- c(s, sprintf("Number of %s(s): %d", names(tbl)[kk], tbl[kk]));
      }
    }
  }
  s <- c(s, sprintf("Content type: %s", getType(x)));
  md <- getMetadata(x, local=FALSE);
  for (key in names(md)) {
    s <- c(s, sprintf("Metadata '%s': '%s'", key, md[[key]]));
  }
  s <- paste(s, collapse="\n");
  cat(s, "\n", sep="");
}, protected=TRUE)




#########################################################################/**
# @RdocMethod getType
#
# @title "Gets the type of the RspDocument"
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
setMethodS3("getType", "RspDocument", function(object, default=NA, as=c("text", "IMT"), ...) {
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
# @title "Gets the source reference of an RSP document"
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
setMethodS3("getSource", "RspDocument", function(object, ...) {
  getAttribute(object, "source", default=NA_character_);
}, protected=TRUE, createGeneric=FALSE)



#########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path to the source reference of an RSP string"
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
setMethodS3("getPath", "RspDocument", function(object, ...) {
  pathname <- getSource(object, ...);
  if (is.na(pathname)) {
    path <- getwd();
  } else {
    path <- getParent(pathname);
  }
  path;
}, protected=TRUE, createGeneric=FALSE)



#########################################################################/**
# @RdocMethod dropEmptyText
#
# @title "Drops all empty RSP text constructs"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "RspDocument".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("dropEmptyText", "RspDocument", function(object, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Nothing to do?
  if (length(object) == 0L) return(object);

  verbose && enter(verbose, "Dropping empty RSP text constructs");

  isEmptyText <- sapply(object, FUN=function(expr) {
    (inherits(expr, "RspText") && (nchar(getContent(expr)) == 0L));
  })
  idxs <- which(isEmptyText);
  n <- length(idxs);
  verbose && cat(verbose, "Number of empty RSP texts: ", n);

  # Anything to drop?
  if (n > 0L) {
    # If dropping everything, at least keep one empty RspText
    # so there will some output
    if (n == length(object)) {
      idxs <- idxs[-n];
    }
    object <- object[-idxs];
  }

  verbose && exit(verbose);

  object;
}) # dropEmptyText()



#########################################################################/**
# @RdocMethod trimNonText
#
# @title "Trims all non-text RSP constructs"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "RspDocument".
# }
#
# \details{
#   For this method to work properly, the RspDocument should not contain
#   any @see "RspUnparsedDirective":s or @see "RspUnparsedExpression":s,
#   i.e. they should already have been parsed.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("trimNonText", "RspDocument", function(object, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  tailString <- function(s, n=10L) {
    len <- nchar(s);
    s <- substring(s, first=max(1L, len-n+1, n));
    s <- gsub("\n", "\\n", s, fixed=TRUE);
    s <- gsub("\r", "\\r", s, fixed=TRUE);
    s;
  } # tailString()

  headString <- function(s, n=10L) {
    s <- substring(s, first=1L, last=n);
    s <- gsub("\n", "\\n", s, fixed=TRUE);
    s <- gsub("\r", "\\r", s, fixed=TRUE);
    s;
  } # headString()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Trimming non-text RSP constructs");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (1) Drop empty text and merge neighboring texts
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop empty RSP texts
  object <- dropEmptyText(object);

  # Merge neighboring RSP texts
  object <- mergeTexts(object);

  isText <- sapply(object, FUN=inherits, "RspText");
  idxsText <- unname(which(isText));
  idxsNonText <- unname(which(!isText));
  idxsSilentNonText <- idxsNonText[!sapply(object[idxsNonText], FUN=getInclude)];
  verbose && cat(verbose, "Number of text RSP constructs: ", length(idxsText));
  verbose && cat(verbose, "Number of non-text RSP constructs: ", length(idxsNonText));
  verbose && cat(verbose, "Number of \"silent\" non-text RSP constructs: ", length(idxsSilentNonText));

  # Nothing todo?
  if (length(idxsNonText) == 0L) {
    verbose && exit(verbose);
    return(object);
  }
  if (length(idxsSilentNonText) == 0L) {
    verbose && exit(verbose);
    return(object);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (2) Drop "empty" RSP text inbetween (non-text) RSP constructs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Dropping 'empty' RSP text inbetween other RSP constructs");
  idxsInbetweenText <- idxsText[1L < idxsText & idxsText < length(object)];
  if (length(idxsInbetweenText) > 0L) {
    for (idx in idxsInbetweenText) {
      # (a) Does the preceeding non-text RSP construct include content?
      item <- object[[idx-1L]];
      if (getInclude(item)) {
         # ... then don't do anything.
         next;
      }

      # (b) Otherwise...
      expr <- object[[idx]];
  ##    verbose && enter(verbose, sprintf("RSP inbetween text #%d ('%s') of %d", idx, class(expr)[1L], length(idxsInbetweenText)));

      text <- getContent(expr);
      # Is text a single line break?
      # (with optional whitespace before and after)?
      isSingleLineBreak <- (regexpr("^[ \t]*(\n|\r|\r\n)[ \t]*$", text) != -1L);
      if (isSingleLineBreak) {
        object[[idx]] <- NA;
      }

    ##    verbose && exit(verbose);
    } # for (idx ...)

    # Cleanup
    excl <- which(sapply(object, FUN=identical, NA));
    if (length(excl) > 0L) {
      object <- object[-excl];

      verbose && cat(verbose, "Number of 'empty' RSP text dropped: ", length(excl));

      isText <- sapply(object, FUN=inherits, "RspText");
      idxsText <- unname(which(isText));
      idxsNonText <- unname(which(!isText));
      idxsSilentNonText <- idxsNonText[!sapply(object[idxsNonText], FUN=getInclude)];
      verbose && cat(verbose, "Number of text RSP constructs: ", length(idxsText));
      verbose && cat(verbose, "Number of non-text RSP constructs: ", length(idxsNonText));
      verbose && cat(verbose, "Number of \"silent\" non-text RSP constructs: ", length(idxsSilentNonText));
    } else {
      verbose && cat(verbose, "No 'empty' RSP text found.");
    }
  } else {
    verbose && cat(verbose, "No inbetween RSP text. Skipping.");
  }
  idxsInbetweenText <- NULL; # Not needed anymore

  verbose && exit(verbose);

  # Sanity checks
  stopifnot(all(idxsText <= length(object)));
  stopifnot(all(idxsNonText <= length(object)));
  stopifnot(all(idxsSilentNonText <= length(object)));
  stopifnot(inherits(object, "RspDocument"));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (3) Drop "empty" line break after non-text RSP constructs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  idxsTextLTrimmed <- NULL;
  idxsTextRTrimmed <- NULL;

  for (kk in seq_along(idxsNonText)) {
    idx <- idxsNonText[kk];
    expr <- object[[idx]];
    verbose && enter(verbose, sprintf("Trimming non-text RSP construct #%d ('%s') of %d", kk, class(expr)[1L], length(idxsNonText)));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (a) Is the RSP construct on its own line?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (i) Find preceeding RSP text
    idxTextL <- idxsText[idxsText < idx];
    if (length(idxTextL) == 0L) {
      textL <- NULL;
      emptyL <- TRUE;
    } else {
      idxTextL <- idxTextL[length(idxTextL)];
      exprL <- object[[idxTextL]];
      textL <- getContent(exprL);
      verbose && printf(verbose, "The text to the left is: '%s'\n", textL);
      emptyL <- (regexpr("\n[ \t\v]*$", textL) != -1L);
    }

    # Not on an empty line?
    if (!emptyL) {
      # We know that this RSP text is non-empty and to the left,
      # so we don't need to consider it again.
      idxsText <- setdiff(idxsText, idxTextL);

      verbose && printf(verbose, "The text to the left is non-empty: '[...]%s'\n", tailString(textL));
      verbose && exit(verbose);
      next;
    }

    # (ii) Find succeeding RSP text
    idxTextR <- idxsText[idxsText > idx];
    if (length(idxTextR) == 0L) {
      textR <- NULL;
      emptyR <- TRUE;
    } else {
      idxTextR <- idxTextR[1L];
      if (idxTextR == idx + 1L) {
        exprR <- object[[idxTextR]];
        textR <- getContent(exprR);
        emptyR <- (regexpr("^[ \t\v]*\n", textR) != -1L);
      } else {
        textR <- NULL;
        emptyR <- TRUE;
      }
    }

    # Not on an empty line?
    if (!emptyR) {
      verbose && printf(verbose, "The text to the right is non-empty: '%s'\n", headString(textR));
      verbose && exit(verbose);
      next;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (b) Now we are working with an non-text RSP construct
    #     that is on an line by itself
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (getInclude(expr)) {
      verbose && printf(verbose, "RSP construct is on its own line but itself includes content (just as RSP text does).\n");
      verbose && exit(verbose);
      next;
    }

    verbose && printf(verbose, "RSP construct is on its own line.\n");


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (c) Trim whitespace and trailing newline
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (i) Trim white space (excluding newline) to the left
    #     (this white space is on the same line as the RSP construct)
    if (!is.null(textL)) {
      textL2 <- gsub("[ \t\v]*$", "", textL);
      if (textL2 != textL) {
        verbose && printf(verbose, "Trimmed %d white-space characters to the left: '%s' -> '%s'\n", nchar(textL)-nchar(textL2), tailString(textL), tailString(textL2));
        exprL2 <- RspText(textL2);
        object[[idxTextL]] <- exprL2;

        # Prevent this RSP text from being trimmed again
        idxsTextLTrimmed <- c(idxsTextLTrimmed, idxTextL);
      }
    }

    # (ii) Trim white space (including newline) to the right
    #     (this white space is on the same line as the RSP construct)
    if (!is.null(textR)) {
      if (regexpr("^[ \t\v]*\n", textR) != -1L) {
        textR2 <- gsub("^[ \t\v]*", "", textR);
        if (textR2 != textR) {
          verbose && printf(verbose, "Trimmed %d white-space characters to the right: '%s' -> '%s'\n", nchar(textR)-nchar(textR2), headString(textR), headString(textR2));
        }

        # Postspone dropping the newline until processing?
        specs <- getSuffixSpecs(expr);
        if (!is.null(specs)) {
          verbose && printf(verbose, "Postponing newline trimming due to suffix specifications: '%s'\n", specs);
        } else {
          textR3 <- gsub("^\n", "", textR2);
          if (textR3 != textR2) {
            verbose && printf(verbose, "Dropped newline to the right: '%s' -> '%s'\n", headString(textR2), headString(textR3));
            textR2 <- textR3;
          }
        }

        exprR2 <- RspText(textR2);
        object[[idxTextR]] <- exprR2;

        # Prevent this RSP text from being trimmed again
        idxsTextRTrimmed <- c(idxsTextRTrimmed, idxTextR);
      }
    }

    verbose && exit(verbose);
  } # for (kk ...)

  verbose && exit(verbose);

  object;
}) # trimNonText()


#########################################################################/**
# @RdocMethod trim
#
# @title "Trims each of the RSP constructs"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{envir}{The @environment where the RSP document is evaluated.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the trimmed @see "RspDocument".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("trim", "RspDocument", function(object, ..., verbose=FALSE) {
  doc <- object;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Trimming RSP text based on surrounding RSP constructs");

  # Identify RSP-only lines by looking at the preceeding
  # and succeeding text parts of each RSP part

  # All RSP text constructs
  isText <- sapply(object, FUN=inherits, "RspText");
  idxs <- which(isText);
  verbose && cat(verbose, "Number of RSP texts: ", length(idxs));

  # Nothing todo?
  if (length(idxs) == 0L) {
    verbose && exit(verbose);
    return(object);
  }

  # Extract RSP texts as plain text
  docT <- unlist(doc[idxs], use.names=FALSE);
  verbose && cat(verbose, "RSP texts as plain text: ");
  verbose && str(verbose, docT);

  # This code assumes that the first and the last part in 'doc'
  # is always a "text" part.
  stopifnot(idxs[1L] == 1L);
#  stopifnot(idxs[length(idxs)] == length(doc));

  # Find text parts that ends with a new line
  endsWithNewline <- (regexpr("(\n|\r|\r\n)[ \t\v]*$", docT[-length(docT)]) != -1L);
  endsWithNewline <- which(endsWithNewline);
  verbose && cat(verbose, "Number of RSP texts ending with an empty line: ", length(endsWithNewline));

  # Don't trim the last RSP text if it is the second last RSP construct
  endsWithNewline <- setdiff(endsWithNewline, length(doc)-1L);

  # Total count of RSP texts trimmed
  count <- 0L;

  # Any candidates?
  if (length(endsWithNewline) > 0L) {
    # Check the following text part
    nextT <- endsWithNewline + 1L;
    docTT <- docT[nextT];

    # Among those, which starts with an empty line?
    startsWithNewline <- (regexpr("^[ \t\v]*(\n|\r|\r\n)", docTT) != -1L);
    startsWithNewline <- nextT[startsWithNewline];
    count <- length(startsWithNewline);
    verbose && cat(verbose, "Number of those RSP texts starting with an empty line: ", count);

    # Any remaining candidates?
    if (count > 0L) {
      # Trim matching text blocks
      endsWithNewline <- startsWithNewline - 1L;

      # Trim to the right (excluding new line because it belongs to text)
      docT[endsWithNewline] <- sub("[ \t\v]*$", "", docT[endsWithNewline]);

      # Trim to the left (drop also any new line because it then
      # belongs to preceeding RSP construct)
      docT[startsWithNewline] <- sub("^[ \t\v]*(\n|\r|\r\n)", "", docT[startsWithNewline]);

      for (kk in seq_along(docT)) {
        value <- RspText(docT[kk]);
        doc[[idxs[kk]]] <- value;
      }
    }
  } # if (length(endsWithNewline) > 0L)

  verbose && cat(verbose, "Number of RSP texts trimmed: ", count);

  verbose && exit(verbose);

  doc;
}, protected=TRUE, createGeneric=FALSE) # trim()



#########################################################################/**
# @RdocMethod mergeTexts
#
# @title "Merge neighboring 'text' elements"
#
# \description{
#  @get "title" by pasting them together.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns an @see "RspDocument" with equal or fever number of elements.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("mergeTexts", "RspDocument", function(object, trim=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Nothing to do?
  if (length(object) <= 1L) return(object);

  # All RSP text constructs
  isText <- sapply(object, FUN=inherits, "RspText");
  idxs <- which(isText);

  # Nothing todo?
  if (length(idxs) == 0L) return(object);

  verbose && enter(verbose, "Merging RSP texts");

  # Locate neighboring RSP text constructs
  while (length(nidxs <- which(diff(idxs) == 1L)) > 0L) {
    idx <- idxs[nidxs[1L]];
    # Merge (idx,idx+1)
    texts <- object[c(idx,idx+1L)];
    text <- paste(texts, collapse="");
    class(text) <- class(texts[[1L]]);
    object[[idx]] <- text;

    # Drop
    object <- object[-(idx+1L)];
    isText <- isText[-(idx+1L)];
    idxs <- which(isText);
  }

  if (trim) {
    verbose && enter(verbose, "Trimming RSP texts");
    object <- trim(object, verbose=verbose);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);


  object;
}, protected=TRUE) # mergeTexts()



#########################################################################/**
# @RdocMethod flatten
#
# @title "Flattens an RspDocument"
#
# \description{
#  @get "title" by expanding and inserting the @list of
#  @see "RspConstruct"s for any @see "RspDocument".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an @see "RspDocument" that contains only @see "RspConstruct":s
#   (and no @see "RspDocument").
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("flatten", "RspDocument", function(object, ..., verbose=FALSE) {
  # Merge neighboring RspText objects
  object <- mergeTexts(object);

  # Nothing to do?
  if (length(object) == 0L) return(object);

  # Nothing todo?
  idxs <- which(sapply(object, FUN=inherits, "RspDocument"));
  if (length(idxs) == 0L) return(object);

  res <- list();

  keys <- names(object);
  for (kk in seq_along(object)) {
    key <- keys[kk];
    expr <- object[[kk]];
    if (inherits(expr, "RspDocument")) {
      expr <- flatten(expr, ..., verbose=verbose);
    } else {
      expr <- list(expr);
      names(expr) <- key;
    }
    res <- append(res, expr);
  } # for (kk ...)

  class(res) <- class(object);

  # Preserve attributes
  res <- setAttributes(res, getAttributes(object));

  # RSP text cleanup
  object <- dropEmptyText(object);
  object <- mergeTexts(object);

  res;
}, protected=TRUE) # flatten()



#########################################################################/**
# @RdocMethod "["
# @aliasmethod "[<-"
#
# @title "Subsets an RspDocument"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{i}{Indices of the RSP elements to extract.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "RspDocument".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("[", "RspDocument", function(x, i) {
  # Preserve the class and other attributes
  res <- .subset(x, i);
  class(res) <- class(x);
  # Preserve attributes
  res <- setAttributes(res, getAttributes(x));
  res;
}, protected=TRUE)

setMethodS3("[<-", "RspDocument", function(x, i, value) {
  # Preserve the class and other attributes
  res <- unclass(x);
  res[i] <- unclass(value);
  class(res) <- class(x);
  # Preserve attributes
  res <- setAttributes(res, getAttributes(x));
  res;
}, protected=TRUE)



#########################################################################/**
# @RdocMethod "subset"
#
# @title "Subsets an RspDocument"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{subset}{An @expression used for subsetting.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "RspDocument".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("subset", "RspDocument", function(x, subset, ...) {
  # To please R CMD check
  doc <- x;

  if (missing(subset)) {
  } else {
    expr <- substitute(subset);
    env <- new.env();
    env$types <- env$names <- names(doc);
    subset <- eval(expr, envir=env, enclos=parent.frame());
    doc <- doc[subset];
  }

  doc;
}, protected=TRUE) # subset()



setMethodS3("asRspString", "RspDocument", function(doc, ...) {
##  isText <- (names(doc) == "text");
##  if (!all(isText)) {
##    throw("Currently it is not possible to coerce an RspDocument to an RspString if it contains elements of other types than 'text': ", hpaste(unique(names(doc))));
##  }

  text <- lapply(doc, FUN=asRspString);
  text <- unlist(text, use.names=FALSE);
  text <- paste(text, collapse="");
  res <- RspString(text, attrs=getAttributes(doc));
  res;
}, protected=TRUE) # asRspString()




setMethodS3("parseIfElseDirectives", "RspDocument", function(object, firstIdx=1L, ..., verbose=FALSE) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'object' & 'firstIdx':
    idx <- firstIdx;
    ifdirective <- object[[idx]];
    if (!inherits(ifdirective, "RspIfDirective")) {
      throw(RspPreprocessingException("First RSP construct is not an RSP 'if' directive", item=ifdirective));
    }

    # Already done?
    value <- getAttribute(ifdirective, ".TRUE");
    if (!is.null(value)) {
      verbose && cat(verbose, "Already parsed. Skipping.");
      verbose && exit(verbose);
      return(object);
    }

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    title <- as.character(asRspString(ifdirective));
    verbose && enter(verbose, "Extracting the 'TRUE' and 'FALSE' statements of ", title);

    verbose && printf(verbose, "RSP 'if-then-else' directive (#%d): %s\n", idx, asRspString(ifdirective));

    idx <- idx + 1L;

    docT <- docF <- list();

    # Build TRUE statement # (find else or endif)
    verbose && enter(verbose, "Collecting 'TRUE' statements for ", title);
    endFound <- FALSE;
    while (idx <= length(object)) {
      item <- object[[idx]];
      if (verbose) itemStr <- gsub("\n", "\\\\n", asRspString(item));

      if (inherits(item, "RspEndifDirective")) {
        verbose && printf(verbose, "Detected ENDIF (#%d: %s)\n", idx, itemStr);
        endFound <- TRUE;
        idx <- idx + 1L;
        break;
      }

      if (inherits(item, "RspElseDirective")) {
        verbose && printf(verbose, "Detected ELSE (#%d: %s)\n", idx, itemStr);
        idx <- idx + 1L;
        break;
      }

      if (inherits(item, "RspIfDirective")) {
        verbose && enter(verbose, sprintf("Detected nested IF (#%d: %s)", idx, itemStr));
        item <- parseIfElseDirectives(object, firstIdx=idx, verbose=verbose);
        if (verbose) {
          for (what in c(".TRUE", ".FALSE")) {
            printf(verbose, "%s statement: {\n", gsub(".", "", what, fixed=TRUE));
            value <- getAttribute(item, what);
            if (!is.null(value)) {
              cat(verbose, asRspString(value));
            }
            printf(verbose, "}\n");
          }
        }

        # Consume indices
        idxs <- getAttribute(item, ".idxs");
        verbose && printf(verbose, "Item range #%d-#%d\n", min(idxs), max(idxs));
        idx <- max(idxs);
        verbose && exit(verbose);
      } else {
        verbose && printf(verbose, "Adding item #%d: '%s'\n", idx, itemStr);
      }

      docT <- c(docT, list(item));

      idx <- idx + 1L;
    } # while()
    verbose && exit(verbose);


    # Build FALSE statement? (find endif)
    if (!endFound) {
      verbose && enter(verbose, "Collecting 'FALSE' statement for ", title);

      while (idx <= length(object)) {
        item <- object[[idx]];
        if (verbose) itemStr <- gsub("\n", "\\\\n", asRspString(item));

        if (inherits(item, "RspEndifDirective")) {
          verbose && printf(verbose, "Detected ENDIF (#%d: %s)\n", idx, itemStr);
          endFound <- TRUE;
          idx <- idx + 1L;
          break;
        }

        if (inherits(item, "RspElseDirective")) {
          throw(RspPreprocessingException(sprintf("Syntax error. Stray RSP 'else' directive (#%d)", idx), item=item));
        }

        if (inherits(item, "RspIfDirective")) {
          verbose && enter(verbose, sprintf("Detected nested IF (#%d: %s)", idx, itemStr));
          item <- parseIfElseDirectives(object, firstIdx=idx, verbose=verbose);
          if (verbose) {
            for (what in c(".TRUE", ".FALSE")) {
              printf(verbose, "%s statement: {\n", gsub(".", "", what, fixed=TRUE));
              value <- getAttribute(item, what);
              if (!is.null(value)) {
                cat(verbose, asRspString(value));
              }
              printf(verbose, "}\n");
            }
          }

          # Consume indices
          idxs <- getAttribute(item, ".idxs");
          verbose && printf(verbose, "Item range #%d-#%d\n", min(idxs), max(idxs));
          idx <- max(idxs);

          verbose && exit(verbose);
        } else {
          verbose && printf(verbose, "Adding item #%d: '%s'\n", idx, itemStr);
        }

        docF <- c(docF, list(item));

        idx <- idx + 1L;
      } # while()

      verbose && exit(verbose);
    }


    if (!endFound) {
      throw(RspPreprocessingException(sprintf("Syntax error. Unclosed RSP 'IF' directive (#%d)", idx), item=ifdirective));
    }

    verbose && printf(verbose, "Consumed items #%d-#%d\n", firstIdx, idx);

    res <- ifdirective;
    attr(res, ".idxs") <- firstIdx:(idx-1L);

    if (length(docT) > 0L) {
      attr(res, ".TRUE") <- RspDocument(docT, attrs=getAttributes(object));
    }
    if (length(docF) > 0L) {
      attr(res, ".FALSE") <- RspDocument(docF, attrs=getAttributes(object));
    }

    verbose && exit(verbose);

    res;
}, protected=TRUE) # parseIfElseDirectives()





setMethodS3("parseCutNPasteDirectives", "RspDocument", function(object, firstIdx=1L, ..., verbose=FALSE) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'object' & 'firstIdx':
    idx <- firstIdx;
    directive <- object[[idx]];

    if (!inherits(directive, "RspCutDirective")) {
      throw(RspPreprocessingException("First RSP construct is not an RSP 'cut' directive", item=directive));
    }

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    title <- as.character(asRspString(directive));
    verbose && enter(verbose, "Extracting the statements of ", title);

    verbose && printf(verbose, "RSP '%s' directive (#%d): %s\n", directive, idx, asRspString(directive));

    idx <- idx + 1L;

    content <- list();

    # Build cut statement # (find endcut)
    verbose && enter(verbose, "Collecting statements for ", title);
    endFound <- FALSE;
    while (idx <= length(object)) {
      item <- object[[idx]];
      if (verbose) itemStr <- gsub("\n", "\\\\n", asRspString(item));

      if (inherits(item, "RspEndcutDirective")) {
        verbose && printf(verbose, "Detected END%s (#%d: %s)\n", toupper(directive), idx, itemStr);
        endFound <- TRUE;
        idx <- idx + 1L;
        break;
      }

      if (inherits(item, "RspCutDirective")) {
        verbose && enter(verbose, sprintf("Detected nested %s (#%d: %s)", toupper(directive), idx, itemStr));
        throw("Nested CUT'N'PASTE directives are not yet supported!");
        item <- parseCutNPasteDirectives(object, firstIdx=idx, verbose=verbose);

        # Consume indices
        idxs <- getAttribute(item, ".idxs");
        verbose && printf(verbose, "Item range #%d-#%d\n", min(idxs), max(idxs));
        idx <- max(idxs);
        verbose && exit(verbose);
      } else {
        verbose && printf(verbose, "Adding item #%d: '%s'\n", idx, itemStr);
      }

      content <- c(content, list(item));

      idx <- idx + 1L;
    } # while()
    verbose && exit(verbose);

    if (!endFound) {
      throw(RspPreprocessingException(sprintf("Syntax error. Unclosed RSP '%s' directive (#%d)", toupper(directive), idx), item=directive));
    }

    verbose && printf(verbose, "Consumed items #%d-#%d\n", firstIdx, idx);

    res <- directive;
    attr(res, ".content") <- content;
    attr(res, ".idxs") <- firstIdx:(idx-1L);

    verbose && exit(verbose);

    res;
}, protected=TRUE) # parseCutNPasteDirectives()



#########################################################################/**
# @RdocMethod preprocess
# @aliasmethod parseIfElseDirectives
# @alias parseIfElseDirectives
# @aliasmethod parseCutNPasteDirectives
# @alias parseCutNPasteDirectives
#
# @title "Processes all RSP preprocessing directives"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{recursive}{If @TRUE, any @see "RspDocument"s introduced via
#      preprocessing directives are recursively parsed and preprocessed
#      as well.}
#   \item{flatten}{If @TRUE, any @see "RspDocument" introduced is
#      replaced (inserted and expanded) by its @list of
#      @see "RspConstruct"s.}
#   \item{envir}{The @environment where the preprocessing is evaluated.}
#   \item{clipboard}{An @environment hold cut'n'paste directives during
#      preprocessing.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "RspDocument".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("preprocess", "RspDocument", function(object, recursive=TRUE, flatten=TRUE, envir=parent.frame(), clipboard=new.env(), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  wrapText <- function(text, wrap=NULL) {
    if (is.null(wrap)) return(text);
    text <- paste(text, collapse="\n");
    text <- gsub("(\r|\r\n)", "\n", text);
    text <- unlist(strsplit(text, split="\n", fixed=TRUE), use.names=FALSE);
    text <- lapply(text, FUN=function(line) {
      first <- seq(from=1L, to=nchar(line), by=wrap);
      last <- first + wrap - 1L;
      substring(line, first=first, last=last);
    });
    text <- unlist(text, use.names=FALSE);
    text <- paste(text, collapse="\n");
    text;
  } # wrapText()


  suffixSpecToCounts <- function(spec, default=1L, specOrg=spec, ...) {
    if (is.null(spec)) {
      count <- 0L;
    } else if (spec == "") {
      count <- default;
    } else if (spec == "*") {
      count <- Inf;
    } else {
      count <- as.numeric(spec);
      if (is.na(count)) {
        if (!identical(spec, specOrg)) {
          spec <- specOrg;
        }
        throw(RspPreprocessingException(sprintf("Invalid/unknown count specifier ('%s') in RSP comment (#%d)", spec, idx)));
      }
    }
    count;
  } # suffixSpecToCounts()

  getFileT <- function(expr, path=".", ..., index=NA, verbose=FALSE) {
    file <- getFile(expr);
    # Sanity check
    stopifnot(!is.null(file));

    verbose && cat(verbose, "Attribute 'file': ", file);

    # URL?
    if (isUrl(file)) {
      verbose && cat(verbose, "URL: ", file);
      fh <- url(file);
      return(fh);
    }

    if (isAbsolutePath(file)) {
      throw(RspPreprocessingException(sprintf("Attribute 'file' specifies an absolute pathname ('%s'). Only relative pathnames are allowed.", file), item=expr));
    }

    if (!is.null(path)) {
      verbose && cat(verbose, "Path: ", path);
      file <- file.path(path, file);
    }
    verbose && cat(verbose, "File: ", file);
    # Sanity check
    stopifnot(!is.null(file));

    if (isUrl(file) || FALSE) {
    } else {
      tryCatch({
        file <- Arguments$getReadablePathname(file);
      }, error = function(ex) {
        throw(RspPreprocessingException(sprintf("File not found (%s), because '%s'", file, gsub("Pathname not found: ", "", ex$message)), item=expr));
      });
    }

    ext <- tolower(file_ext(file));;
    attr(file, "ext") <- ext;

    verbose && cat(verbose, "File: ", file);

    file;
  } # getFileT()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'envir':
  stopifnot(!is.null(envir));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Preprocessing RSP document");
  verbose && cat(verbose, "Number of RSP constructs: ", length(object));

  path <- getPath(object);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (0) Restructure according cut'n'paste
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Parsing cut'n'paste statements");

  idx <- 1L;
  while (idx <= length(object)) {
    item <- object[[idx]];
    if (inherits(item, "RspCutDirective")) {
      cut <- parseCutNPasteDirectives(object, firstIdx=idx, verbose=verbose);
      name <- attr(cut, "name");
      content <- attr(cut, ".content");
      assign(name, content, envir=clipboard, inherits=FALSE);
      # Not needed anymore
      name <- content <- NULL;

      # RSP expressions to be dropped
      idxs <- getAttribute(cut, ".idxs");
      if (item == "copy") {
        # If a copy directive, then keep the content.
        idxs <- range(idxs);
      }

      # Drop
      for (ii in idxs) {
        object[[ii]] <- RspVoid();
      }
    }
    idx <- idx + 1L;
  } # for (idx ...)
  item <- NULL; # Not needed anymore

  verbose && exit(verbose);

  # Assert that all 'cut' statements are consumed
  isCut <- sapply(object, FUN=inherits, "RspCutDirective");
  stopifnot(!any(isCut));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (1) Restructure according to IF-ELSE-THEN directives
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!isTRUE(getAttribute(object, ".ifElseParsed"))) {
    verbose && enter(verbose, "Parsing if-else-then statements");

    items <- list();
    idx <- 1L;
    while (idx <= length(object)) {
      item <- object[[idx]];
      if (inherits(item, "RspIfDirective")) {
        item <- parseIfElseDirectives(object, firstIdx=idx, verbose=verbose);
        idx <- max(getAttribute(item, ".idxs"));
      }
      items <- c(items, list(item));
      idx <- idx + 1L;
    } # for (idx ...)
    item <- NULL; # Not needed anymore

    # Assert that all ELSE and ENDIF directives are gone
    isElse <- sapply(items, FUN=inherits, "RspElseDirective");
    stopifnot(!any(isElse));
    isEndif <- sapply(items, FUN=inherits, "RspEndifDirective");
    stopifnot(!any(isEndif));

    res <- object[c()];
    res[seq_along(items)] <- items;
    object <- res;

    verbose && exit(verbose);
  }


  # Number of empty lines to drop from RSP texts
  nbrOfEmptyTextLinesToDropNext <- 0L;



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (2) Process directives
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Process directives");

  for (idx in seq_along(object)) {
    item <- object[[idx]];
    verbose && enter(verbose, sprintf("RSP construct #%d ('%s') of %d", idx, class(item)[1L], length(object)));

    verbose && cat(verbose, asRspString(item));

    # Number of empty lines to drop from RSP texts
    nbrOfEmptyTextLinesToDrop <- 0L;
    if (nbrOfEmptyTextLinesToDropNext != 0L) {
      nbrOfEmptyTextLinesToDrop <- nbrOfEmptyTextLinesToDropNext;
      nbrOfEmptyTextLinesToDropNext <- 0L;
    }

    # Get the suffix specifications
    spec <- getSuffixSpecs(item);
    if (is.null(spec)) {
      verbose && cat(verbose, "Suffix specifications: <none>");
    } else {
      verbose && printf(verbose, "Suffix specifications: '%s'\n", spec);

      # Don't drop line breaks?
      if (spec == "+") {
        nbrOfEmptyTextLinesToDropNext <- 0L;
      } else if (spec == "-") {
        nbrOfEmptyTextLinesToDropNext <- 1L;
      } else if (regexpr("-\\[(.*)\\]", spec) != -1L) {
        spec <- gsub("-\\[(.*)\\]", "\\1", spec);
        # Expand specifications
        specT <- gstring(spec, envir=envir);
        if (specT != spec) {
          verbose && printf(verbose, "Expanded suffix specifications: '%s'\n", specT);
        }

        # Trim following RSP 'text' construct according to suffix specs?
        nbrOfEmptyTextLinesToDropNext <- suffixSpecToCounts(specT, specOrg=spec);
      } else {
        throw(sprintf("Unknown suffix specification: '%s'", spec));
      }
      verbose && printf(verbose, "Max number of empty lines to drop in next RSP text: %g\n", nbrOfEmptyTextLinesToDropNext);

      # Reset suffix specifications
      attr(item, "suffixSpecs") <- NULL;
      object[[idx]] <- item;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RSP void
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspVoid")) {
      # Drop void RSP items
      object[[idx]] <- NA;
      verbose && exit(verbose);
      next;
    } # RspVoid


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RSP comments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspComment")) {
      # Drop comment
      ## Should we keep this around for improved/better
      ## trimming of newlines? /HB 2014-09-02
      object[[idx]] <- NA;
      verbose && exit(verbose);
      next;
    } # RspComment



    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Keep RSP code expression as is
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspCode")) {
      verbose && exit(verbose);
      next;
    } # RspCode


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Keep RSP text as is, unless empty lines should be dropped
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspText")) {
      # Drop empty lines?
      if (nbrOfEmptyTextLinesToDrop != 0L) {
        content <- getContent(item);

        count <- nbrOfEmptyTextLinesToDrop;

        # Drop all but 'count' empty rows
        if (count < 0) {
          verbose && cat(verbose, "Number of empty lines to drop from the end: ", -count);
          # Count max number of empty rows
          patternR <- "([ \t\v]*(\n|\r|\r\n))*";
          posT <- regexpr(patternR, content);
          if (posT == 1L) {
            nT <- attr(posT, "match.length");
            bfrT <- substring(content, first=1L, last=nT);
            bfrT <- gsub("[ \t\v]*", "", bfrT);
            bfrT <- gsub("\r\n", "\n", bfrT);
            max <- nchar(bfrT);
          } else {
            max <- 0L;
          }

          count <- max + count;
          if (count < 0) count <- 0;
        }

        verbose && cat(verbose, "Number of empty lines to drop: ", count);

        # Drop lines?
        if (count != 0) {
          if (count == 1) {
            patternC <- "?";
          } else if (is.infinite(count)) {
            patternC <- "*";
          } else if (count > 1) {
            patternC <- sprintf("{0,%d}", count);
          }

          # Row pattern
          patternR <- sprintf("([ \t\v]*(\n|\r|\r\n))%s", patternC);

          # Drop empty lines
          content <- sub(patternR, "", content);

          if (nchar(content) > 0L) {
            # Update RspText object
            item2 <- content;
            class(item2) <- class(item);
            object[[idx]] <- item2;
            item2 <- NULL; # Not needed anymore
          } else {
            # ...or drop it if empty
            object[[idx]] <- NA;
          }
        }
      } # if (nbrOfEmptyTextLinesToDrop != 0L)

      verbose && exit(verbose);
      next;
    } # RspText & RspCode


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Support GString-style attribute values for all RSP directives.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspDirective")) {
      attrs <- getAttributes(item);
      for (key in names(attrs)) {
        value <- attrs[[key]];
        value <- gstring(value, envir=envir);
        attr(item, key) <- value;
      }
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Paste
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspPasteDirective")) {
      attrs <- getAttributes(item);
      name <- attrs$name;

      doc <- get(name, envir=clipboard, inherits=FALSE);
      doc <- RspDocument(doc);

      verbose && enter(verbose, "Recursively preprocessing pasted RSP document");
      doc <- preprocess(doc, recursive=TRUE, flatten=flatten, envir=envir, clipboard=clipboard, ..., verbose=verbose);
      metaChild <- getMetadata(doc);
      if (length(metaChild) > 0L) {
        object <- setMetadata(object, metaChild);
      }
      metaChild <- NULL;
      verbose && exit(verbose);

      # Paste content
      object[[idx]] <- doc;
      verbose && exit(verbose);
      next;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RspMetaDirective => ...
    # Setters:
    # <@meta name="<name>" content="<content>"%>
    # <@meta <name>="<content>"%>
    # <@meta <name>="<content>" default="<content>"%>
    # <@meta content="<expr>" lang="<language>"%>
    # Getters:
    # <@meta name="<name>"%>
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspMetaDirective")) {
      attrs <- getNameContentDefaultAttributes(item, doc=object);
      name <- attrs$name;
      content <- attrs$content;

      res <- NA;

      if (!is.null(name) && !is.null(content)) {
        # <@meta name="<name>" content="<content>"%>
        object <- setMetadata(object, name=name, value=content);
      } else if (is.null(name) && !is.null(content)) {
        # <@meta content="<expr>" lang="<language>"%>
        lang <- getAttribute(item, "language");
        if (is.null(lang)) {
          throw(RspPreprocessingException("Attribute 'language' must be specified when parsing metadata from 'content'", item=item));
        }
        if (lang == "R-vignette") {
          metadata <- .parseRVignetteMetadata(content);
        } else {
          throw(RspPreprocessingException(sprintf("Unknown 'language' ('%s')", lang), item=item));
        }
        object <- setMetadata(object, metadata);
      } else if (!is.null(name) && is.null(content)) {
        # <@meta name="<name>"%>
        default <- attrs$default;
        content <- getMetadata(object, name=name, default=default, local=TRUE);
        if (is.null(content)) {
          throw(RspPreprocessingException(sprintf("No such metadata variable ('%s')", name), item=item));
        }
        res <- RspText(content, attrs=getAttributes(object));
      }

      # Drop/insert RSP result
      object[[idx]] <- res;

      verbose && exit(verbose);
      next;
    } # RspMetaDirective


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RspVariableDirective => ...
    # Setters:
    # <@string name="<name>" content="<content>"%>
    # <@string name="<name>" content="<content>" default="<default>"%>
    # <@string <name>="<content>"%>
    # <@string <name>="<content>" default="<default>"%>
    # Getters:
    # <@string name="<name>"%>
    # <@string name="<name>" default="<content>"%>
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspVariableDirective")) {
      attrs <- getNameContentDefaultAttributes(item, doc=object);

      name <- attrs$name;
      if (is.null(name)) {
        throw(RspPreprocessingException("Missing attribute 'name'", item=item));
      }
      value <- attrs$value;

      if (!is.null(value)) {
        # Coerce value
        if (inherits(item, "RspStringDirective")) {
          value <- as.character(value);
        } else if (inherits(item, "RspLogicalDirective")) {
          value <- as.logical(value);
        } else if (inherits(item, "RspIntegerDirective")) {
          value <- as.integer(value);
        } else if (inherits(item, "RspNumericDirective")) {
          value <- as.numeric(value);
        }
      }
      res <- NA;
      if (!is.null(name) && !is.null(value)) {
        # <@string name="<name>" content="<content>"%>
        assign(name, value, envir=envir);
      } else if (!is.null(name) && is.null(value)) {
        # <@string name="<name>"%>
        if (exists(name, envir=envir, inherits=FALSE)) {
          value <- get(name, envir=envir, inherits=FALSE);
        } else {
          value <- attrs$default;
          # <@string name="<name>" default="<content>"%>?
          if (is.null(value)) {
            throw(RspPreprocessingException(sprintf("No such variable ('%s')", name), item=item));
          }
        }
        # Coerce value
        if (inherits(item, "RspStringDirective")) {
          value <- as.character(value);
        } else if (inherits(item, "RspLogicalDirective")) {
          value <- as.logical(value);
        } else if (inherits(item, "RspIntegerDirective")) {
          value <- as.integer(value);
        } else if (inherits(item, "RspNumericDirective")) {
          value <- as.numeric(value);
        }
        res <- RspText(value, attrs=getAttributes(object));
      }

      # Drop/insert RSP result
      object[[idx]] <- res;

      verbose && exit(verbose);
      next;
    } # RspVariableDirective


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RspEvalDirective => ...
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspEvalDirective")) {
      file <- getFile(item);
      content <- getContent(item);
      language <- getAttribute(item, "language", default=NA_character_);
      if (is.null(content) && is.null(file)) {
        throw(RspPreprocessingException("Either attribute 'file' or 'content' must be given", item=item));
      }

      if (!is.null(file)) {
        file <- getFileT(item, path=getPath(object), index=idx, verbose=verbose);
        content <- .readText(file);
      }

      verbose && print(verbose, getAttributes(item));

      if (language == "R") {
        # Parse
        tryCatch({
          expr <- base::parse(text=content);
        }, error = function(ex) {
          throw(sprintf("Failed to parse RSP '%s' directive (%s): %s", item[1L], asRspString(item), ex$message));
        })

        # Evaluate
        tryCatch({
          value <- eval(expr, envir=envir);
        }, error = function(ex) {
          throw(sprintf("Failed to process RSP '%s' directive (%s): %s", item[1L], asRspString(item), ex$message));
        })

        # Drop RSP construct
        object[[idx]] <- NA;

        verbose && exit(verbose);
        next;
      } # if (language == "R")

      if (language == "system") {
        # Evaluate code using system()
        tryCatch({
          value <- system(content, intern=TRUE);
        }, error = function(ex) {
          throw(sprintf("Failed to process RSP '%s' directive (%s): %s", item[1L], asRspString(item), ex$message));
        })

        # Drop RSP construct
        object[[idx]] <- NA;

        verbose && exit(verbose);
        next;
      } # if (language == "system")


      if (language == "shell") {
        # Evaluate code using shell()
        tryCatch({
          value <- shell(content, intern=TRUE);
        }, error = function(ex) {
          throw(sprintf("Failed to process RSP '%s' directive (%s): %s", item[1L], asRspString(item), ex$message));
        })

        # Drop RSP construct
        object[[idx]] <- NA;

        verbose && exit(verbose);
        next;
      } # if (language == "shell")


      throw(RspPreprocessingException(sprintf("Cannot evaluate preprocessing code. Unsupported 'language' ('%s')", language), item=item));
    } # RspEvalDirective


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RspIncludeDirective => ...
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspIncludeDirective")) {
      contentType <- getAttribute(item, "type");

      # Backward compatibility
      if (is.null(contentType)) {
        verbatim <- getAttribute(item, "verbatim");
        if (!is.null(verbatim)) {
          warning("Attribute 'verbatim' for RSP 'include' preprocessing directives is deprecated. Use attribute 'type' instead.");
          if (isTRUE(as.logical(verbatim)))
            contentType <- "text/plain";
        }
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (a) Get content types of host and include document
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      hostContentType <- getType(object, default="text/plain");

      content <- getContent(item);
      if (!is.null(content)) {
        file <- getSource(object);

        # The default content-type of the 'content' attribute is always "text/plain".
        if (is.null(contentType)) {
          contentType <- "text/plain";
        }
      } else {
        file <- getFileT(item, path=getPath(object), index=idx, verbose=verbose);

        # Assert that an endless loop of including the same
        # file over and over does not occur.  This is tested
        # by the number of call frames, which is grows with
        # the number of nested files included.
        if (sys.nframe() > 300L) {
          # For now, don't use throw() because it outputs a very
          # long traceback list.
          stop("Too many nested RSP 'include' preprocessing directives. This indicates an endless recursive loop of including the same file over and over. This was detected while trying to include ", sQuote(file), " (file=", sQuote(getFile(item)), " with type='application/x-rsp') in RSP document ", sQuote(getSource(object)), ".");
        }

        content <- .readText(file);

        # The default content type for the 'file' attribute is
        # inferred from the filename extension, iff possible
        if (is.null(contentType)) {
          ext <- attr(file, "ext");
          if (is.null(ext)) {
            throw(RspPreprocessingException(sprintf("Attribute 'type' must be given because it can not be inferred from the 'file' attribute ('%s') which has no filename extension.", file), item=item));
          }
          contentType <- extensionToIMT(ext=ext, default="text/plain");
        }
      }
      content <- paste(content, collapse="\n");

      # Sanity check
      stopifnot(!is.null(contentType));

      # Parse content types
      hostCT <- parseInternetMediaType(hostContentType);
      inclCT <- parseInternetMediaType(contentType);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (b) Wrap content, iff argument 'wrap' is specified
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      content <- wrapText(content, wrap=getWrap(item));


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (c) Escape content from source and host content type
      #     (This is still very shaky and because it is rather
      #      complicated and there are so many cases to support
      #      it may be dropped in the future.  The 'escaping'
      #      between include to host content types should be
      #      considered a hidden feature. /HB 2013-03-12)
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      defaultEscape <- getAttribute(object, "escape", default=FALSE);
      defaultEscape <- isTRUE(as.logical(defaultEscape));
      escape <- defaultEscape;
      if (escape) {
        content <- escapeRspContent(content, srcCT=inclCT, targetCT=hostCT, verbose=verbose);
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (d) Escape any remaining RSP tags (hide from RSP parser)
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (inclCT$contentType != "application/x-rsp") {
        content <- escapeRspTags(content);
      }


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (e) Parse into an RspText or and RspDocument
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (inclCT$contentType == "application/x-rsp") {
        # "Child" RspDocument:s should "inherit" meta data
        # from the "parent" RspDocument. /HB 2013-11-03
        meta <- getMetadata(object);
        rstr <- RspString(content, type=hostContentType, source=file);
        rstr <- setMetadata(rstr, meta);
        rstr <- setMetadata(rstr, name="source", value=file);
        meta <- NULL; # Not needed anymore

        until <- inclCT$args["until"];
        if (is.null(until)) until <- "*";
        verbose && printf(verbose, "Parsing RSP document until '%s'\n", until);

        # Parse RSP string to RSP document
        doc <- parse(rstr, envir=envir, until=until, verbose=verbose);
        verbose && cat(verbose, "Included RSP document:");
        verbose && print(verbose, doc);

        # Update meta data (child to parent)
        metaChild <- getMetadata(doc);
        if (length(metaChild) > 0L) {
          object <- setMetadata(object, metaChild);
        }
        metaChild <- NULL;

        if (recursive && until == "*") {
          verbose && enter(verbose, "Recursively preprocessing included RSP document");
          doc <- preprocess(doc, recursive=TRUE, flatten=flatten, envir=envir, clipboard=clipboard, ..., verbose=verbose);
          metaChild <- getMetadata(doc);
          if (length(metaChild) > 0L) {
            object <- setMetadata(object, metaChild);
          }
          metaChild <- NULL;
          verbose && exit(verbose);
          item <- doc;
        } else {
          content <- asRspString(doc);
          content <- escapeRspTags(content);
          item <- RspText(content, escape=FALSE, type=hostContentType, source=file);
        }
        rstr <- doc <- NULL; # Not needed anymore
      } else {
        item <- RspText(content, escape=FALSE, type=hostContentType, source=file);
      }

      # Replace RSP directive with imported RSP document
      object[[idx]] <- item;

      content <- item <- NULL; # Not needed anymore

      verbose && exit(verbose);
      next;
    } # RspIncludeDirective


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RspPageDirective => ...
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspPageDirective")) {
      # Update host RSP document attributes
      for (name in c("type", "escape", "language")) {
        value <- getAttribute(item, name, default=getAttribute(object, name));
        object <- setAttribute(object, name, value);
      }

      for (name in c("title", "author", "keywords")) {
        if (!hasAttribute(item, name)) next;
        object <- setMetadata(object, name=name,
                              value=getAttribute(item, name));
      }

      # Drop RSP construct
      object[[idx]] <- NA;

      verbose && exit(verbose);
      next;
    } # RspPageDirective


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RspIfDirective => ...
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspIfDirective")) {
      test <- getAttribute(item, "test");
      attrs <- getNameContentDefaultAttributes(item, known=c("test", "negate"), doc=object);
      name <- attrs$name;

      # Special case <%@if name="<logical>"%>
      if (!hasAttribute(item, "test")) {
        if (!exists(name, envir=envir)) {
          throw(RspPreprocessingException(sprintf("Variable (%s) not found", name), item=item));
        }
        value <- get(name, envir=envir);
        if (is.logical(value)) {
          test <- "equal-to";
        } else {
          throw(RspPreprocessingException("Failed to evaluate IF statement, because attribute 'test' is not specified", item=item));
        }
      }

      # Check for existance of variable
      exist <- FALSE;
      for (mode in c("character", "numeric", "integer", "logical")) {
        exist <- exists(name, mode=mode, envir=envir);
        if (exist) {
          value <- get(name, mode=mode, envir=envir);
          break;
        }
      } # for (mode ...)

      if (test == "exists") {
        result <- exist;
      } else {
        if (!exist) {
          throw(RspPreprocessingException(sprintf("Variable (%s) does not exist", name), item=item));
        }

        otherValue <- attrs$value;
        if (is.null(otherValue)) {
          if (is.logical(value)) {
            otherValue <- TRUE;
          } else {
            otherValue <- "";
          }
        }
        storage.mode(otherValue) <- storage.mode(value);

        if (test == "equal-to") {
          result <- isTRUE(all.equal(value, otherValue));
        } else if (test == "not-equal-to") {
          result <- !isTRUE(all.equal(value, otherValue));
        } else if (test == "greater-than") {
          result <- isTRUE(value > otherValue);
        } else if (test == "greater-than-or-equal-to") {
          result <- isTRUE(value >= otherValue);
        } else if (test == "less-than") {
          result <- isTRUE(value < otherValue);
        } else if (test == "less-than-or-equal-to") {
          result <- isTRUE(value <= otherValue);
        } else {
          throw(RspPreprocessingException(sprintf("Unknown test (%s)", test), item=item));
        }
      }

      # Negate test result?
      negate <- as.logical(getAttribute(item, "negate", FALSE));
      if (negate) {
        result <- !result;
      }

      verbose && enter(verbose, sprintf("Inserting %s statements", result));

      # Extract TRUE or FALSE statements?
      if (result) {
        doc <- getAttribute(item, ".TRUE");
      } else {
        doc <- getAttribute(item, ".FALSE");
      }

      # Recursively pre-process these statements
      if (!is.null(doc)) {
        verbose && print(verbose, doc);
        doc <- setAttribute(doc, ".ifElseParsed", TRUE);
        doc <- preprocess(doc, recursive=TRUE, flatten=flatten, envir=envir, clipboard=clipboard, ..., verbose=verbose);
        # Sanity check
        isIf <- sapply(doc, FUN=inherits, "RspIfDirective");
        stopifnot(!any(isIf));
      } else {
        verbose && print(verbose, "<not available>\n");
        doc <- NA;
      }

      verbose && exit(verbose);

      # Drop/insert RSP result
      object[[idx]] <- doc;

      verbose && exit(verbose);
      next;
    } # RspIfDirective


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RspErrorDirective => ...
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspErrorDirective")) {
      content <- getAttribute(item, "content");
      throw(RspPreprocessingException(content, item=item));
    } # RspErrorDirective


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Stray RSP 'unknown' directive?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspUnknownDirective")) {
      throw(RspPreprocessingException(sprintf("Unknown preprocessing directive (#%d)", idx), item=item));
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Unknown RSP directive?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (inherits(item, "RspDirective")) {
      throw(RspPreprocessingException(sprintf("Do not know how to process an RSP '%s' preprocessing directive (#%d)", item[1L], idx), item=item));
    }

    verbose && exit(verbose);
  } # for (idx ...)

  verbose && exit(verbose);

  # Cleanup (remove NAs, e.g. former RSP comments)
  excl <- which(sapply(object, FUN=identical, NA));
  if (length(excl) > 0L) {
    object <- object[-excl];
  }

##   # Sanity check:
##   # Here all objects are expected to be all RSP text items
##   isOk <- sapply(object, FUN=function(x) {
##     inherits(x, "RspText") || inherits(x, "RspComment")
##   });
##   if (!all(isOk)) {
##     classes <- sapply(object[!isOk], FUN=function(x) class(x)[1L]);
##     tbl <- table(classes);
##     msg <- sprintf("%s [n=%d]", names(tbl), tbl);
##     warning("INTERNAL ERROR: Unexpected classes of RSP items after preprocessing: ", paste(msg, collapse=", "));
##   }

  if (flatten) {
    verbose && enter(verbose, "Flatten RSP document");
    object <- flatten(object, verbose=less(verbose, 10));
    verbose && exit(verbose);
  }

  # RSP text cleanup
  object <- dropEmptyText(object);
  object <- mergeTexts(object);

  if (verbose) {
    if (length(object) > 0L) {
      classes <- sapply(object, FUN=function(x) class(x)[1L]);
      if (length(classes) > 0L) {
        tbl <- table(classes);
        msg <- sprintf("%s [n=%d]", names(tbl), tbl);
        printf(verbose, "Returning RSP document with %d RSP constructs: %s\n", length(object), paste(msg, collapse=", "));
      }
    }
  }


  verbose && exit(verbose);

  object;
}, protected=TRUE) # preprocess()


##############################################################################
# HISTORY:
# 2014-09-03
# o ROBUSTNESS: Now RspIfDirective is more conservative in how it
#   locates RSP variables, i.e. it will only search for RSP variables
#   of certain modes.
# 2014-09-02
# o Clarified some verbose output; useful for troubleshooting.
# 2014-07-02
# o Added support for copy directives.
# o Now the cut'n'paste clipboard is available across files.
# 2014-07-01
# o Added support for cut'n'paste preprocessing directives.
# 2014-05-30
# o RSP directives <%@meta ...%>, <%@string ...%>, ... <%@integer ...%> for
#   getting values gained attribute 'default'.
# 2013-11-03
# o Now "child" RSP documents imported into a "parent" RSP document,
#   sees all meta data of the parent, and any meta data set by the
#   child document are also set in the parent one.  Added a system
#   test for this.
# 2013-10-14
# o BUG FIX: If an RspEvalDirective for language="R" had a parse or an
#   evaluation error, the intended error message was not generated because
#   it in turn would give another error.
# o Grammar correction of a few error messages.
# 2013-06-30
# o Harmonized get- and setMetadata().
# 2013-03-26
# o Now trimNonText() for RspDocument only drops following "empty" text
#   of an RSP construct iff it does not include content itself.
# 2013-03-25
# o BUG FIX: trimNonText() for RspDocument would cause constructs to also
#   drop newlines in text that is following the next construct.
# 2013-03-24
# o Added support for <%@include file="foo.rsp"
#   type="application/x-rsp; until=expressions"%>.
# 2013-03-15
# o Now fully supporting the RSP eval directive.
# 2013-03-13
# o Now preprocess() handles nested if-then-else preprocessing directives.
# o Added protected parseIfElseStatements().
# 2013-03-12
# o Renamed annotations to metadata.
# 2013-03-08
# o Added 'language' attribute to RspIncludeDirective.
# 2013-03-07
# o Added annotation attributes to RspString and RspDocument.
# o Added support for language = "R-vignette" to the RSP 'eval' directive.
#   It parses \Vignette*{} entries to infer RSP title and keywords.
#   The can also be set by the RSP 'page' directive.
# 2013-02-23
# o Added dropEmptyText() and trimNonText() for RspDocument.
# 2013-02-22
# o Added subset() and asRspString() for RspDocument.
# 2013-02-19
# o Now support suffix comment specifications for all RSP expressions.
# o Added mergeTexts() for RspDocument.
# o Added support for <%@ifeq ...%> ... <%@else%> ... <%@endif%> directives.
# 2013-02-14
# o Now RspDocument can include URLs as well.
# 2013-02-13
# o Added getType() for RspDocument.
# o Added support for language:s 'system' and 'shell' for RspEvalDirective.
# o Added print(), preprocess() and flatten() for RspDocument.
# 2013-02-09
# o Created.
##############################################################################
