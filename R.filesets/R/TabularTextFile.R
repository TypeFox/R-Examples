###########################################################################/**
# @RdocClass TabularTextFile
#
# @title "The TabularTextFile class"
#
# \description{
#  @classhierarchy
#
#  A TabularTextFile is an object refering to a tabular text file
#  on a file system containing data in a tabular format.
#  Methods for reading all or a subset of the tabular data exist.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "GenericTabularFile".}
#   \item{sep}{A @character specifying the symbol used to separate the
#     cell entries.  If more than one symbol is specified, it will try to
#     select the correct one by peeking into the file.}
#   \item{quote}{A @character specifying the quote symbol used, if any.}
#   \item{fill}{As in @see "utils::read.table".}
#   \item{skip}{As in @see "utils::read.table".}
#   \item{columnNames}{A @logical or a @character @vector. If @TRUE,
#      then column names are inferred from the file.  If a @character
#      @vector, then the column names are given by this argument.}
#   \item{commentChar}{A single @character specifying which symbol
#      should be used for comments, cf. @see "utils::read.table".}
#   \item{.verify, verbose}{(Internal only) If @TRUE, the file is
#      verified while the object is instantiated by the constructor.
#      The verbose argument is passed to the verifier function.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @examples "../incl/TabularTextFile.readDataFrame.Rex"
#
# @author
#
# \seealso{
#   An object of this class is typically part of an @see "TabularTextFileSet".
# }
#*/###########################################################################
setConstructorS3("TabularTextFile", function(..., sep=c("\t", ","), quote="\"", fill=FALSE, skip=0L, columnNames=NA, commentChar="#", .verify=TRUE, verbose=FALSE) {
  # Argument 'columnNames':
  if (is.logical(columnNames)) {
    readColumnNames <- columnNames;
    columnNames <- NULL;
  } else if (is.character(columnNames)) {
    readColumnNames <- FALSE;
  } else {
    throw("Argument 'columnNames' must be either a logical or a character vector: ", class(columnNames)[1]);
  }

  # Argument 'commentChar':
  if (identical(commentChar, "") || identical(commentChar, FALSE)) {
    commentChar <- NULL
  } else if (!is.null(commentChar)) {
    commentChar <- Arguments$getCharacter(commentChar, nchar=c(1,1));
  }


  this <- extend(GenericTabularFile(..., .verify=FALSE), "TabularTextFile",
    "cached:.fileHeader" = NULL,
    "cached:.nbrOfLines" = NULL,
    .columnNameTranslator = NULL,
    sep = sep,
    quote = quote,
    fill = fill,
    skip = skip,
    .commentChar = commentChar,
    .columnNames = columnNames,
    readColumnNames = readColumnNames
  );

  if (.verify) {
    verify(this, ..., verbose=verbose);
    # Clear temporary settings
    this$.fileHeader <- NULL;
  }

  this;
})


setMethodS3("as.character", "TabularTextFile", function(x, ...) {
  this <- x;
  s <- NextMethod("as.character");
  colnames <- getColumnNames(this);
  if (length(colnames) > 0L) {
    columns <- paste("'", colnames, "'", sep="");
    s <- c(s, sprintf("Columns [%d]: %s", length(columns), paste(columns, collapse=", ")));
  } else {
    s <- c(s, sprintf("Columns [NA]: <not reading column names>"));
  }
  s <- c(s, sprintf("Number of text lines: %d", nbrOfLines(this, fast=TRUE)));
  s;
}, protected=TRUE)


setMethodS3("verify", "TabularTextFile", function(this, ..., verbose=FALSE) {
  # Nothing to do?
  pathname <- getPathname(this);
  if (is.null(pathname) || is.na(pathname))
    return(invisible(this));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Validating file contents");

  tryCatch({
    data <- readDataFrame(this, nrow=10L, verbose=verbose);
  }, error = function(ex) {
    throw("File format error of the tabular file ('", getPathname(this), "'): ", ex$message);
  })

  verbose && exit(verbose);

  invisible(this);
}, private=TRUE)


setMethodS3("getCommentChar", "TabularTextFile", function(this, ...) {
  this$.commentChar;
}, protected=TRUE)


setMethodS3("setCommentChar", "TabularTextFile", function(this, ch, ...) {
  if (!is.null(ch)) {
    ch <- Arguments$getCharacter(ch, nchar=c(1,1));
  }
  this$.commentChar <- ch;
  invisible(this);
}, protected=TRUE)


setMethodS3("readColumnNames", "TabularTextFile", function(this, ...) {
  res <- as.logical(this$readColumnNames);

  # No need to infer from header?
  if (!is.na(res)) {
    return(res);
  }

  hdr <- readRawHeader(this, ...);
  namesHdr <- hdr$commentArgs$columnNames;
  if (is.null(namesHdr)) {
    msg <- sprintf("Cannot infer whether the data table has column names or not, because header comment argument 'columnNames' is missing. Will assume there are column names (readColumnNames=TRUE): %s", getPathname(this));
#    warning(msg);
    return(TRUE);
  }

  # If argument 'columnNames' equals the first data row, then yes.
  namesData <- hdr$topRows[[1]];
  res <- all(namesData == namesHdr);

  # Store results(?!?)
#  this$readColumnNames <- res;

  res;
}, protected=TRUE)


###########################################################################/**
# @RdocMethod hasColumnHeader
#
# @title "Checks if there are column names in the header"
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
#   Returns a @logical.
# }
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("hasColumnHeader", "TabularTextFile", function(this, ...) {
  readColumnNames(this);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getDefaultColumnNames
#
# @title "Gets the default column names"
#
# \description{
#  @get "title" by inferring it from the file header.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Optional arguments passed @seemethod "getHeader".}
# }
#
# \value{
#   Returns @character @vector,
#   or @NULL if there are no column names in the file header.
# }
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("getDefaultColumnNames", "TabularTextFile", function(this, ...) {
  # (a) Hardcoded/user-specified column names?
  names <- this$.columnNames;
  if (!is.null(names)) {
    return(names);
  }

  hdr <- getHeader(this, ...);

  # (b) Infer column names from data table
  if (hasColumnHeader(this)) {
    names <- hdr$columns;
    return(names);
  }

  # (c) Infer column names from header argument 'columnNames'?
  useHeaderArgs <- this$.useHeaderArgs;
  if (is.null(useHeaderArgs)) useHeaderArgs <- TRUE;
  if (useHeaderArgs) {
    args <- hdr$commentArgs;
    if (length(args) > 0L) {
      names <- args$columnNames;
      return(names);
    }
  }

  # (c) There are no column names
  NULL;
}, protected=TRUE)


setMethodS3("getDefaultColumnClasses", "TabularTextFile", function(this, ...) {
  ncol <- nbrOfColumns(this);

  hdr <- getHeader(this, ...);

  # (a) Infer column names from header argument 'columnNames'?
  useHeaderArgs <- this$.useHeaderArgs;
  if (is.null(useHeaderArgs)) useHeaderArgs <- TRUE;
  if (useHeaderArgs) {
    args <- hdr$commentArgs;
    if (length(args) > 0L) {
      colClasses <- args$columnClasses;
      if (!is.null(colClasses)) {
        # Sanity check
        stopifnot(length(colClasses) == ncol);
        return(colClasses);
      }
    }
  }

  # (b) Column classes are not specified in the file
  NULL;
}, protected=TRUE)


setMethodS3("getDefaultColumnClassPatterns", "TabularTextFile", function(this, ...) {
  names <- getColumnNames(this, ...);
  colClasses <- getDefaultColumnClasses(this, ...);
  if (is.null(colClasses)) return(NULL);

  names <- sprintf("^%s$", names);
  names(colClasses) <- names;
  colClasses;
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getHeader
#
# @title "Gets the file header"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Passed to internal @seemethod "readRawHeader".}
#   \item{header}{A @logical specifying whether there are column
#    headers or not.}
#   \item{force}{If @TRUE, an already retrieved header will be ignored.}
# }
#
# \value{
#   Returns a named @list.
# }
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("getHeader", "TabularTextFile", function(this, ..., header=TRUE, force=FALSE) {
  hdr <- this$.fileHeader;
  if (force || is.null(hdr) || hasBeenModified(this, update=FALSE)) {
    hdr <- readRawHeader(this, ...);
    if (header && hasColumnHeader(this)) {
      hdr$columns <- hdr$topRows[[1]];
    }
    this$.fileHeader <- hdr;
  }
  hdr;
})




setMethodS3("readRawHeader", "TabularTextFile", function(this, con=NULL, skip=this$skip, sep=this$sep, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'skip':
  if (is.null(skip)) {
    skip <- 0L;
  } else if (is.character(skip)) {
    skip <- Arguments$getRegularExpression(skip);
  } else {
    skip <- Arguments$getInteger(skip, range=c(0,Inf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading tabular file header from ", class(this)[1]);

  # Open a file connection?
  if (is.null(con)) {
    pathname <- getPathname(this);
    verbose && cat(verbose, "Pathname: ", pathname);

    # Open file connection
    con <- file(pathname, open="r");
    on.exit({
      if (!is.null(con)) {
        close(con);
        con <- NULL;
      }
    })
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # How to skip
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  skipMax <- 0L;

  # Skip a fixed number of rows?
  if (is.numeric(skip)) {
    skippedLines <- readLines(con, n=skip);
    skipMax <- skipMax + skip;
  }

  # Skip until regular expression?
  skipUntil <- NULL;
  if (is.character(skip)) {
    skipUntil <- skip;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Header comments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read header comments
  comments <- c();
  ch <- getCommentChar(this);
  ch <- if (is.null(ch)) "#" else ch;
  pattern <- sprintf("^%s", ch);
  ready <- FALSE;
  while (!ready) {
    line <- readLines(con, n=1L);
    isEmpty <- (regexpr("^$", line) != -1L);
    if (!isEmpty) {
      if (!is.null(skipUntil)) {
        if (regexpr(skipUntil, line) != -1L) {
          break;
        }
      }
      isComments <- (regexpr(pattern, line) != -1L);
      if (isComments) {
        comments <- c(comments, line);
        skipMax <- skipMax + 1L;
      } else if (is.null(skipUntil)) {
        break;
      }
    }
  } # while(!ready)

  verbose && cat(verbose, "Header comments:", level=-20);
  verbose && str(verbose, comments, level=-20);


  # Parse header comments
  pattern <- sprintf("^%s[ ]*([^:]+):[ ]*(.*)", ch);
  commentArgs <- grep(pattern, comments, value=TRUE);
  if (length(commentArgs) > 0L) {
    keys <- gsub(pattern, "\\1", commentArgs);
    keys <- trim(keys);
    commentArgs <- gsub(pattern, "\\2", commentArgs);
    commentArgs <- trim(commentArgs);
    commentArgs <- strsplit(commentArgs, split="\t", fixed=TRUE);
    names(commentArgs) <- keys;
  } else {
    commentArgs <- NULL;
  }
  verbose && cat(verbose, "Header comment arguments:", level=-20);
  verbose && str(verbose, commentArgs, level=-20);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Column separator
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Infer column separator from the first line after the header comments?
  if (length(sep) > 1L) {
    verbose && enter(verbose, "Identifying the separator that returns most columns");
    verbose && cat(verbose, "Line:");
    verbose && print(verbose, line);
    verbose && cat(verbose, "Separators:");
    verbose && str(verbose, sep);
    columns <- base::lapply(sep, FUN=function(split) {
      strsplit(line, split=split)[[1]];
    });
    nbrOfColumns <- sapply(columns, FUN=length);
    max <- which.max(nbrOfColumns);
    sep <- sep[max];
    verbose && printf(verbose, "Choosen separator: '%s' (0x%s)\n", sep, charToRaw(sep));
    verbose && exit(verbose);
  }

  lines <- c(line, readLines(con, n=9));
  verbose && print(verbose, line);
  topRows <- strsplit(lines, split=sep);
  topRows <- lapply(topRows, trim);
  verbose && print(verbose, topRows);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove quotes?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  quote <- this$quote;
  if (!is.null(quote)) {
    for (pattern in c(sprintf("^%s", quote), sprintf("%s$", quote))) {
      topRows <- lapply(topRows, FUN=function(row) {
        gsub(pattern, "", row);
      })
    }
  }

  verbose && cat(verbose, "Columns: ", paste(paste("'", topRows, "'", sep=""), collapse=", "), level=-10);

  hdr <- list(
    comments=comments,
    commentArgs=commentArgs,
    sep=sep,
    quote=quote,
    skip=skip,
    skipMax=skipMax,
    topRows=topRows
  );

  verbose && str(verbose, hdr);

  verbose && exit(verbose);

  hdr;
}, protected=TRUE); # readRawHeader()


setMethodS3("getReadArguments", "TabularTextFile", function(this, fileHeader=NULL, colClasses=c("*"=NA, getDefaultColumnClassPatterns(this)), defColClass="NULL", stringsAsFactors=FALSE, skip=this$skip, ..., verbose=FALSE) {

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fileHeader':
  if (is.null(fileHeader)) {
    fileHeader <- getHeader(this, skip=skip, ...);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Building arguments for read.table()");


  # Optional user arguments
  userArgs <- list(stringsAsFactors=stringsAsFactors, ...);
  verbose && cat(verbose, "User arguments:");
  verbose && str(verbose, userArgs);

  # Backward compatibility
  if (is.element("colClassPatterns", names(userArgs))) {
    .Deprecated(msg="Argument 'colClassPatterns' has been renamed to 'colClasses'. Please update your code accordingly.");
    colClasses <- userArgs[["colClassPatterns"]];
    userArgs[["colClassPatterns"]] <- NULL;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Infer column classes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Default column classes
  columns <- getColumnNames(this);
  hasColClassPatterns <- (!is.null(names(colClasses)));
  if (!is.null(columns) && hasColClassPatterns) {
    nbrOfColumns <- length(columns);
    defColClasses <- rep(defColClass, times=nbrOfColumns);
    defColClassPatterns <- defColClasses;

    colClassPatterns <- colClasses;
    names <- names(colClassPatterns);
    colClasses <- NULL ## Not needed anymore

    # If the same column class name is specified more than ones, then
    # let the latter override the former
    dups <- rev(duplicated(rev(names)))
    if (any(dups)) {
      colClassPatterns <- colClassPatterns[!dups]
      names <- names(colClassPatterns)
    }

    # Default columns?
    pos <- which(names == "*");
    if (length(pos) > 0) {
      # Exclude extra '*':s
      if (length(pos) > 1) {
        colClassPatterns <- colClassPatterns[-(pos[-1])];
        pos <- pos[1];
      }

      # Insert defaults
      colClass <- colClassPatterns[pos];
      names <- names(colClassPatterns);
      if (length(colClassPatterns) > 1) {
        names <- insert(names[-pos], at=pos, values=rep("*", times=nbrOfColumns));
        idxs <- which(names == "*");
        names[idxs] <- sprintf("^%s$", columns);

        colClassPatterns <- insert(colClassPatterns[-pos], at=pos,
                                   values=rep("*", times=nbrOfColumns));
        names(colClassPatterns) <- names;
        colClassPatterns[idxs] <- colClass;
      } else {
        colClassPatterns <- rep(colClass, times=nbrOfColumns);
        names(colClassPatterns) <- sprintf("^%s$", columns);
      }
    } # if (length(pos) > 0)

    verbose && cat(verbose, "Pattern used to identify column classes:", level=-20);
    verbose && print(verbose, colClassPatterns, level=-20);

    verbose && cat(verbose, "Generate column classes:");
    # Read everything by default
    colClasses <- defColClasses;
    names(colClasses) <- columns;

    # Update column classes according to patterns
    for (kk in seq_along(colClassPatterns)) {
      pattern <- names(colClassPatterns)[kk];
      idxs <- which(regexpr(pattern, columns) != -1);
      if (length(idxs) > 0) {
        colClass <- colClassPatterns[kk];
        colClasses[idxs] <- colClass;
      }
    } # for (kk ...)
    colClassPatterns <- NULL ## Not needed anymore
  } else {
    # If column names cannot be inferred, use the default
    nbrOfColumns <- length(fileHeader$topRows[[1]]);
    # Where argument 'colClasses' explicitly specified?
    if (!missing(colClasses) && !is.null(colClasses)) {
      # Sanity check
      if (length(colClasses) != nbrOfColumns) {
        throw(sprintf("The number of elements in argument 'colClasses' does not match the number of column in the file: %d != %d", length(colClasses),  nbrOfColumns));
      }
    } else {
      colClasses <- rep(colClasses["*"], times=nbrOfColumns);
    }
  } # if (!is.null(columns) && hasColClassPatterns)

  verbose && cat(verbose, "Column classes:", level=-20);
  verbose && print(verbose, colClasses, level=-20);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup read.table() arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Inferred arguments
  ch <- getCommentChar(this)
  args <- list(
    header       = hasColumnHeader(this),
    colClasses   = colClasses,
    skip         = fileHeader$skipMax,
    sep          = fileHeader$sep,
    quote        = fileHeader$quote,
    fill         = this$fill,
    comment.char = ifelse(is.null(ch), "", ch),
    check.names  = FALSE,
    na.strings   = c("---", "NA")
  );

  # Overwrite with user specified arguments, if any
  if (length(userArgs) > 0) {
    verbose && enter(verbose, "Overwriting inferred arguments with user arguments");
    for (key in names(userArgs)) {
      args[[key]] <- userArgs[[key]];
    }
    verbose && exit(verbose);
  }


  # Drop NULL arguments
  keep <- !sapply(args, FUN=is.null);
  args <- args[keep];

  verbose && exit(verbose);


  args;
}, protected=TRUE) # getReadArguments()




###########################################################################/**
# @RdocMethod readDataFrame
#
# @title "Reads the tabular data as a data frame"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{con}{(Internal) If a @connection, then it is used, otherwise
#   a new file connection is temporarly opened and used.}
#   \item{rows}{(Optional) An @integer @vector specifying which rows to
#    be read.}
#   \item{nrow}{(Optional) An @integer specifying how many rows to read.
#    If specified, it corresponds to specifying \code{rows=seq_len(nrow)}.}
#   \item{trimQuotes}{(Optional) If @TRUE, quotes are trimmed from numeric
#    columns before parsing them as numerics.  This makes it possible to
#    read quoted numeric values.}
#   \item{...}{Passed to internal @seemethod "getReadArguments".}
#   \item{debug}{If @TRUE, additional details on the file and how it was
#    read is returned as part of the attributes.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a @data.frame.
# }
#
# \section{Reading quoted numerics}{
#   If a specific data column is specified as being numeric in
#   argument \code{colClasses} and that column contains quoted values
#   it is necessary to use argument \code{trimQuotes=TRUE}, otherwise
#   @see "base::scan" throws an exception similar to:
#   \code{scan() expected 'a real', got '"1.0"'}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("readDataFrame", "TabularTextFile", function(this, con=NULL, rows=NULL, nrow=NULL, trimQuotes=FALSE, ..., debug=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'rows' and 'nrow':
  if (!is.null(rows) && !is.null(nrow)) {
    throw("Only one of arguments 'rows' and 'nrow' can be specified.");
  }

  if (!is.null(rows)) {
    rows <- Arguments$getIndices(rows);
    nrow <- max(rows);
  }

  if (!is.null(nrow)) {
    nrow <- Arguments$getInteger(nrow, range=c(1,Inf));
  }

  # Argument 'trimQuotes':
  trimQuotes <- Arguments$getLogical(trimQuotes);

  # Argument 'debug':
  debug <- Arguments$getLogical(debug);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading ", class(this)[1L]);


  # Attributes to be added
  attributes <- list();


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading header to infer read.table() arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hdr <- getHeader(this, ..., verbose=less(verbose, 5));
  attributes$header <- hdr;

  # Get read arguments
  args <- getReadArguments(this, fileHeader=hdr, nrow=nrow, ...,
                                               verbose=less(verbose, 5));

  verbose && cat(verbose, "Arguments inferred from file header:");
  verbose && print(verbose, args);
  attributes$readArguments <- args;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify names of columns read
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  columns <- getColumnNames(this);
  verbose && printf(verbose, "Column names (%d):\n", length(columns));
  verbose && cat(verbose, paste(columns, collapse=", "));

  if (!is.null(columns)) {
    verbose && enter(verbose, "Matching column names:");
    verbose && printf(verbose, "Column classes (%d):\n",
                                                length(args$colClasses));
    verbose && cat(verbose, paste(args$colClasses, collapse=", "));
    columns[args$colClasses == "NULL"] <- NA;
    columns <- columns[!is.na(columns)];
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Open a file connection?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(con)) {
    pathname <- getPathname(this);
    verbose && cat(verbose, "Pathname: ", pathname);

    # Open file connection
    con <- file(pathname, open="r");
    on.exit({
      if (!is.null(con)) {
        close(con);
        con <- NULL;
      }
    })
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calling read.table()");


  # SPECIAL CASE/WORKAROUND: read.table()/scan() will give an error
  # if a numeric value is quoted and 'colClasses' specifies it as
  # a numeric value.  In order to read such values, we need to remove
  # the quotes first. /HB 2011-07-13

  # Check if we need to trim quotes
  trimQuotes <- (trimQuotes && nchar(args$quote) > 0L);
  trimQuotes <- (trimQuotes && length(args$colClasses) > 0L);
  if (trimQuotes) {
    classesToPatch <- c("integer", "numeric", "double", "complex");
    toPatch <- is.element(args$colClasses, classesToPatch);
    trimQuotes <- c(trimQuotes && any(toPatch));
  }

  # Read by first trimming quotes?
  if (trimQuotes) {
    verbose && enter(verbose, "Trimming quotes from numeric columns before parsing them as numbers");
    # This is a workaround for reading quoted numerics. /HB 2011-07-13
    verbose && cat(verbose, "Columns that need to have quotes trimmed:");
    verbose && str(verbose, which(toPatch));
    args0 <- args;

    # (1) Read to-be-patched columns as characters.
    colClasses <- args$colClasses;
    colClasses[toPatch] <- "character";
    args$colClasses <- colClasses;
  }

  verbose && cat(verbose, "Arguments used to read tabular file:");
  ## Since R v3.2.1 rev 68837, read.table() now utilized (fixed)
  ## named 'colClasses' similar to readDataFrame().  For now, make
  ## sure to pass non-named 'colClasses' to read.table().
  ## FIX ME: Should we rewrite readDataFrame() to better utilize
  ## that read.table() now also looks at fixed names? /HB 2015-08-06
  args$colClasses <- unname(args$colClasses);
  args <- c(list(con), args);
  verbose && print(verbose, args);
  data <- do.call(read.table, args=args);
  nbrOfRowsRead <- nrow(data);
  verbose && cat(verbose, "Raw data read by read.table():");
  verbose && str(verbose, data);
  verbose && cat(verbose, "Number of rows read: ", nbrOfRowsRead);
  verbose && cat(verbose, "Number of columns read: ", ncol(data));

  # Extract subset of rows?
  if (!is.null(rows)) {
    if (max(rows) > nbrOfRowsRead) {
      rows <- intersect(rows, 1:nbrOfRowsRead);
      warning("Argument 'rows' was out of range [1,", nbrOfRowsRead,
                                "]. Ignored rows beyond this range.");
    }
    data <- data[rows,,drop=FALSE];
  } else {
    rownames(data) <- NULL;
  }

  # Was data read by first trimming quotes?
  if (trimQuotes) {
    # (2) Re-read numeric columns one by one
    colClasses <- args0$colClasses;

    # Note that some columns may have been ignored
    keep <- (colClasses != "NULL");
    colClasses <- colClasses[keep];
    toPatch <- toPatch[keep];

    # Sanity check
    stopifnot(length(colClasses) == ncol(data));

    na.strings <- args0$na.strings;

    verbose && enter(verbose, "Parsing numeric columns");
    toPatch <- which(toPatch);
    for (kk in seq_along(toPatch)) {
      col <- toPatch[kk];
      colClass <- colClasses[col];
      verbose && enter(verbose, sprintf("Parsing #%d (column #%d as '%s') of %d", kk, col, colClass, length(toPatch)));
      values <- data[[col]];

      # Quick/bad way, but will turn non-valid values into NA
      ##  storage.mode(values) <- colClass;

      verbose && cat(verbose, "Non-parsed values:");
      verbose && str(verbose, values);

      bfr <- paste(values, collapse="\n");
      conT <- textConnection(bfr, open="r");
      on.exit({
        if (!is.null(conT)) close(conT);
      }, add=TRUE);

      # Try to read the values as the correct type.
      # (Could scan(), which has less overhead, be used here? /HB 2013-09-30)
      valuesT <- read.table(file=conT, quote="", colClasses=colClass, na.strings=na.strings, blank.lines.skip=FALSE)[[1L]];

      verbose && cat(verbose, "Parsed values:");
      verbose && str(verbose, valuesT);

      # Sanity check
      stopifnot(length(valuesT) == length(values));

      close(conT); conT <- NULL;
      data[[col]] <- valuesT;

      # Not needed anymore
      bfr <- values <- valuesT <- NULL;
      verbose && exit(verbose);
    } # for (kk ...)
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # if (trimQuotes)


  # Sanity check
  if (length(columns) > 0L) {
    if (ncol(data) != length(columns)) {
      throw("Number of read data columns does not match the number of column headers: ", ncol(data), " != ", length(columns));
    }
    colnames(data) <- columns;
  }

  verbose && str(verbose, data);
  verbose && exit(verbose);

  if (debug) {
    attr(data, "fileHeader") <- hdr;
    for (key in names(attributes)) {
      attr(data, key) <- attributes[[key]];
    }
  }

  verbose && exit(verbose);

  data;
}) # readDataFrame()



setMethodS3("readColumns", "TabularTextFile", function(this, columns=seq_len(ncol(this)), colClasses=NULL, ..., check.names=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'columns':
  maxNbrOfColumns <- ncol(this);
  if (is.null(columns)) {
    columnNames <- getColumnNames(this);
  } else if (is.numeric(columns)) {
    columns <- Arguments$getIndices(columns, max=maxNbrOfColumns);
    columnNames <- getColumnNames(this);
    if (!is.null(columnNames)) {
      columnNames <- columnNames[columns];
    }
  } else {
    columnNames <- Arguments$getCharacters(columns);
  }

  # Argument 'colClasses':
  if (!is.null(colClasses)) {
    names <- names(colClasses);
    colClasses <- Arguments$getCharacters(colClasses);
    names(colClasses) <- names;
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading columns");
  verbose && cat(verbose, "Argument 'columns': ", hpaste(columns));
  verbose && cat(verbose, "Argument 'colClasses':");
  verbose && print(verbose, colClasses);

  verbose && cat(verbose, "Column names': ", hpaste(columnNames));
  # Setup column classes, iff missing
  if (is.null(colClasses)) {
    if (is.null(columnNames)) {
      colClasses[-columns] <- "NULL";
    } else {
      colClasses <- rep("character", times=length(columnNames));
      names(colClasses) <- sprintf("^%s$", columnNames);
    }
  } else {
    if (is.null(names(colClasses))) {
      names(colClasses) <- sprintf("^%s$", columnNames);
    }
  }

  verbose && cat(verbose, "Column classes:");
  verbose && print(verbose, colClasses);

  data <- readDataFrame(this, colClasses=colClasses, ..., verbose=less(verbose, 50));

  # Subset
  if (ncol(data) > 0L) {
    verbose && enter(verbose, "Subsetting columns");
    verbose && str(verbose, data);

    if (!is.null(columnNames)) {
      columns <- match(columnNames, names(data));
    }
    verbose && cat(verbose, "Columns to keep: ", hpaste(columns));

    # Sanity check
    columns <- Arguments$getIndices(columns, max=ncol(data));

    # Need to rearrange?
    if (any(diff(columns) != 1L) || ncol(data) > length(columns)) {
      data <- data[,columns,drop=FALSE];
      if (!check.names) {
        colnames(data) <- columnNames;
      }
    }

    verbose && exit(verbose);
  }

  # Sanity check
  stopifnot(ncol(data) == length(columns));

  verbose && exit(verbose);

  data;
}, protected=TRUE)



# AD HOC fix to speed up ll(), which calls dimension() on each object,
# which in turn calls dim() and dim() is really slow for this class,
# because it has to infer the number of rows by reading the complete
# file. The fix is to return NA for the number of rows if the file size
# is larger than 10MB unless nbrOfLines() has already been called.
# /HB 2008-07-22
setMethodS3("dimension", "TabularTextFile", function(this, ...) {
  c(nbrOfRows(this, fast=TRUE), nbrOfColumns(this, fast=TRUE));
}, private=TRUE);



###########################################################################/**
# @RdocMethod nbrOfRows
#
# @title "Counts the number of data rows"
#
# \description{
#  @get "title".  The count does not include the header rows or comments.
# }
#
# @synopsis
#
# \arguments{
#   \item{fast}{Argument passed to @seemethod "nbrOfLines".}
#   \item{...}{Optional arguments passed to @seemethod "getHeader".}
# }
#
# \value{
#   Returns a @character @vector.
# }
# @author
#
# \seealso{
#    The number of data rows is always less or equal to the number of lines
#    as returned by @seemethod "nbrOfLines".
#   Internally, @see "R.utils::countLines" is used.
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("nbrOfRows", "TabularTextFile", function(this, fast=FALSE, ...) {
  hdr <- getHeader(this, ...);

  nbrOfCommentRows <- length(hdr$comments);
  nbrOfRowsToSkip <- hdr$skip;
  hasColumnNames <- (length(hdr$columns) > 0L);
  nbrOfNonDataRows <- nbrOfCommentRows + nbrOfRowsToSkip +
                      as.integer(hasColumnNames);

  n <- nbrOfLines(this, fast=fast);
  n <- n - nbrOfNonDataRows;
  n <- as.integer(n);

  n;
})



###########################################################################/**
# @RdocMethod nbrOfLines
#
# @title "Counts the number of lines in the data file"
#
# \description{
#  @get "title".  The count include header rows, comments and more.
# }
#
# @synopsis
#
# \arguments{
#   \item{fast}{If @TRUE, @NA is returned for large data files (>1Mb),
#     unless the number of lines has already been counted.}
#   \item{...}{Optional arguments passed to @see "R.utils::countLines".}
# }
#
# \value{
#   Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#    To count the number of data rows is the data table,
#    use @seemethod "nbrOfRows".
#   Internally, @see "R.utils::countLines" is used.
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("nbrOfLines", "TabularTextFile", function(this, fast=FALSE, ...) {
  pathname <- getPathname(this);

  n <- this$.nbrOfLines;
  mtime <- getLastModifiedOn(this);
  if (is.null(n) || is.na(mtime) || !identical(mtime, attr(n, "mtime"))) {
    if (fast) {
      if (getFileSize(this) < 10e6)
        fast <- FALSE;
    }

    if (fast) {
      n <- as.integer(NA);
    } else {
      n <- countLines(pathname, ...);
      attr(n, "mtime") <- mtime;
      this$.nbrOfLines <- n;
    }
  }

  n;
})




###########################################################################/**
# @RdocMethod readLines
#
# @title "Reads the lines of the data file as strings"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Optional arguments passed to @see "base::readLines".}
# }
#
# \value{
#   Returns a @character @vector.
# }
# @author
#
# \seealso{
#   @seemethod "readDataFrame".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("readLines", "TabularTextFile", function(con, ...) {
  # To please R CMD check
  this <- con;
  pathname <- getPathname(this);
  readLines(pathname, ...);
})




############################################################################
# HISTORY:
# 2015-05-04
# o BUG FIX: Now getReadArguments() for TabularTextFile let duplicated
#   named 'colClasses' entries override earlier ones, e.g.
#   colClasses=c("*"=NA, "*"="NULL", a="integer") is effectively the
#   same as colClasses=c("*"="NULL", a="integer").
# 2014-08-25
# o Now readColumns(..., colClasses) adds regexpr names to 'colClasses',
#   iff missing.
# o BUG FIX: readDataFrame() would ignore argument 'colClasses'
#   iff it had no names.
# 2014-08-23
# o BUG FIX: Using commentChar=NULL for TabularTextFile:s failed.
# 2014-08-02
# o Added argument 'sep' to readRawHeader().
# 2014-01-24
# o Now readColumns() for TabularTextFile handles also header-less files.
# 2013-12-18
# o BUG FIX: Now getReadArguments() for TabularTextFile returns a
#   'colClasses' vector of the correct length also in the case when
#   there are no column names.
# 2013-09-30
# o Now readDataFrame() for TabularTextFile subsets by row, if requested,
#   before reparsing numerical columns that were quoted.
# 2013-09-23
# o SPEEDUP/CLEANUP: Package no longer uses R.utils::whichVector(), which
#   use to be 10x faster, but since R 2.11.0 which() is 3x times faster.
# 2013-01-17
# o In addition to a fixed integer, argument 'skip' for readDataFrame()
#   (default and for TabularTextFile) may also specify a regular
#   expression matching the first row of the data section.
# o SPEEDUP: Now getReadArguments() for TabularTextFile returns 'skip'
#   as the maximum number of rows possible for read.table() to skip,
#   which includes the original skip plus all skipped comment rows.
# 2013-01-16
# o Now 'skip' simply skips the first 'skip' lines before parsing
#   the header or the data.  This is how read.table() works.
# o Added argument 'skip' to readRawHeader(..., skip=this$skip).
# o Now arguments '...' to readDataFrame() for TabularTextFile are passed
#   to getHeader().
# o Now readDataFrame(..., debug=TRUE) returns also the read arguments.
# 2012-12-20
# o Renamed argument 'colClassPatterns' of getReadArguments() for
#   TabularTextFile to 'colClasses'.  However, if the old name is
#   still supported.
# 2012-12-08
# o BUG FIX: nbrOfRows() for TabularTextFile forgot to exclude comment
#   rows in the file header.
# o BUG FIX: readColumns() for TabularTextFile would not preserve the
#   order of the requested 'columns'.
# 2012-12-02
# o BUG FIX: getDefaultColumnNames() for TabularTextFile did not use
#   'columnNames' if it was set when creating the TabularTextFile object.
# o BUG FIX: Now getReadArguments() for TabularTextFile drops arguments
#   that are NULL, because they could cause errors downstreams, e.g.
#   readDataFrame() calling read.table(..., colClasses=NULL) =>
#   rep_len(NULL, x) => "Error in rep_len(colClasses, cols) :
#   cannot replicate NULL to a non-zero length".
# 2012-11-28
# o Declaring '.fileHeader' and '.nbrOfLines' as 'cached' fields.
# 2012-11-15
# o Made it possible for TabularTextFile to ignore header comment
#   arguments when inferring column names and classes.
# 2012-11-08
# o Now getReadArguments() for TabularTextFile also includes column
#   class patterns from getDefaultColumnClassPatterns().
# o Added getDefaultColumnClass() and getDefaultColumnClassPatterns()
#   to TabularTextFile.
# o Now getDefaultColumnNames() for TabularTextFile falls back to
#   header comment argument 'columnNames', if there are no column
#   names in the actual data table.
# o Now readRawHeader() for TabularTextFile also parses and returns
#   header comment arguments.
# 2012-11-02
# o Added getDefaultColumnNames() for TabularTextFile and dropped
#   getColumnNames(), which is implemented by ColumnNamesInterface.
# 2012-10-31
# o CLEANUP: Now readDataFrame() for TabularTextFile no longer returns
#   attribute 'fileHeader', unless argument 'debug' is TRUE.
# 2012-09-27
# o ROBUSTNESS: Now getHeader() for TabularTextFile checks if the file
#   has been modified before returned cached results.
# o Added argument 'stringsAsFactors=FALSE' to getReadArguments() for
#   TabularTextFile such that the default is to read strings as character
#   rather than as factors.
# 2011-09-26
# o Added methods set- and getCommentChar() to TabularTextFile and
#   argument 'commentChar' to its constructor.  This allows to use
#   custom comment characters other than just "#".
# 2011-07-13
# o GENERALIZATION: Now readDataFrame() of TabularTextFile can read
#   numeric columns that are quoted and for which 'colClasses' in non-NA.
#   This is done by first reading the as quoted character strings, and
#   dropping the quotes, and then rereading them as numeric values.
# 2011-05-12
# o Added more verbose output to readDataFrame().
# 2011-03-14
# o Improved the verbose output of getReadArguments().
# 2010-08-14
# o Added Rdoc comments.
# 2010-04-22
# o Added "NA" to the default 'na.strings' returned by getReadArguments()
#   for TabularTextFile.
# 2009-10-06
# o Added subsetting via [() to TabularTextFile.
# 2009-05-09
# o Added argument 'translate' to getColumnNames() of TabularTextFile.
# 2008-07-23
# o Now nbrOfLines() cache the results and only recount if the file has
#   been modified since last time (or the file system does not provide
#   information on last modification time).  It also uses the new
#   countLines() in R.utils.
# 2008-07-22
# o Added ad hoc dimension() to speed up ll(). It may return NA for the
#   number of rows.
# 2008-06-12
# o Added readRawHeader().  Removed readHeader() [now in getHeader()].
# o Added hasColumnHeader() to TabularTextFile.
# 2008-05-16
# o Took the text-file based parts of GenericTabularFile and placed them
#   in new subclass TabularTextFile.
# 2008-05-12
# o Added extractMatrix().
# o BUG FIX: getReadArguments() did not infer column classes if there was
#   no header to read but the column names was manually set.
# o BUG FIX: readDataFrame() did not read the first data row if there was
#   no column header; it was eaten up by a preceeding readHeader().
# 2008-04-29
# o Added readLines(), nbrOfLines(), nbrOfRows() and dim().
# o Now readDataFrame() keeps the row names if arguments rows != NULL.
# 2008-04-25
# o Now argument 'verbose' of the constructor is passed to verfity().
# 2008-04-24
# o Added argument 'rows' to readDataFrame() for TabularTextFile.
# 2008-04-14
# o Renamed readData() to readDataFrame() for TabularTextFile.
# 2008-03-22
# o Added {get|set}ColumnNameTranslator().
# 2008-03-18
# o Now any '...' arguments to getReadArguments() override the inferred
#   read arguments, e.g. na.strings="NA".
# 2008-02-27
# o Since 'affy' defines standardGeneric("colnames") and because S3 methods
#   are not found by such S4 generic functions, we avoid that method name,
#   and instead use getColumnNames().
# 2007-09-16
# o Removed all 'camelCaseNames' arguments.  Now column names are decided
#   by getColumnNames() and translateColumnNames(), which can be overridden.
# 2007-09-14
# o Extracted from AffymetrixTabularFile.
# 2007-09-10
# o Created from AffymetrixCsvGenomeInformation.R.
############################################################################
