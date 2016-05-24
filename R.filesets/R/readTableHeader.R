setMethodS3("readTableHeader", "default", function(con, sep=c("\t", ","), quote=FALSE, comment.char="#", skip=0, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'con':
  if (is.character(con)) {
  } else if (inherits(con, "connection")) {
  } else {
    throw("Argument 'con' must be a connection or a pathname: ", 
                                                          class(con)[1]);
  }

  # Argument 'sep':
  sep <- Arguments$getCharacters(sep);

  # Argument 'quote':
  if (is.null(quote)) {
  } else {
    quote <- Arguments$getCharacter(quote);
  }

  # Argument 'skip':
  skip <- Arguments$getInteger(skip, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading header of tabular file");

  # Open a file connection?
  if (is.character(con)) {
    pathname <- con;
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

  ready <- FALSE;
  comments <- c();
  pattern <- sprintf("^%s", comment.char);
  while (!ready) {
    line <- readLines(con, n=1);
    isComments <- (regexpr(pattern, line) != -1);
    if (!isComments) {
      if (skip == 0)
        break;
      skip <- skip - 1;
    }
    comments <- c(comments, line);
  }

  verbose && cat(verbose, "Header comments:", level=-20);
  verbose && str(verbose, comments, level=-20);

  # Infer column separator?
  if (length(sep) > 1) {
    verbose && enter(verbose, "Identifying the separator that returns most columns");
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
#  topRows <- lapply(topRows, trim);
  verbose && print(verbose, topRows);

  # Remove quotes?
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
    sep=sep,
    quote=quote,
    skip=skip,
    topRows=topRows
  );

  verbose && str(verbose, hdr);

  verbose && exit(verbose);

  hdr;
}) # readTableHeader()



############################################################################
# HISTORY:
# 2008-06-14
# o Created from readRawHeader() of TabularTextFile in aroma.core.
############################################################################ 
