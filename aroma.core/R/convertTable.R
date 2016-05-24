setMethodS3("convertTable", "default", function(src, dest=NULL, sep="\t", ..., verbose=FALSE) {
  # Argument 'src':
  src <- Arguments$getReadablePathname(src, mustExist=TRUE);

  # Argument 'dest':
  dest <- Arguments$getWritablePathname(dest);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Convert tabular file");
  verbose && cat(verbose, "Source: ", src);
  verbose && cat(verbose, "Destination: ", dest);

  inCon <- file(src, open="r");
  on.exit(close(inCon));
  outCon <- file(dest, open="w");
  on.exit(close(outCon), add=TRUE);

  # Header
  header <- readLines(con=inCon, n=1);
  writeLines(con=outCon, header);

  # Number of columns
  verbose && cat(verbose, "Header: ", header);
  nbrOfColumns <- length(strsplit(header, split=sep)[[1]]);
  verbose && cat(verbose, "Number of columns: ", nbrOfColumns);

  verbose && enter(verbose, "Padding the rows for equal number of columns");
  # Predefine possible paddings
  upad <- sapply(1:nbrOfColumns, FUN=function(n) {
    paste(rep(sep, n), collapse="");
  })
  names(upad) <- 1:nbrOfColumns;

  # Convert file in blocks of 100,000 rows.
  blockSize <- 1e5;
  count <- 1;
  while (TRUE) {
    verbose && printf(verbose, "Block #%d\n", as.integer(count));
    lines <- readLines(con=inCon, n=blockSize);
    if (length(lines) == 0)
      break;

    # Count the number of column in each row
    widths <- gregexpr(sep, lines, fixed=TRUE);
    widths <- unlist(lapply(widths, FUN=length), use.names=FALSE);
    widths <- widths + 1;

    # Fill with TABs if too few
    missing <- nbrOfColumns - widths;
    nok <- which(missing > 0);
    if (length(nok) > 0) {
      lines[nok] <- paste(lines[nok], upad[missing[nok]], sep="");
    }
    # Not needed anymore
    missing <- nok <- NULL;

    writeLines(con=outCon, lines);
    # Not needed anymore
    lines <- NULL;

    count <- count + 1;
  }
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(dest);
})



############################################################################
# HISTORY:
# 2006-09-15
# o This methods will probably end up in R.utils at some stage.
# o Added convertTable().  The different dChip files are ragged rows! :(
############################################################################
