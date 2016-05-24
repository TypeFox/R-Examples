ps2pdf <- function(filename, path=NULL, opts=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'filename' & 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);

  # Argument 'opts':
  opts <- Arguments$getCharacters(opts);
  opts <- paste(opts, collapse=" ");

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "ps2pdf");
  verbose && cat(verbose, "PS pathname: ", pathname);

  filename <- basename(pathname);
  path <- dirname(pathname);

  filename2 <- gsub("[.]ps", ".pdf", filename);
  pathname2 <- Arguments$getWritablePathname(filename2, path=path);
  verbose && cat(verbose, "PDF pathname: ", pathname2);

  # Is output file up to date?
  isUpToDate <- FALSE;
  if (!force && isFile(pathname2)) {
    date <- file.info(pathname)$mtime;
    verbose && cat(verbose, "Source file modified on: ", date);
    outDate <- file.info(pathname2)$mtime;
    verbose && cat(verbose, "Output file modified on: ", outDate);
    if (is.finite(date) && is.finite(outDate)) {
      isUpToDate <- (outDate >= date);
    }
    verbose && printf(verbose, "Output file is %sup to date.\n", ifelse(isUpToDate, "", "not "));
  }

  if (!isUpToDate || !isFile(pathname2)) {
    verbose && enter(verbose, "Calling ps2pdf");
    # Change working directory
    opwd <- getwd();
    on.exit(setwd(opwd));
    setwd(path);

    verbose && cat(verbose, "Working directory: ", getwd());
    verbose && cat(verbose, "Filename: ", filename);

    cmd <- sprintf("ps2pdf %s %s", opts, filename);
    system(cmd);

    verbose && exit(verbose);
  }

  # Sanity check
  pathname2 <- Arguments$getReadablePathname(pathname2, mustExist=TRUE);

  verbose && exit(verbose);

  invisible(pathname2);
} # ps2pdf()
