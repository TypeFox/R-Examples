setMethodS3("compileRsp0", "default", function(..., envir=parent.frame(), force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'envir':
  if (!is.environment(envir)) {
    throw("Argument 'envir' is not an environment: ", class(envir)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "compileRsp0()");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Translate an RSP file to an R RSP source file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname2 <- translateRsp(..., force=force, verbose=less(verbose,5));
  verbose && cat(verbose, "Translated RSP source code: ", pathname2);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Run R RSP file to generate output document
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path2 <- dirname(pathname2);
  filename2 <- basename(pathname2);

  path3 <- path2;
  filename3 <- gsub("[.](rsp|RSP)[.](r|R)$", "", filename2);
  pathname3 <- Arguments$getWritablePathname(filename3,
                                       path=path3, overwrite=TRUE);

  # Is output file up to date?
  isUpToDate <- FALSE;
  if (!force && isFile(pathname3)) {
    date <- file.info(pathname2)$mtime;
    verbose && cat(verbose, "Source file modified on: ", date);
    outDate <- file.info(pathname3)$mtime;
    verbose && cat(verbose, "Output file modified on: ", outDate);
    if (is.finite(date) && is.finite(outDate)) {
      isUpToDate <- (outDate >= date);
    }
    verbose && printf(verbose, "Output file is %sup to date.\n", ifelse(isUpToDate, "", "not "));
  }

  if (force || !isUpToDate) {
    response <- FileRspResponse(pathname3, overwrite=TRUE);
    envir$response <- response;
    sourceTo(pathname2, envir=envir);
  }

  verbose && cat(verbose, "Output RSP document: ", pathname3);
  verbose && exit(verbose);


  invisible(pathname3);
}) # compileRsp0()


###########################################################################
# HISTORY:
# 2013-03-31
# o One issue with translateRsp() and compileRsp0() is that they create
#   files in the source directory, which may be read-only.
# 2013-03-29
# o Renamed to compileRsp0(). May be dropped rather soon.
# 2009-02-23
# o Updated to use parseRsp().
# 2009-02-22
# o Created.
###########################################################################
