###########################################################################/**
# @RdocDefault translateRsp
#
# @title "Translates an RSP file to an R RSP source file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{A filename to be read.}
#   \item{path}{An optional path to the file.}
#   \item{...}{Not used.}
#   \item{force}{A @logical.}
#   \item{verbose}{@see "R.utils::Verbose".}
# }
#
# \value{
#   Returns (invisibly) the pathname to the R RSP source code.
# }
#
# @author
#
# \seealso{
#   Internally @see "parseRsp" parses the RSP file into an R code string.
#   @see "sourceRsp".
# }
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("translateRsp", "default", function(filename, path=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ruler <- function(char="#", width=70, ...) {
    ruler <- rep(char, length=width);
    ruler <- paste(ruler, collapse="");
    ruler <- substring(ruler, 1, width);
    ruler;
  } # ruler()

  comment <- function(..., prefix="#", suffix="\n", collapse="") {
    text <- c(...);
    comment <- paste(prefix, text, suffix, sep="");
    comment <- paste(comment, collapse=collapse);
    comment;
  } # comment();

  banner <- function(..., char="#") {
    comment(ruler(char), ..., ruler(char));
  } # banner()

  packageRCode <- function(rCode, pathname=NULL, collapse="", ...) {
    rCodeOrg <- rCode;

    rCode <- banner(
      " DO NOT EDIT!  DO NOT EDIT!  DO NOT EDIT!  DO NOT EDIT!  DO NOT EDIT!",
      "",
      " This R code was parsed from RSP by the R.rsp package."
    );

    code <- "# Assert that write() of R.rsp is used below\n";
    rCode <- c(rCode, code);
    code <- "write <- R.rsp:::write;\n";
    rCode <- c(rCode, code);

    code <- "# Sets the public RspPage 'page' object\n";
    rCode <- c(rCode, code);
    if (is.null(pathname)) {
      code <- "page <- RspPage(pathname=NULL);\n";
      code <- sprintf(fmtstr, "NULL");
    } else {
      fmtstr <- "page <- RspPage(pathname=\"%s\");\n";
      code <- sprintf(fmtstr, pathname);
    }
    rCode <- c(rCode, code);

    code <- "# Gets the output connection (or filename) for the response [OBSOLETE]\n";
    rCode <- c(rCode, code);
    code <- "out <- getOutput(response);\n";
    rCode <- c(rCode, code);

    code <- banner(char="=",
      " BEGIN RSP CODE"
    );
    rCode <- c(rCode, code);
    rCodeHeader <- rCode;

    code <- banner(char="=",
      " END RSP CODE"
    );
    rCodeFooter <- code;

    rCode <- c(rCodeHeader, rCodeOrg, rCodeFooter);

    rCode <- paste(rCode, collapse=collapse);

    attr(rCode, "pathname") <- pathname;

    rCode;
  } # packageRCode()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' & 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path,
                                                         mustExist=TRUE);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Translates an RSP document in to an R RSP source file");

  # Assert correct filename extension
  pattern <- "[.](rsp|RSP)$";
  if (regexpr(pattern, pathname) == -1) {
    throw("Cannot compile file. File does not have extension rsp: ", pathname);
  }

  path <- dirname(pathname);
  filename <- basename(pathname);

  # Setup output pathname name
  overwrite <- TRUE;
  outPath <- path;
  outFilename <- sprintf("%s.R", filename);
  outPathname <- Arguments$getWritablePathname(outFilename, path=outPath, overwrite=overwrite);

  # Check if an up-to-date output file already exists
  isUpToDate <- FALSE;
  if (!force && isFile(outPathname)) {
    date <- file.info(pathname)$mtime;
    verbose && cat(verbose, "Source file modified on: ", date);
    outDate <- file.info(outPathname)$mtime;
    verbose && cat(verbose, "Output file modified on: ", outDate);
    if (is.finite(date) && is.finite(outDate)) {
      isUpToDate <- (outDate >= date);
    }
    verbose && printf(verbose, "Output file is %sup to date.\n", if(!isUpToDate) "not " else "");
  }

  errPathname <- sprintf("%s.ERROR", outPathname);

  if (force || !isUpToDate) {
    verbose && enter(verbose, "Parses the RSP file");
    verbose && cat(verbose, "Pathname: ", pathname);

    # Read RSP code
    rspCode <- readLines(pathname, warn=FALSE);
    verbose && str(verbose, rspCode);

    # Compile RSP to output file
    tryCatch({
      rCode <- parseRsp(rspCode, verbose=less(verbose, 10));
    }, error = function(ex) {
      rCode <- parseRsp(rspCode, validate=FALSE, verbose=less(verbose, 10));
      cat(file=errPathname, rCode);
      ex$message <- paste(ex$message, "\nThe translated RSP code is not valid: ", errPathname, sep="");
      verbose && cat(verbose, "The translated RSP code is not valid: ", errPathname);
      stop(ex);
    })

    rCode <- packageRCode(rCode, pathname=pathname);
    cat(file=outPathname, rCode);
    rCode <- NULL; # Not needed anymore

    verbose && cat(verbose, "R RSP file generated: ", outPathname);
    verbose && exit(verbose);
  }

  if (isFile(errPathname)) {
    file.remove(errPathname);
  }

  verbose && exit(verbose);

  invisible(outPathname);
}) # translateRsp()


###########################################################################
# HISTORY:
# 2014-10-18
# o CLEANUP/ROBUSTNESS: translateRsp() and translateRspV1(), which are
#   both deprecated, no longer assume that write() is exported from R.rsp.
# 2011-11-17
# o Now the generated R script adds 'write <- R.rsp::write' at the
#   beginning, to assure that it is used instead of base::write().
# 2009-02-23
# o Renamed from compileRsp() to translateRsp().
# o Updated to use parseRsp().
# 2009-02-22
# o Created.
###########################################################################
