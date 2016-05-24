###########################################################################/**
# @RdocDefault compileAsciiDoc
#
# @title "Compiles an AsciiDoc file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename, path}{The filename and (optional) path of the
#      document to be compiled.}
#   \item{...}{Additional arguments passed to executable \code{ascii}
#     (which must be on the system search path)
#     called via @see "base::system2".}
#   \item{outPath}{The output and working directory.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the pathname of the generated document.
# }
#
# @author
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("compileAsciiDoc", "default", function(filename, path=NULL, ..., outPath=".", postprocess=TRUE, verbose=FALSE) {
  use("ascii", quietly=TRUE);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' & 'path':
  pathname <- if (is.null(path)) filename else file.path(path, filename);
  if (!isUrl(pathname)) {
    pathname <- Arguments$getReadablePathname(pathname);
  }

  # Arguments 'outPath':
  outPath <- Arguments$getWritablePath(outPath);
  if (is.null(outPath)) outPath <- ".";

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Compiling AsciiDoc noweb document");
  # Download URL?
  if (isUrl(pathname)) {
    verbose && enter(verbose, "Downloading URL");
    url <- pathname;
    verbose && cat(verbose, "URL: ", url);
    pathname <- downloadFile(url, verbose=less(verbose,50));
    verbose && cat(verbose, "Local file: ", pathname);
    verbose && exit(verbose);
  }

  pathname <- getAbsolutePath(pathname);
  verbose && cat(verbose, "Pathname (absolute): ", pathname);
  verbose && printf(verbose, "Input file size: %g bytes\n", file.info(pathname)$size);
  verbose && cat(verbose, "Output and working directory: ", outPath);

  bin <- findAsciiDoc(mustExist=TRUE, verbose=less(verbose, 10));

  opwd <- ".";
  on.exit(setwd(opwd), add=TRUE);
  if (!is.null(outPath)) {
    opwd <- setwd(outPath);
  }

  args <- c("-v");
##  args <- c(args, "-a data-uri");

  # Output file
  fullnameR <- gsub("[.][^.]*$", "", basename(pathname));
  extR <- "html";
  filenameR <- sprintf("%s.%s", fullnameR, extR);
  pathnameR <- filenameR;
  arg <- sprintf("-o %s", pathnameR);
  args <- c(args, arg);

  # Source file
  pathname <- normalizePath(pathname);
  args <- c(args, pathname);

  if (verbose && isVisible(less(verbose, 5))) {
     verbose && cat(verbose, "AsciiDoc executable:");
     verbose && print(verbose, bin);
     verbose && cat(verbose, "Command-line arguments:");
     verbose && print(verbose, args);
  }

  res <- system2(bin, args=args, stderr=TRUE);
  res <- trim(res);
  if (length(res) == 0L) {
    throw("Failed to run external 'asciidoc'. No output is available.");
  }

  if (verbose && isVisible(less(verbose, 10))) {
     verbose && cat(verbose, res, collapse="\n");
  }

  # Locate asciidoc warnings
  pattern <- "^asciidoc: WARNING: ";
  warns <- grep(pattern, res, value=TRUE);
  if (length(warns) > 0L) {
    warns <- gsub(pattern, "", warns);
    warns <- trim(warns);
    verbose && enter(verbose, sprintf("Detected %d AsciiDoc warnings", length(warns)));
    verbose && cat(verbose, warns, collapse="\n");
    for (warn in warns) warning("AsciiDoc WARNING: ", warn);
    verbose && exit(verbose);
  }

  # Locate output filename
  pattern <- "^asciidoc: writing: ";
  pathname2 <- grep(pattern, res, value=TRUE);
  pathname2 <- gsub(pattern, "", pathname2);
  pathname2 <- trim(pathname2);
  pathname2 <- getAbsolutePath(pathname2);
  setwd(opwd); opwd <- ".";

  res <- RspFileProduct(pathname2);
  verbose && print(verbose, res);

  # Postprocess?
  if (postprocess) {
    res <- process(res, outPath=outPath, recursive=TRUE, verbose=verbose);
  }

  verbose && exit(verbose);

  res;
}) # compileAsciiDoc()


############################################################################
# HISTORY:
# 2013-03-29
# o Created (from compileAsciiDocNoweb.R).
############################################################################
