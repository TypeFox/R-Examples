###########################################################################/**
# @RdocDefault compileAsciiDocNoweb
#
# @title "Compiles an AsciiDoc noweb file"
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
#   \item{...}{Additional arguments passed to @see "ascii::Asciidoc".}
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
setMethodS3("compileAsciiDocNoweb", "default", function(filename, path=NULL, ..., outPath=".", postprocess=TRUE, verbose=FALSE) {
  use("ascii", quietly=TRUE);
  # To please R CMD check
  Asciidoc <- NULL; rm(list="Asciidoc");

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

  opwd <- ".";
  on.exit(setwd(opwd), add=TRUE);
  if (!is.null(outPath)) {
    opwd <- setwd(outPath);
  }

  pathname2 <- Asciidoc(pathname);
  pathname2 <- getAbsolutePath(pathname2);
  setwd(opwd); opwd <- ".";

  res <- RspFileProduct(pathname2, type="application/x-asciidoc");
  verbose && print(verbose, res);

  # Postprocess?
  if (postprocess) {
    res <- process(res, outPath=outPath, recursive=TRUE, verbose=verbose);
  }

  verbose && exit(verbose);

  res;
}) # compileAsciiDocNoweb()


############################################################################
# HISTORY:
# 2013-03-29
# o Created (from compileMarkdown.R).
############################################################################
