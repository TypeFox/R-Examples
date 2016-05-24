###########################################################################/**
# @RdocDefault compileSweave
#
# @title "Compiles a Sweave file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename, path}{The filename and (optional) path of the
#      Sweave document to be compiled.}
#   \item{...}{Additional arguments passed to @see "compileLaTeX".}
#   \item{outPath}{The output and working directory.}
#   \item{postprocess}{If @TRUE, and a postprocessing method exists for
#      the generated product, it is postprocessed as well.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the pathname of the generated document.
# }
#
# @author
#
# \seealso{
#   Internally, @see "utils::Sweave" is used.
# }
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("compileSweave", "default", function(filename, path=NULL, ..., outPath=".", postprocess=TRUE, verbose=FALSE) {
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

  verbose && enter(verbose, "Compiling Sweave document");
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
  verbose && cat(verbose, "Sweave pathname (absolute): ", pathname);
  verbose && printf(verbose, "Input file size: %g bytes\n", file.info(pathname)$size);
  verbose && cat(verbose, "Output and working directory: ", outPath);

  opwd <- ".";
  on.exit(setwd(opwd), add=TRUE);
  if (!is.null(outPath)) {
    opwd <- setwd(outPath);
  }

  pathname2 <- Sweave(pathname);
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
}) # compileSweave()


############################################################################
# HISTORY:
# 2013-02-18
# o Added argument 'fake' to compileSweave().
# 2012-12-06
# o Added argument 'outPath' to compileSweave(), which is also the
#   working directory.
# 2011-04-14
# o Now compileSweave() only calls compileLaTeX() if Sweave outputs
#   a file with filename extension *.tex (non-case sensitive).
# 2011-04-12
# o Added compileSweave().
# o Created.
############################################################################
