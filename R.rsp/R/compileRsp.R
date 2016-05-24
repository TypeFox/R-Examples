###########################################################################/**
# @RdocDefault compileRsp
#
# @title "Compiles an RSP file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename, path}{The filename and (optional) path of the
#      RSP document to be compiled.}
#   \item{...}{Additional arguments passed to @see "R.rsp::rfile".}
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
setMethodS3("compileRsp", "default", function(filename, path=NULL, ..., outPath=".", envir=parent.frame(), verbose=FALSE) {

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

  verbose && enter(verbose, "Compiling RSP document");

  # A local file?
  if (!isUrl(pathname)) {
    pathname <- getAbsolutePath(pathname);
    verbose && cat(verbose, "RSP pathname (absolute): ", pathname);
    verbose && printf(verbose, "Input file size: %g bytes\n", file.info(pathname)$size);
  }

  verbose && cat(verbose, "Output and working directory: ", outPath);

  opwd <- ".";
  on.exit(setwd(opwd), add=TRUE);
  if (!is.null(outPath)) {
    opwd <- setwd(outPath);
  }

  pathname2 <- rfile(pathname, ..., envir=envir);
  setwd(opwd); opwd <- ".";

  res <- pathname2;
  verbose && print(verbose, res);

  verbose && exit(verbose);

  res;
}) # compileRsp()


###########################################################################
# HISTORY:
# 2013-12-14
# o BUG FIX: compileRsp() did not acknowledge argument 'postprocess'.
# 2013-12-13
# o BUG FIX: compileRsp(..., postprocess=TRUE) would throw 'Error in
#   UseMethod("process"): no applicable method for 'process' applied to
#   an object of class "character"'.
# 2013-03-29
# o Created.
###########################################################################
