###########################################################################/**
# @RdocDefault rsp
#
# @title "Compiles an RSP document"
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
#   \item{text}{A @character @vector of RSP code to be processed,
#      iff argument \code{filename} is not given.}
#   \item{response}{Specifies where the final output should be sent.
#      If argument \code{text} is given, then @see "base::stdout" is used.
#      Otherwise, the output defaults to that of the type-specific compiler.}
#   \item{...}{Additional arguments passed to the type-specific compiler.}
#   \item{envir}{The @environment in which the RSP document is evaluated.}
#   \item{outPath}{The output and working directory.}
#   \item{postprocess}{If @TRUE, and a postprocessing method exists for
#      the generated document type, it is postprocessed as well.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   If argument \code{response} specifies a file output, then the
#   absolute pathname of the generated file is returned.
#   If argument \code{text} is specified, then the generated string
#   is returned (invisibly).
# }
#
# @examples "../incl/rsp.Rex"
#
# @author
#
# \section{Postprocessing}{
#   For some document types, the \code{rsp()} method automatically
#   postprocesses the generated document as well.
# }
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("rsp", "default", function(filename=NULL, path=NULL, text=NULL, response=NULL, ..., envir=parent.frame(), outPath=".", postprocess=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(filename) && !is.null(text)) {
    throw("Only one of arguments 'filename' and 'text' can be specified.");
  }

  # Argument 'text':
  if (!is.null(text)) {
    text <- Arguments$getCharacter(text);
  }

  # Arguments 'filename' & 'path':
  if (!is.null(filename)) {
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
    pathname <- getAbsolutePath(pathname);
  } else {
    pathname <- NULL;
  }

  if (is.null(text) && is.null(pathname)) {
    throw("Either argument 'filename' or 'text' must be given.");
  }

  # Arguments 'outPath':
  if (is.null(outPath)) {
    outPath <- ".";
  } else {
    outPath <- Arguments$getWritablePath(outPath);
    if (is.null(outPath)) outPath <- getwd();
  }
  outPath <- getAbsolutePath(outPath);

  # Argument 'envir':
#  envir <- Arguments$getEnvironment(envir);

  # Argument 'postprocess':
  postprocess <- Arguments$getLogical(postprocess);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Compile an RSP string
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(text)) {
    verbose && enter(verbose, "Parsing and evaluating RSP string");

    if (is.null(response)) {
      response <- stdout();
    }

    # Change working directory?
    opwd <- getwd();
    on.exit({ if (!is.null(opwd)) setwd(opwd) }, add=TRUE);
    setwd(outPath);

    res <- rcat(text, output=response, envir=envir, ...);

    # Reset working directory
    if (!is.null(opwd)) {
      setwd(opwd);
      opwd <- NULL;
    }

    verbose && exit(verbose);
    return(invisible(res));
  } # if (!is.null(text))


  verbose && enter(verbose, "Parsing and evaluating RSP file");
  if (is.character(response)) {
    response <- getAbsolutePath(response);
  }

  verbose && enter(verbose, "Processing RSP file");
  verbose && cat(verbose, "Current directory: ", getwd());
  res <- rfile(pathname, output=response, workdir=outPath, envir=envir, postprocess=postprocess, ..., verbose=verbose);
  verbose && print(verbose, res);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(res);
}) # rsp()


############################################################################
# HISTORY:
# 2013-02-13
# o CLEANUP: Removed obsolete local tempfile() function.
# 2013-02-12
# o Now rsp(text=...) uses rcat() and rsp(file=...) uses rfile().
# 2013-02-08
# o Made internal rspPlain() its own function.
# 2011-11-14
# o Added argument 'envir' to rsp(..., envir=parent.frame()).
# 2011-04-16
# o BUG FIX: On R v2.12.x, rsp(text="...") would throw 'Error ...: unused
#   argument(s) (fileext = ".txt.rsp")'.  Solved by providing a patched
#   tempfile() with this feature for R v2.12.x.  Thanks Uwe Ligges for
#   spotting this.
# 2011-04-12
# o Added rsp().
# o Created.
############################################################################
