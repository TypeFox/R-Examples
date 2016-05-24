###########################################################################/**
# @RdocDefault rcat
# @alias rcat.RspString
# @alias rcat.RspDocument
# @alias rcat.RspRSourceCode
# @alias rcat.function
# @alias rcat.expression
# @alias rsource
# @alias rsource.default
# @alias rsource.RspString
# @alias rsource.RspDocument
# @alias rsource.RspRSourceCode
# @alias rsource.function
# @alias rsource.expression
#
# @title "Evaluates an RSP string and outputs the generated string"
#
# \description{
#  @get "title".
# }
#
# \usage{
#  @usage rcat,default
#  @usage rsource,default
# }
#
# \arguments{
#   \item{...}{A @character string with RSP markup.}
#   \item{file, path}{Alternatively, a file, a URL or a @connection from
#      with the strings are read.
#      If a file, the \code{path} is prepended to the file, iff given.}
#   \item{envir}{The @environment in which the RSP string is
#     preprocessed and evaluated.}
#   \item{args}{A named @list of arguments assigned to the environment
#     in which the RSP string is parsed and evaluated.
#     See @see "R.utils::cmdArgs".}
#   \item{output}{A @connection, or a pathname where to direct the output.
#               If \code{""}, the output is sent to the standard output.}
#   \item{buffered}{If @TRUE, and \code{output=""}, then the RSP output is
#     outputted as soon as possible, if possible.}
#   \item{append}{Only applied if \code{output} specifies a pathname;
#     If @TRUE, then the output is appended to the file, otherwise
#     the files content is overwritten.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns (invisibly) the outputted @see "RspStringProduct".
# }
#
# \section{Processing RSP strings from the command line}{
#   Using @see "Rscript" and \code{rcat()}, it is possible to process
#   an RSP string and output the result from the command line.  For example,
#
#   \code{Rscript -e "R.rsp::rcat('A random integer in [1,<\%=K\%>]: <\%=sample(1:K, size=1)\%>')" --args --K=50}
#
#   parses and evaluates the RSP string and outputs the result to
#   standard output.
# }
#
# \section{rsource()}{
#   The \code{rsource(file, ...)} is a convenient wrapper
#   for \code{rcat(file=file, ..., output="", buffered=FALSE)}.
#   As an analogue, \code{rsource()} is to an RSP file what
#   \code{source()} is to an R script file.
# }
#
# @examples "../incl/rcat.Rex"
#
# @author
#
# \seealso{
#  To store the output in a string (instead of displaying it), see
#  @see "rstring".
#  For evaluating and postprocessing an RSP document and
#  writing the output to a file, see @see "rfile".
# }
#
# @keyword print
# @keyword IO
# @keyword file
#*/###########################################################################
setMethodS3("rcat", "default", function(..., file=NULL, path=NULL, envir=parent.frame(), args="*", output="", buffered=TRUE, append=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'file' & 'path':
  if (inherits(file, "connection")) {
  } else if (is.character(file)) {
    if (!is.null(path)) {
      file <- file.path(path, file);
    }
    if (!isUrl(file)) {
      file <- Arguments$getReadablePathname(file, absolute=TRUE);
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "rcat() for default");

  if (is.null(file)) {
    s <- RspString(...);
  } else {
    verbose && cat(verbose, "Input file: ", file);
    s <- .readText(file);
    s <- RspString(s, source=file, ...);
    s <- setMetadata(s, name="source", value=file);
  }
  verbose && cat(verbose, "Length of RSP string: ", nchar(s));

  res <- rcat(s, output=output, buffered=buffered, append=append, envir=envir, args=args, verbose=verbose);

  verbose && exit(verbose);

  invisible(res);
}) # rcat()


setMethodS3("rcat", "RspString", function(..., envir=parent.frame(), args="*", output="", buffered=TRUE, append=FALSE, verbose=FALSE) {
  # Argument 'buffered':
  if (!buffered) {
    isStdout <- FALSE;
    if (is.character(output) && output == "") {
      isStdout <- TRUE;
    } else if (inherits(output, "connection")) {
      ci <- summary(output);
      isStdout <- identical(ci$class, "terminal") &&
                  identical(ci$description, "stdout");
    }
    if (!isStdout) {
      throw("Argument 'buffered' must be TRUE unless 'output' directs to the standard output.");
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "rcat() for RspString");

  outputP <- ifelse(buffered, "RspStringProduct", "stdout");
  verbose && printf(verbose, "Buffered: %s\n", buffered)
  verbose && printf(verbose, "Type of output: %s\n", outputP)
  verbose && cat(verbose, "Arguments:")
  verbose && str(verbose, args)

  s <- rstring(..., envir=envir, args=args, output=outputP);

  verbose && cat(verbose, "Result:");
  verbose && str(verbose, s);

  if (!is.null(s)) {
    verbose && enter(verbose, "Outputting");
    outputT <- output;
    if (is.character(output)) {
      if (output == "")
        outputT <- "<stdout>"
    } else {
      outputT <- "<connection>"
    }
    verbose && printf(verbose, "Output: %s\n", outputT)
    cat(s, file=output, append=append);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  invisible(s);
}) # rcat()


setMethodS3("rcat", "RspDocument", rcat.RspString)
setMethodS3("rcat", "RspRSourceCode", rcat.RspString)
setMethodS3("rcat", "function", rcat.RspString)
setMethodS3("rcat", "expression", rcat.RspString)


##############################################################################
# HISTORY:
# 2014-01-02
# o CLEANUP: Now rcat() methods for several classes uses the exact
#   same function definition.  Also harmonized the ordering of arguments.
# o Added rstring(), rcat() and rfile() for expression:s too.
# 2013-07-16
# o Added rstring(), rcat() and rfile() for function:s.
# 2013-05-08
# o Explicitly added arguments 'file' & 'path' to rcat() [although they're
#   just passed as is to rstring()].
# 2013-02-20
# o Renamed argument 'file' for rcat() to 'output', cf. rfile().  This
#   automatically makes argument 'file' & 'path' work also for rcat()
#   just as it works for rstring() and rfile().
# 2013-02-13
# o Added rcat() for several RSP-related classes.
# 2013-02-11
# o Added Rdoc help.
# 2013-02-09
# o Created.
##############################################################################
