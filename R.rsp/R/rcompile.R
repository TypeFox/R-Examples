###########################################################################/**
# @RdocDefault rcompile
# @alias rcompile.RspString
# @alias rcompile.RspDocument
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
#   \item{...}{@character strings with RSP markup.}
#   \item{file, path}{Alternatively, a file, a URL or a @connection from
#      with the strings are read.
#      If a file, the \code{path} is prepended to the file, iff given.}
#   \item{output}{A @character string or a @connection specifying where
#      output should be directed.
#      The default is a file with a filename where the input filename
#      has been prepended by \code{compiled-} and saved
#      in the directory given by the \code{workdir} argument.}
#   \item{workdir}{The working directory to use after parsing and
#      preprocessing.
#      If argument \code{output} specifies an absolute pathname,
#      then the directory of \code{output} is used, otherwise the
#      current directory is used.}
#   \item{envir}{The @environment in which the RSP string is
#      preprocessed and evaluated.}
#   \item{args}{A named @list of arguments assigned to the environment
#     in which the RSP string is parsed and evaluated.
#     See @see "R.utils::cmdArgs".}
#   \item{until}{Specifies how far the compilation should proceed.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an @see "RspString", @see "RspDocument" or
#   an @see "RspFileProduct" (depending on argument \code{output}).
# }
#
# @author
#
# \seealso{
#  @see "rcat" and @see "rfile".
# }
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("rcompile", "default", function(..., file=NULL, path=NULL, output=NULL, workdir=NULL, envir=parent.frame(), args="*", until="*", verbose=FALSE) {
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

  # Argument 'workdir':
  if (is.null(workdir)) {
    if (isAbsolutePath(output)) {
      workdir <- getParent(output);
    } else {
      workdir <- ".";
    }
  }
  workdir <- Arguments$getWritablePath(workdir);
  if (is.null(workdir)) workdir <- ".";

  # Argument 'output':
  if (is.null(output)) {
    # Default is to return an RSP source code object
    output <- RspString();

    if (inherits(file, "connection")) {
      throw("When argument 'file' is a connection, then 'output' must be specified.");
    } else if (is.character(file)) {
      # Is the input a filename or an URI?
      if (isUrl(file)) {
        # If URI, drop any URI arguments
        url <- splitUrl(file);
        filename <- basename(url$path);
        filename <- Arguments$getReadablePathname(filename, adjust="url", mustExist=FALSE);
        filename <- basename(file);
      } else {
        filename <- basename(file);
      }

      outputF <- sprintf("compiled-%s", filename);
      withoutGString({
        output <- Arguments$getWritablePathname(outputF, path=workdir);
      })
      output <- getAbsolutePath(output);
      # Don't overwrite the input file
      if (output == file) {
        throw("Cannot process RSP file. The inferred argument 'output' is the same as argument 'file' & 'path': ", output, " == ", file);
      }
    }
  } else if (inherits(output, "connection")) {
  } else if (identical(output, "")) {
    output <- RspString();
  } else if (inherits(output, "RspString")) {
  } else if (inherits(output, "RspDocument")) {
  } else if (is.character(output)) {
    withoutGString({
      if (isAbsolutePath(output)) {
        output <- Arguments$getWritablePathname(output);
      } else {
        output <- Arguments$getWritablePathname(output, path=workdir);
        output <- getAbsolutePath(output);
      }
    })
    if (is.character(file) && (output == file)) {
      throw("Cannot process RSP file. Argument 'output' specifies the same file as argument 'file' & 'path': ", output, " == ", file);
    }
  } else {
    throw("Argument 'output' of unknown type: ", class(output)[1L]);
  }

  # Argument 'until':
##  until <- match.arg(until);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "rcompile() for default");
  verbose && cat(verbose, "Compile until: ", sQuote(until));
  verbose && cat(verbose, "Workdir: ", workdir);
  verbose && cat(verbose, "Output: ", output);

  if (is.null(file)) {
    s <- RspString(...);
  } else {
    verbose && cat(verbose, "Input: ", file);
    s <- .readText(file);
    s <- RspString(s, source=file);
    s <- setMetadata(s, name="source", value=file);
  }
  verbose && cat(verbose, "Length of RSP string: ", nchar(s));

  res <- rcompile(s, envir=envir, output=output, workdir=workdir, args=args, until=until, verbose=verbose);

  verbose && exit(verbose);

  res;
}) # rcompile()


setMethodS3("rcompile", "RspString", function(object, envir=parent.frame(), output=NULL, workdir=NULL, args="*", ..., until="*", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'workdir':
  if (is.null(workdir)) {
    workdir <- ".";
    if (inherits(output, "RspString")) {
    } else if (inherits(output, "RspDocument")) {
    } else if (isAbsolutePath(output)) {
      workdir <- getParent(output);
    }
  }

  # Argument 'output':
  if (!is.null(output)) {
    if (inherits(output, "connection")) {
    } else if (inherits(output, "RspString")) {
    } else if (inherits(output, "RspDocument")) {
    } else {
      withoutGString({
        if (isAbsolutePath(output)) {
          output <- Arguments$getWritablePathname(output);
        } else {
          output <- Arguments$getWritablePathname(output, path=workdir);
          output <- getAbsolutePath(output);
        }
      })
    }
  }

  # Argument 'args':
  args <- cmdArgs(args=args);

  # Argument 'until':
##  until <- match.arg(until);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "rcompile() for ", class(object)[1L]);
  verbose && cat(verbose, "Compile until: ", sQuote(until));

  if (length(args) > 0L) {
    verbose && enter(verbose, "Assigning RSP arguments to processing environment");
    verbose && cat(verbose, "Environment: ", getName(envir));

    verbose && cat(verbose, "RSP arguments:");
    verbose && str(verbose, args);

    # Assign arguments to the parse/evaluation environment
    names <- attachLocally(args, envir=envir);
    if (verbose) {
      if (length(names) > 0L) {
        printf(verbose, "Variables assigned: [%d] %s\n", length(names), hpaste(names));
        member <- NULL; rm(list="member"); # To please R CMD check
        ll <- subset(ll(envir=envir), member %in% names);
        print(verbose, ll);
      }
    }
    verbose && exit(verbose);
  } else {
    names <- NULL;
  }

  if (verbose) {
    enter(verbose, "Parse RSP string");
    cat(verbose, "Parse environment: ", getName(envir));
    if (length(names) > 0L) {
      ll <- subset(ll(envir=envir), member %in% names);
      print(verbose, ll);
    }
  }

  # Class to parse to
  as <- if (inherits(output, "RspDocument")) "RspDocument" else "RspString";
  res <- parse(object, envir=envir, ..., until=until, as=as, verbose=verbose);
  verbose && print(verbose, res);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return as RspSourceCode, write to file, or ...?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (inherits(output, "RspDocument")) {
  } else if (inherits(output, "RspString")) {
  } else if (!is.null(output)) {
    # Write to file
    verbose && enter(verbose, "Writing to output");
    cat(res, file=output);
    verbose && exit(verbose);
    res <- RspFileProduct(output, type=getType(res), metadata=getMetadata(res, local=TRUE), mustExist=FALSE);
  }

  verbose && exit(verbose);

  res;
}) # rcompile()



setMethodS3("rcompile", "RspDocument", function(object, envir=parent.frame(), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "rcompile() for ", class(object)[1L]);

  verbose && enter(verbose, "Coercing RSP document to RSP string");
  s <- asRspString(object);
  verbose && exit(verbose);

  res <- rcompile(s, ..., envir=envir, verbose=verbose);

  verbose && exit(verbose);

  res;
}) # rcompile()



##############################################################################
# HISTORY:
# 2014-05-30
# o Now rcompile() can write to file and if the input is a file, then
#   it does so by default.
# o Added arguments 'output' and 'workdir' to rcompile() and dropped
#   argument 'as'.
# 2014-01-26
# o CLEANUP: Now R.oo::ll() is only called if 'verbose' is enabled, because
#   calling ll() still triggers attachment of R.oo as of R.oo (>= 1.17.0).
# 2013-03-10
# o Added rcompile().
# 2013-02-23
# o Now rstring() captures standard output such that all user output to
#   stdout will be part of the output document in the order they occur.
# 2013-02-16
# o Now rstring() only takes a character vector; it no longer c(...) the
#   '...' arguments.
# o Added argument 'args' to rstring() for RspString and RspRSourceCode.
# o Added argument 'envir' to rstring() for RspString.
# 2013-02-13
# o Added argument 'file' to rstring().
# 2013-02-11
# o Added Rdoc help.
# 2013-02-09
# o Created.
##############################################################################
