###########################################################################/**
# @RdocDefault rcode
# @alias rcode.RspString
# @alias rcode.RspDocument
# @alias rscript
#
# @title "Compiles an RSP document and returns the generated source code script"
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
#      The default is a file with a filename where the file extension(s)
#      (typically \code{".*.rsp"}) has been replaced by \code{".R"}
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
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an @see "RspFileProduct" if possible,
#   otherwise an @see "RspSourceCode".
# }
#
# @examples "../incl/rcode.Rex"
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
setMethodS3("rcode", "default", function(..., file=NULL, path=NULL, output=NULL, workdir=NULL, envir=parent.frame(), args="*", verbose=FALSE) {
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
    output <- RspSourceCode();

    if (inherits(file, "connection")) {
      throw("When argument 'file' is a connection, then 'output' must be specified.");
    } else if (is.character(file)) {
      # Is the input a filename or an URI?
      if (isUrl(file)) {
        # If URI, drop any URI arguments
        url <- splitUrl(file);
        filename <- basename(url$path);
        filename <- Arguments$getReadablePathname(filename, adjust="url", mustExist=FALSE);
      } else {
        filename <- basename(file);
      }

      pattern <- "((.*)[.]([^.]+)|([^.]+))[.]([^.]+)$";
      outputF <- gsub(pattern, "\\1.R", filename, ignore.case=TRUE);
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
    output <- stdout();
  } else if (inherits(output, "RspSourceCode")) {
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

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "rcode() for default");

  if (is.null(file)) {
    s <- RspString(...);
  } else {
    verbose && cat(verbose, "Input file: ", file);
    s <- .readText(file);
    s <- RspString(s, source=file, ...);
    s <- setMetadata(s, name="source", value=file);
  }
  verbose && cat(verbose, "Length of RSP string: ", nchar(s));

  res <- rcode(s, output=output, workdir=workdir, envir=envir, args=args, verbose=verbose);

  verbose && exit(verbose);

  res;
}) # rcode()


setMethodS3("rcode", "RspString", function(object, output=NULL, workdir=NULL, envir=parent.frame(), args="*", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'args':
  args <- cmdArgs(args=args);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "rcode() for ", class(object)[1L]);

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
    enter(verbose, "Parse RSP string to RSP document");
    cat(verbose, "Parse environment: ", getName(envir));
    if (length(names) > 0L) {
      ll <- subset(ll(envir=envir), member %in% names);
      print(verbose, ll);
    }
  }
  expr <- parse(object, envir=envir, ..., verbose=verbose);
  verbose && print(verbose, expr);
  verbose && exit(verbose);

  res <- rcode(expr, output=output, workdir=workdir, envir=envir, args=NULL, ..., verbose=verbose);

  verbose && exit(verbose);

  res;
}) # rcode()


setMethodS3("rcode", "RspDocument", function(object, output=NULL, workdir=NULL, envir=parent.frame(), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'workdir':
  if (is.null(workdir)) {
    workdir <- ".";
    if (inherits(output, "RspSourceCode")) {
    } else if (isAbsolutePath(output)) {
      workdir <- getParent(output);
    }
  }

  # Argument 'output':
  if (!is.null(output)) {
    if (inherits(output, "connection")) {
    } else if (inherits(output, "RspSourceCode")) {
    } else {
      withoutGString({
        output <- Arguments$getWritablePathname(output, path=workdir);
      })
      output <- getAbsolutePath(output);
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "rcode() for ", class(object)[1L]);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Coerce
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Coerce RSP document to source code");
  language <- getAttribute(object, "language", default="R");
  language <- capitalize(tolower(language));
  className <- sprintf("Rsp%sSourceCodeFactory", language);
  ns <- getNamespace("R.rsp");
  clazz <- Class$forName(className, envir=ns);
  factory <- newInstance(clazz);
  verbose && cat(verbose, "Language: ", getLanguage(factory));
  code <- toSourceCode(factory, object, ..., verbose=verbose);

  if (verbose) {
    enter(verbose, "Generated source code:");
    codeS <- c(head(code, n=10L), "", "[...]", "", tail(code, n=10L));
    cat(verbose, codeS);
    exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return as RspSourceCode, write to file, or ...?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (inherits(output, "RspSourceCode")) {
    # Return as RspSourceCode
    output <- code;
  } else if (!is.null(output)) {
    # Write to file
    verbose && enter(verbose, "Writing to output");
    writeLines(code, con=output);
    verbose && exit(verbose);

    output <- RspFileProduct(output, type=getType(code), metadata=getMetadata(code, local=TRUE), mustExist=FALSE);
  }

  verbose && exit(verbose);

  output;
}) # rcode()

## BACKWARD COMPATIBILITY:
setMethodS3("rscript", "default", function(...) {
  .Deprecated(new="rcode")
  rcode(...)
}, deprecated=TRUE)



##############################################################################
# HISTORY:
# 2015-05-12
# o Renamed rscript() to rcode().  rscript() is not deprecated.
# 2014-05-27
# o Now rscript(file) writes to file by default.
# o Now rscript() adds metadata attributes.
# o Added arguments 'output' and 'workdir' to rscript().
# 2014-01-26
# o CLEANUP: Now R.oo::ll() is only called if 'verbose' is enabled, because
#   calling ll() still triggers attachment of R.oo as of R.oo (>= 1.17.0).
# 2013-03-14
# o Created from rstring.R.
##############################################################################
