###########################################################################/**
# @RdocDefault rclean
# @alias rclean.RspString
# @alias rclean.RspDocument
#
# @title "Compiles an RSP document into a preprocessed and validated RSP document"
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
#   \item{envir}{The @environment in which the RSP string is preprocessed.}
#   \item{args}{A named @list of arguments assigned to the environment
#     in which the RSP document is parsed.
#     See @see "R.utils::cmdArgs".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an @see "RspString".
# }
#
# @examples "../incl/rclean.Rex"
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
setMethodS3("rclean", "default", function(..., file=NULL, path=NULL, envir=parent.frame(), args="*", verbose=FALSE) {
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


  verbose && enter(verbose, "rclean() for default");

  if (is.null(file)) {
    s <- RspString(...);
  } else {
    verbose && cat(verbose, "Input file: ", file);
    s <- .readText(file);
    s <- RspString(s, source=file, ...);
    s <- setMetadata(s, name="source", value=file);
  }
  verbose && cat(verbose, "Length of RSP string: ", nchar(s));

  res <- rclean(s, envir=envir, args=args, verbose=verbose);

  verbose && exit(verbose);

  res;
}) # rclean()


setMethodS3("rclean", "RspString", function(object, envir=parent.frame(), args="*", ..., verbose=FALSE) {
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

  verbose && enter(verbose, "rclean() for ", class(object)[1L]);

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

  res <- rclean(expr, envir=envir, args=NULL, ..., verbose=verbose);

  verbose && exit(verbose);

  res;
}) # rclean()


setMethodS3("rclean", "RspDocument", function(object, envir=parent.frame(), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "rclean() for ", class(object)[1L]);

  verbose && enter(verbose, "Coerce RSP document to RSP string");
  s <- asRspString(object);
  verbose && exit(verbose);

  verbose && exit(verbose);

  s;
}) # rclean()


##############################################################################
# HISTORY:
# 2014-01-26
# o CLEANUP: Now R.oo::ll() is only called if 'verbose' is enabled, because
#   calling ll() still triggers attachment of R.oo as of R.oo (>= 1.17.0).
# 2013-03-14
# o Created from rscript.R.
##############################################################################
