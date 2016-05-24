###########################################################################/**
# @RdocDefault findPandoc
#
# @title "Locates the pandoc executable"
#
# \description{
#  @get "title" on the current system.
# }
#
# @synopsis
#
# \arguments{
#   \item{mustExist}{If @TRUE, an exception is thrown if the executable
#      could not be located.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the pathname to the executable, or @NULL if not found.
# }
#
# \details{
#  The executable is searched for as follows:
#  \enumerate{
#   \item \code{Sys.which("pandoc")}
#  }
# }
#
# @author
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("findPandoc", "default", function(mustExist=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'mustExist':
  mustExist <- Arguments$getLogical(mustExist);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  command <- "pandoc";

  verbose && enter(verbose, "Locating external software");
  verbose && cat(verbose, "Command: ", command);

  bin <- Sys.getenv("R_PANDOC")
  if (identical(bin, "")) bin <- Sys.getenv("RSTUDIO_PANDOC")
  if (identical(bin, "")) bin <- Sys.which(command)
  if (identical(bin, "")) bin <- NULL
  if (!isFile(bin)) bin <- NULL

  verbose && cat(verbose, "Located pathname: ", bin);

  if (mustExist && !isFile(bin)) {
    throw(sprintf("Failed to located external executable: '%s'", command));
  }

  # Validate by retrieving version information
  if (isFile(bin)) {
    res <- tryCatch({
      system2(bin, args="--version", stdout=TRUE)
    }, error = function(ex) {
      NULL
    })

    if (!is.null(res)) {
      pattern <- "pandoc.* ([0-9.-]+).*"
      ver <- grep(pattern, res, value=TRUE)
      ver <- gsub(pattern, "\\1", ver)
      ver <- numeric_version(ver)
      attr(bin, "version") <- ver
    }
  }

  verbose && exit(verbose);

  bin;
}) # findPandoc()


############################################################################
# HISTORY:
# 2013-04-01
# o Created from findAsciiDoc.R.
############################################################################
