###########################################################################/**
# @RdocDefault findAsciiDoc
#
# @title "Locates the asciidoc executable"
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
#  The 'asciidoc' executable is searched for as follows:
#  \enumerate{
#   \item \code{Sys.which("asciidoc")}
#  }
# }
#
# @author
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("findAsciiDoc", "default", function(mustExist=TRUE, ..., verbose=FALSE) {
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

  command <- "asciidoc";

  verbose && enter(verbose, "Locating external software");
  verbose && cat(verbose, "Command: ", command);

  bin <- Sys.getenv("R_ASCIIDOC")
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
      ver <- trim(gsub("asciidoc", "", res))
      ver <- numeric_version(ver)
      attr(bin, "version") <- ver
    }
  }

  verbose && exit(verbose);

  bin;
}) # findAsciiDoc()


############################################################################
# HISTORY:
# 2013-04-01
# o Now findAsciiDoc() only sets 'version' attribute if executable exists.
# 2013-03-29
# o Created.
############################################################################
