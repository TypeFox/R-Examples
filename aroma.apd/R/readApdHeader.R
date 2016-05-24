#########################################################################/**
# @RdocDefault readApdHeader
#
# @title "Reads the header of an Affymetrix probe data (APD) file"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{filename}{The filename of the APD file.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{.checkArgs}{If @TRUE, arguments are validated, otherwise not.}
# }
#
# \value{
#   A named @list.
#   The @numeric element \code{nbrOfProbes} is the number of probe values
#   available in the APD file.
#   The optional @character element \code{name} specifies the name of
#   the APD vector.
#   The optional @character element \code{chipType} specifies the
#   chip type, cf. the same field in @see "affxparser::readCelHeader".
#   The optional @character element \code{maptype} specifies the type of
#   probe-index map for this APD file.  Its value can be used to find
#   the mapping file, see @see "findApdMap" and @see "readApdMap".
#   All other fields are optional and @character values.
# }
#
# \details{
#   The file format of an APD file is identical to the file format of an
#   @see "R.huge::FileVector".  Most elements of the APD header are stored
#   in the \code{comment} @character string of the file vector's header.
#   The APD header \code{nbrOfProbes} is identical to the length of the
#   file vector, and is \emph{not} stored in the above comment string.
# }
#
# @author
#
# \examples{\dontrun{#See ?createApd for an example.}}
#
# \seealso{
#   @see "readApd".
# }
#
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
setMethodS3("readApdHeader", "default", function(filename, ..., verbose=FALSE, .checkArgs=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathname':
  if (inherits(filename, "FileVector")) {
    apd <- filename;
  } else {
    if (.checkArgs) {
      filename <- Arguments$getReadablePathname(filename, mustExists=TRUE);
    }
    apd <- FileVector(filename=filename);
    on.exit(close(apd));
  }

  if (.checkArgs) {
    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
  }

  verbose && enter(verbose, "Reading APD header");

  comment <- getComments(apd);

  header <- unlist(strsplit(comment, split="\n"), use.names=FALSE);
  header <- strsplit(header, split="=");

  keys <- unlist(lapply(header, FUN=function(s) trim(s[1])), use.names=FALSE);
  header <- lapply(header, FUN=function(s) {
    values <- trim(s[2]);
    if (!is.na(values) && regexpr(";", values) != -1) {
      values <- unlist(strsplit(values, split=";"));
      values <- strsplit(values, split=":");
      names <- unlist(lapply(values, FUN=function(s) trim(s[1])), use.names=FALSE);
      values <- unlist(lapply(values, FUN=function(s) trim(s[2])));
      names(values) <- names;
    }
    values;
  });
  names(header) <- keys;

  header[["nbrOfProbes"]] <- length(apd);

  verbose && exit(verbose);

  header;
})


############################################################################
# HISTORY:
# 2006-04-08
# o Removed argument 'path'.  Renamed argument 'pathname' to 'filename'
#   as in readApd().
# 2006-04-02
# o BUG FIX: Forgot to close FileVector on exit.
# 2006-03-15
# o BUG FIX: Forgot to close the FileVector.
# 2006-03-14
# o Added virtual header 'nbrOfProbes' for conveniency.
# 2006-02-27
# o Created by HB.
############################################################################
