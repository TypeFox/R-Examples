#########################################################################/**
# @RdocDefault writeApd
#
# @title "Writes an APD probe data file"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{filename}{The filename of the APD file.}
#   \item{data}{A @numeric @vector of elements to be written.}
#   \item{...}{Arguments passed to @see "createApd", e.g. \code{chipType},
#      \code{mapType} etc.}
#   \item{writeMap}{A @vector of indicies used to change the order how
#      data elements are \emph{written} (by default).}
# }
#
# \value{
#   Returns (invisibly) the pathname to the created file.
# }
#
# @author
#
# \seealso{
#   To create an APD map file, see @see "readApdMap".
# }
#
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
setMethodS3("writeApd", "default", function(filename, data, ..., writeMap=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename':
  pathname <- Arguments$getWritablePathname(filename, mustNotExist=TRUE);

  # Argument 'data':
  data <- Arguments$getNumerics(data);

  nbrOfCells <- length(data);

  # Argument 'writeMap':
  if (!is.null(writeMap)) {
    writeMap <- Arguments$getIndices(writeMap, range=c(1,nbrOfCells));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reorder data elements according to the write map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(writeMap)) {
    data <- data[writeMap];
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write APD file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write to a temporary pathname
  pathnameT <- pushTemporaryFile(pathname);

  # Create the APD file
  createApd(pathnameT, nbrOfCells=nbrOfCells, ...);

  # Write the values
  apd <- FileVector(pathnameT);
  # Close on exit, in case of error/interrupt.
  on.exit(if (!is.null(apd)) close(apd));

  writeAllValues(apd, data);

  # Not needed anymore
  data <- NULL;

  # Close
  close(apd); apd <- NULL;

  # Rename temporary pathname
  pathname <- popTemporaryFile(pathnameT);

  invisible(pathname);
})


############################################################################
# HISTORY:
# 2014-01-05
# o ROBUSTNESS: Now writeApd() is atomic, which is achieved by first
#   writing to a temporary file which is then renamed upon completion.
# 2009-05-16
# o Updated writeApd() to coerce argument 'writeMap' to integer indices.
#   Before it used to coerce to doubles (before updating R.utils).
# 2006-03-28
# o Renamed argument 'map' to 'writeMap'.  Removed argument 'isWriteMap'.
# 2006-03-14
# o Created.
############################################################################
