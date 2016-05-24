#########################################################################/**
# @RdocDefault updateApd
#
# @title "Updates an Affymetrix probe data (APD) file"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{filename}{The filename of the APD file.}
#   \item{indices}{A @numeric @vector of cell (probe) indices specifying
#     which cells to updated.}
#   \item{data}{A @numeric @vector of data elements to be assigned.}
#   \item{writeMap}{A @vector of indicies used to change the order how
#      data elements are \emph{written} (by default).}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{.checkArgs}{If @TRUE, arguments are checked, otherwise not.}
# }
#
# \value{
#   Returns (invisibly) the pathname of the file updated.
# }
#
# @author
#
# \examples{\dontrun{#See ?createApd for an example.}}
#
# \seealso{
#   @see "createApd" and @see "updateApd".
# }
#
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
setMethodS3("updateApd", "default", function(filename, indices=NULL, data, writeMap=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.checkArgs) {
    # Argument 'filename':
    filename <- Arguments$getReadablePathname(filename, mustExist=TRUE);
  }

  # Affymetrix probe signals are stored as floats in CEL files.
  apd <- FileVector(filename);
  on.exit(close(apd));
  nbrOfProbes <- length(apd);

  if (.checkArgs) {
    # Argument 'data':
    data <- Arguments$getNumerics(data);

    # Argument 'indices':
    if (is.null(indices)) {
      nbrOfIndices <- nbrOfProbes;
    } else {
      # An APD file has zero- or one-based indices
      indices <- Arguments$getIntegers(indices, range=c(1, nbrOfProbes));
      nbrOfIndices <- length(indices);
    }

    # Argument 'data':
    if (length(data) != nbrOfIndices) {
      throw("Length of argument 'data' does not match the number of cell indices: ", length(data), " != ", nbrOfIndices);
    }

    # Argument 'writeMap':
    if (!is.null(writeMap)) {
      writeMap <- Arguments$getIndices(writeMap, range=c(1,nbrOfProbes));
    }

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reorder indices and data elements according to the write map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(writeMap)) {
    if (is.null(indices)) {
      data <- data[writeMap];
    } else {
      windices <- writeMap[indices];
      o <- match(windices, indices);
      if (any(is.na(o))) {
        throw("Internal error: NA indices generated when remapping according write map.");
      }
      data <- data[o];
      indices <- windices;
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write data to file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  writeValues(apd, indices=indices, values=data);

  invisible(filename);
})


############################################################################
# HISTORY:
# 2009-05-16
# o Updated updateApd() to coerce argument 'writeMap' to integer indices.
#   Before it used to coerce to doubles (before updating R.utils).
# 2006-03-28
# o Removed argument 'indexOffset'.
# 2006-03-18
# o Added argument 'indexOffset' and made it one by default (as in R).
# 2006-02-27
# o Created by HB.
############################################################################
