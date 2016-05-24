#########################################################################/**
# @RdocDefault readApdRectangle
#
# @title "Reads a spatial subset of probe-level data from Affymetrix APD files"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{filename}{The pathname of the APD file.}
#   \item{xrange}{A @numeric @vector of length two giving the left
#     and right coordinates of the cells to be returned.}
#   \item{yrange}{A @numeric @vector of length two giving the top
#     and bottom coordinates of the cells to be returned.}
#   \item{...}{Additional arguments passed to @see "readApd".}
#   \item{asMatrix}{If @TRUE, the APD data fields are returned as
#     matrices with element (1,1) corresponding to cell
#     (xrange[1],yrange[1]).}
# }
#
# \value{
#   A named @list APD structure similar to what @see "readApd".
#   In addition, if \code{asMatrix} is @TRUE, the APD data fields
#   are returned as matrices, otherwise not.
# }
#
# @author
#
# @examples "../incl/readApdRectangle.Rex"
#
# \seealso{
#   The @see "readApd" method is used internally.
# }
#
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
setMethodS3("readApdRectangle", "default", function(filename, xrange=c(0,Inf), yrange=c(0,Inf), ..., asMatrix=TRUE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser");

  # Get the chip layout from the APD header
  header <- readApdHeader(filename);
  chipType <- header$chipType;
  cdfFile <- affxparser::findCdf(chipType);
  if (is.null(cdfFile)) {
    throw("Could not find the CDF file for the APD's chip type: ", chipType);
  }
  header <- affxparser::readCdfHeader(cdfFile);
  nrow <- header$rows;
  ncol <- header$cols;

  xrange[1] <- max(min(xrange[1], ncol-1),0);
  xrange[2] <- max(min(xrange[2], ncol-1),0);
  yrange[1] <- max(min(yrange[1], nrow-1),0);
  yrange[2] <- max(min(yrange[2], nrow-1),0);

  yy <- yrange[1]:yrange[2];
  offsets <- yy * ncol + xrange[1];
  xrange <- xrange - xrange[1];
  xx <- xrange[1]:xrange[2];

  cells <- matrix(offsets, ncol=length(yy), nrow=length(xx), byrow=TRUE);
  # Cell indices are one-based in R
  cells <- cells + xx + 1;
  # Not needed anymore
  xrange <- yrange <- yy <- xx <- offsets <- NULL;

  # Read APD data
  apd <- readApd(filename, indices=cells, ...);

  # Rearrange each field into a matrix?
  if (asMatrix) {
    nrow <- nrow(cells);
    cells <- NULL; # Not needed anymore

    # Cell-value fields
    fields <- c("x", "y", "intensities", "stdvs", "pixels");
    fields <- intersect(fields, names(apd));
    for (field in fields) {
      apd[[field]] <- matrix(apd[[field]], ncol=nrow, byrow=TRUE);
    }
  }

  apd;
}) # readApdRectangle()

############################################################################
# HISTORY:
# 2006-03-30
# o Created from readCelRectangle.R in affxparser.
############################################################################
