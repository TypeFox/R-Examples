#########################################################################/**
# @RdocFunction readDcpRectangle
#
# @title "Reads a spatial subset of probe-level data from a dChip DCP file"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{filename}{The pathname of the DCP file.}
#   \item{fields}{The cell fields to be read.}
#   \item{xrange}{A @numeric @vector of length two giving the left
#     and right coordinates of the cells to be returned.}
#   \item{yrange}{A @numeric @vector of length two giving the top
#     and bottom coordinates of the cells to be returned.}
#   \item{asMatrix}{If @TRUE, the CEL data fields are returned as
#     matrices with element (1,1) corresponding to cell
#     (xrange[1],yrange[1]).}
#   \item{...}{Additional arguments passed to @see "readDcp".}
# }
#
# \value{
#   A named @list CEL structure similar to what @see "readDcp".
#   In addition, if \code{asMatrix} is @TRUE, the CEL data fields
#   are returned as matrices, otherwise not.
# }
#
# @author
#
# @examples "../incl/readDcpRectangle.Rex"
#
# \seealso{
#   The @see "readDcp" method is used internally.
#   This method was inspired by \code{readCelRectangle()} of the
#   \pkg{affxparser} package.
# }
#
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
readDcpRectangle <- function(filename, fields=c("rawIntensities", "normalizedIntensities"), xrange=c(0,Inf), yrange=c(0,Inf), ..., asMatrix=TRUE) {
  # Get the chip layout from the DCP header
  header <- readDcpHeader(filename);
  nrow <- ncol <- header$CellDim;

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

  ## Not needed anymore
  xrange <- yrange <- yy <- xx <- offsets <- NULL

  # Read DCP data
  data <- readDcp(filename, fields=fields, cells=cells, ...);

  # Rearrange each field into a matrix?
  if (asMatrix) {
    nrow <- nrow(cells);
    cells <- NULL   ## Not needed anymore

    # Dcpl-value fields
    fields <- c("rawIntensities", "normalizedIntensities");
    fields <- intersect(fields, names(data));
    for (field in fields) {
      data[[field]] <- matrix(data[[field]], ncol=nrow, byrow=FALSE);
    }
  }

  data;
} # readDcpRectangle()


############################################################################
# HISTORY:
# 2008-08-20
# o Created from affxparser::readCelRectangle.R
############################################################################
