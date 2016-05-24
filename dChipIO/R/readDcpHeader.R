###########################################################################/**
# @RdocFunction readDcpHeader
#
# @title "Reads the file header of a dChip DCP file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{con}{A @connection or a @character filename.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @list structure containing the file header.
# }
#
# @author
#
# \seealso{
#   To read also the DCP file data, see @see "readDcp".
# }
#
# @keyword "file"
# @keyword "IO"
#*/###########################################################################
readDcpHeader <- function(con, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'con':
  if (is.character(con)) {
    pathname <- con;
    if (!file.exists(pathname)) {
      stop("File not found: ", pathname);
    }
    con <- file(con, open="rb");
    on.exit({
      if (!is.null(con))
        close(con);
      con <- NULL;
    });
  }

  if (!inherits(con, "connection")) {
    stop("Argument 'con' must be either a connection or a pathname: ",
                                                            mode(con));
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read the file header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hdr <- list();

  hdr$Header <- .readString(con, n=1000)
  hdr$Format <- .readByte(con, n=1)

  # Assert that the file format is correct/the expected on
  if (!is.element(hdr$Format, 3:4)) {
    stop("File format error: The DCP header format ('Format') is not v3 or v4: ", hdr$Format);
  }

  hdr$Normalized <- as.logical(.readByte(con, n=1))
  hdr$ThetaValid <- as.logical(.readByte(con, n=1))

  # Don't know why we have to add this one
  hdr$dummy <- .readRaw(con, n=1)
# print(seek(con, read="r"));

  hdr$Median <- .readInt(con, n=1)
  hdr$MaxInten <- .readInt(con, n=1)
  hdr$CellDim <- .readInt(con, n=1)
# print(seek(con, read="r"));

  hdr$DatFile <- .readString(con, n=1000)
# print(seek(con, read="r"));
  hdr$BaselineFile <- .readString(con, n=1000)
# print(seek(con, read="r"));

  hdr$ArrayOutlierPct <- .readFloat(con, n=1)
  hdr$SingleOutlierPct <- .readFloat(con, n=1)
  hdr$PresencePct <- .readFloat(con, n=1)
# print(seek(con, read="r"));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Clean up the header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hdr$dummy <- NULL;
  hdr$Header <- hdr$Header
  hdr$DatFile <- hdr$DatFile
  hdr$BaselineFile <- hdr$BaselineFile

  # Sanity checks
  stopifnot(hdr$CellDim >= 0)

  hdr;
} # readDcpHeader()


##############################################################################
# HISTORY:
# 2009-02-13
# o Using .rawToString() instead of rawToChar() to avoid warnings on
#   'truncating string with embedded nul:...'.
# 2008-xx-xx
# o Added Rdoc comments.
# o Added argument 'clean'.
# 2008-01-30
# o Created.
##############################################################################

