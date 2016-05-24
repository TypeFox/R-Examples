###########################################################################/**
# @RdocFunction readCdfBinHeader
#
# @title "Reads the file header of a dChip CDF.bin file"
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
#   To read the CDF.bin file data, see @see "readCdfBin".
# }
#
# @keyword "file"
# @keyword "IO"
#*/###########################################################################
readCdfBinHeader <- function(con, ...) {
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

  hdr <- list();

  LINE_LEN <- 1000;
  UNIT_NAME_LEN <- 50;
  hdr$FileName <- .readString(con, n=LINE_LEN-1)
  hdr$Format <- .readByte(con)
  hdr$ChipType <- .readString(con, n=UNIT_NAME_LEN)

  # Don't know why we have to add this, but in the dChip source code file
  # 'cdf.cpp' there is the following comment:
  # "// 11/11/05, +2 since allocation is in the unit of 4 bytes"
  hdr$dummy <- .readRaw(con, n=2)

  hdr$CellDim <- .readInt(con)
  hdr$NumUnit <- .readInt(con)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Clean up the header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hdr$dummy <- NULL;

  # Sanity checks
  stopifnot(all(hdr$CellDim >= 0))
  stopifnot(all(hdr$NumUnit >= 0))

  hdr;
} # readCdfBinHeader()



##############################################################################
# HISTORY:
# 2009-02-13
# o Using .rawToString() instead of rawToChar() to avoid warnings on
#   'truncating string with embedded nul:...'.
# 2008-08-20
# o Added Rdoc comments.
# 2008-02-03
# o Created.
##############################################################################
