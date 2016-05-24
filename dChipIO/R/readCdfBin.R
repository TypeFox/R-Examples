###########################################################################/**
# @RdocFunction readCdfBin
#
# @title "Reads a dChip CDF.bin file"
#
# \description{
#   @get "title".
#
#   Please note that this method is incomplete as it currently doesn't read
#   all fields.  It is only made available so that someelse can continue
#   the development.
# }
#
# @synopsis
#
# \arguments{
#   \item{con}{A @connection or a @character filename.}
#   \item{units}{An @integer @vector specifying the units to be read.
#     If @NULL, all units are read.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @list structure containing the file header and the unit data.
# }
#
# @examples "../incl/readCdfBin.Rex"
#
# @author
#
# \seealso{
#   To read only the CDF.bin file header, see @see "readCdfBinHeader".
# }
#
# @keyword "file"
# @keyword "IO"
#*/###########################################################################
readCdfBin <- function(con, units=NULL, ...) {
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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read CDF.bin header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hdr <- readCdfBinHeader(con=con);
  nbrOfUnits <- hdr$NumUnit;
  # Sanity checks
  stopifnot(nbrOfUnits >= 0)

  # Arguments 'units':
  if (is.null(units)) {
    units <- 1:nbrOfUnits;
  } else {
    units <- .argAssertRange(as.integer(units), range=c(1, nbrOfUnits));
    nbrOfUnits <- length(units);
  }

  data <- list(header=hdr);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read unit data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  UNIT_NAME_LEN <- 50;
#	char UnitName[UNIT_NAME_LEN]; // as in *.cdf file
#	short NumProbe; // number of probe pairs
#	SHORT_POINT_ARRAY *CellPos;
  CDF_UNIT_SIZE = UNIT_NAME_LEN + 2 + 4;

  firstUnit <- min(units);
  lastUnit <- max(units);
  nbrOfUnitsToRead <- as.integer(lastUnit - firstUnit + 1);

  # Skip to the first unit to be read
  if (firstUnit > 0) {
    seek(con=con, origin="current", where=(firstUnit-1)*CDF_UNIT_SIZE,
                                                            rw="read");
  }

  # Read first to last unit...
  raw <- .readRaw(con, n=nbrOfUnitsToRead*CDF_UNIT_SIZE)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Keep only units of interest
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dim(raw) <- c(CDF_UNIT_SIZE, nbrOfUnitsToRead);
  keep <- match(units, (firstUnit:lastUnit));
  raw <- raw[,keep,drop=FALSE];

  # Extract 'UnitName'
  idxs <- 1:UNIT_NAME_LEN;
  unitNames <- rep("", nbrOfUnits);
  for (idx in idxs) {
    unitNames <- paste(unitNames,
                   rawToChar(raw[idx,,drop=TRUE], multiple=TRUE), sep="");
  }
  raw <- raw[-idxs,,drop=FALSE];
  data$unitNames <- unitNames;
  unitNames <- NULL ## Not needed anymore

  # Extract 'NumProbe'
  idxs <- 1:2;
  data$numProbes <- .readUShort(raw[idxs,,drop=FALSE], n=nbrOfUnits)
  raw <- raw[-idxs,,drop=FALSE];
  # Sanity checks
  stopifnot(data$numProbes >= 0)

  data$CellPos <- .readInt(raw, n=nbrOfUnits)
  # Sanity checks
  stopifnot(data$CellPos >= 0)

  raw <- NULL ## Not needed anymore


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read cell data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Skip to the last unit
  if (lastUnit < nbrOfUnits) {
    seek(con=con, origin="current",
                where=(nbrOfUnits-lastUnit)*CDF_UNIT_SIZE, rw="read");
  }

  # To do

  data;
}


##############################################################################
# HISTORY:
# 2008-08-20
# o BUG FIX: The wrong subset of units was read when 'units' was specified.
# o Added Rdoc comments.
# 2008-02-03
# o Added reader for CDF.bin unit data. Cell data are still to be done.
# o Created.
##############################################################################
