###########################################################################/**
# @RdocFunction readDcp
#
# @title "Reads a dChip DCP file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{con}{A @connection or a @character filename.}
#   \item{fields}{A @character @vector specifying the fields to be read.}
#   \item{cells}{An @integer @vector specifying the indices of the cell
#     data to be read.}
#   \item{units}{An @integer @vector specifying the indices of the unit
#     data to be read.}
#   \item{.nbrOfUnits}{A @integer specifying the number of units available
#     in the file. If @NULL, this is inferred from the file size and the
#     file header. The dChip software itself instead uses the corrsponding
#     value in the CDF.bin file, but that file is specified by the user
#     leaving room for errors.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @list structure containing the file header and the
#   requested data fields.
# }
#
# @examples "../incl/readDcp.Rex"
#
# @author
#
# \seealso{
#   To read only the DCP file header, see @see "readDcpHeader".
# }
#
# @keyword "file"
# @keyword "IO"
#*/###########################################################################
readDcp <- function(con, fields=c("rawIntensities", "normalizedIntensities", "calls", "thetas", "thetaStds", "excludes"),  cells=NULL, units=NULL, .nbrOfUnits=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readElements <- function(con, idxs, nbrOfElements, size=1, skip=FALSE, ..., drop=TRUE) {
    # Local functions
    readBin2 <- function(con, what, size, signed, n) {
      if (mode(what) == "raw") {
        n <- n * size;
        size <- 1;
      }
#print(list(what=what, size=size, signed=signed, n=n));
      readBin(con, what=what, size=size, signed=signed, n=n, endian="little")
    } # readBin2()


    if (skip) {
      seek(con, where=size*nbrOfElements, origin="current", rw="read")
      return(NULL);
    }

    # Read all elements
    values <- readBin2(con=con, size=size, ..., n=nbrOfElements);

    # Turn into a matrix?
    if (mode(values) == "raw") {
      # Sanity check
      if (length(values) %% size != 0) {
        stop("File format error/read error: The number of bytes read is not a multiple of ", size, ": ", length(values));
      }

      dim(values) <- c(size, (length(values) %/% size));

      # Keep only the subset of elements?
      if (!is.null(idxs)) {
        values <- values[,idxs,drop=FALSE];
      }

      # Drop singleton dimensions?
      if (drop)
        values <- drop(values);
    } else {
      if (!is.null(idxs))
        values <- values[idxs];
    }

    values;
  } # readElements()


  readCells <- function(con, cells, what=integer(), size=2, signed=FALSE, ...) {
    readElements(con=con, idxs=cells, nbrOfElements=nbrOfCells,
                                   what=what, size=size, signed=signed, ...);
  }


  readUnits <- function(con, units, what=double(), size=2, signed=FALSE, ...) {
    readElements(con=con, idxs=units, nbrOfElements=nbrOfUnits,
                                   what=what, size=size, signed=signed, ...);
  }

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

  # Argument 'fields':
  fields <- match.arg(fields, several.ok=TRUE);


  # Reading file header
  res <- list();
  res$header <- readDcpHeader(con=con);
  nbrOfCells <- res$header$CellDim^2;
  # Sanity checks
  stopifnot(nbrOfCells >= 0)

  # Argument '.nbrOfUnits':
  if (is.null(.nbrOfUnits)) {
    # Inferring number of units from the file size.
    conInfo <- summary(con);
    if (conInfo$class != "file") {
      stop("Cannot infer the value of '.nbrOfUnits' from the connection, because it is not a file: ", conInfo$class);
    }

    pathname <- conInfo$description;
    if (!is.element(res$header$Format, c(3,4))) {
      stop("Cannot infer the value of '.nbrOfUnits' from the file. The file format is not v3 or v4 but v", res$header$Format, ": ", pathname);
    }

    nbrOfBytes <- file.info(pathname)$size;
    fileHeaderSize <- 3028;
    nbrOfUnitBytes <- nbrOfBytes - fileHeaderSize - 2*2*nbrOfCells;
    .nbrOfUnits <- as.integer(nbrOfUnitBytes / 13);

    # Sanity check
    if (nbrOfUnitBytes %% 13 != 0) {
      stop("Internal file format assumption error: Cannot infer the value of '.nbrOfUnits' from the file. The number of inferred bytes storing unit data is not a multiple of 13: ", nbrOfUnitBytes);
    }
  } else {
    .nbrOfUnits <- .argAssertRange(as.integer(.nbrOfUnits),
                                                    range=c(1, nbrOfCells));
    if (length(.nbrOfUnits) != 1)
      stop("Argument '.nbrOfUnits' must be a single integer.");
  }
  nbrOfUnits <- .nbrOfUnits;

  # Argument 'cells':
  if (!is.null(cells)) {
    cells <- .argAssertRange(as.integer(cells), range=c(1, nbrOfCells));
  }

  # Argument 'units':
  if (!is.null(units)) {
    units <- .argAssertRange(as.integer(units), range=c(1, nbrOfCells));
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading data section
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading/skipping raw and normalized intensities
  for (field in c("rawIntensities", "normalizedIntensities")) {
    res[[field]] <- readCells(con, cells=cells, skip=(!field %in% fields), drop=TRUE);
#    str(res);
  }

  # Reading/skipping calls
  field <- "calls";
  res[[field]] <- readUnits(con, units=units, what=raw(), size=1,
                                       skip=(!field %in% fields), drop=TRUE);
#  str(res);

  # Reading/skipping probe summaries
  skip <- !any(c("thetas", "thetaStds", "excludes") %in% fields);
  raw <- readUnits(con, units=units, what=raw(), size=12, skip=skip);

  # Parse probe summaries
  if (!skip) {
    n <- ncol(raw);

    field <- "thetas";
    if (field %in% fields) {
      res[[field]] <- .readFloat(con=raw[1:4,], n=n);
    }
    raw <- raw[-(1:4),,drop=FALSE];

    field <- "thetaStds";
    if (field %in% fields) {
      res[[field]] <- .readFloat(con=raw[1:4,], n=n);
    }
    raw <- raw[-(1:4),,drop=FALSE];

    field <- "excludes";
    if (field %in% fields) {
      res[[field]] <- .readInt(raw[1:4,], n=n);
    }
  }

  raw <- NULL ## Not needed anymore

# print(seek(con, read="r"));
#print(file.info(pathname2)$size);

  res;
} # readDcp()


##############################################################################
# HISTORY:
# 2008-08-20
# o BUG FIX: The wrong subset of units was read when 'units' was specified.
# o Added Rdoc comments.
# o If '.nbrOfUnits=NULL', then it is inferred automatically from the file
#   size, if possible.
# 2008-01-30
# o Created.
##############################################################################
