###########################################################################/**
# @RdocClass CnagCfhFile
#
# @title "The CnagCfhFile class"
#
# \description{
#  @classhierarchy
#
#  A CnagCfhFile object represents a single CNAG CFH file.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "aroma.core::AromaMicroarrayDataFile".}
#   \item{cdf}{An optional @see "AffymetrixCdfFile"}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#   An object of this class is typically part of an @see "CnagCfhSet".
# }
#*/###########################################################################
setConstructorS3("CnagCfhFile", function(..., cdf=NULL) {
  this <- extend(AffymetrixFile(...), "CnagCfhFile",
    "cached:.header" = NULL,
    .cdf = NULL
  )

  if (!is.null(cdf))
    setCdf(this, cdf);

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("clone", "CnagCfhFile", function(this, ..., verbose=TRUE) {
  # Clone itself (and clear the cached fields)
  object <- NextMethod("clone", clear=TRUE);

  # Clone the CDF here.
  if (!is.null(object$.cdf))
    object$.cdf <- clone(object$.cdf);

  object;
}, protected=TRUE)


setMethodS3("as.character", "CnagCfhFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(this))));
  s <- c(s, sprintf("Timestamp: %s", as.character(getTimestamp(this))));
  s;
}, protected=TRUE)


setMethodS3("getExtensionPattern", "CnagCfhFile", function(static, ...) {
  "[.](cfh|CFH)$";
}, static=TRUE, protected=TRUE)




setMethodS3("getIdentifier", "CnagCfhFile", function(this, ..., force=FALSE) {
  id <- this$.identifier;
  if (force || is.null(id)) {
    id <- getChecksum(this);
    this$.identifier <- id;
  }
  id;
}, private=TRUE)


###########################################################################/**
# @RdocMethod fromFile
#
# @title "Defines an CnagCfhFile object from a CFH file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{filename}{The filename of to the file.}
#  \item{path}{The path to the file.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "CnagCfhFile" object.
#  If the file is not found or if it is of the wrong file format, an
#  error is thrown.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("fromFile", "CnagCfhFile", function(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (.checkArgs) {
    # Argument 'filename' and 'path':
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
  } else {
    pathname <- filename;
  }


  readString <- function(con, ...) {
    len <- readBin(con, what=integer(), size=1, n=1);
    s <- readBin(con, what=raw(), n=len);
    s <- rawToChar(s);
    s;
  } # readString();

  con <- file(pathname, open="rb");
  on.exit(close(con));

  magic <- readString(con);
  if (!identical(magic, "1.001")) {
    throw("Could not read CPH file. File format error: ", pathname);
  }

  # Create a new instance of the same class
  newInstance(static, pathname);
}, static=TRUE, protected=TRUE)



###########################################################################/**
# @RdocMethod getCdf
#
# @title "Gets the CDF structure for this CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCdfFile" object.
# }
#
# \seealso{
#   @seemethod "setCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCdf", "CnagCfhFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    chipType <- getHeader(this)$chipType;
    cdf <- AffymetrixCdfFile$byChipType(chipType);
    this$.cdf <- cdf;
  }
  cdf;
})


###########################################################################/**
# @RdocMethod setCdf
#
# @title "Sets the CDF structure for this CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{cdf}{An @see "AffymetrixCdfFile" object.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# \seealso{
#   @seemethod "getCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("setCdf", "CnagCfhFile", function(this, cdf, ..., .checkArgs=TRUE) {
  if (.checkArgs) {
    # Argument 'cdf':
    cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

    # Assure that the CDF is compatible with the CEL file
#    if (nbrOfCells(cdf) != nbrOfCells(this)) {
#      throw("Cannot set CDF. The specified CDF structure ('", getChipType(cdf), "') is not compatible with the chip type ('", getChipType(this), "') of the CEL file. The number of cells do not match: ", nbrOfCells(cdf), " != ", nbrOfCells(this));
#    }

    # Nothing to do?
#    oldCdf <- getCdf(this);
#    if (equals(cdf, oldCdf))
#      return(invisible(this));
  }

  # Have to clear the cache
  clearCache(this);

  this$.cdf <- cdf;

  invisible(this);
})



###########################################################################/**
# @RdocMethod getHeader
#
# @title "Gets the header of the CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @list structure as returned by @see "affxparser::readCelHeader".
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getHeader", "CnagCfhFile", function(this, ...) {
  header <- this$.header;
  if (is.null(header)) {
    header <- readCfhHeader(getPathname(this), ...);
    this$.header <- header;
  }
  header;
}, private=TRUE)




setMethodS3("getHeaderLength", "CnagCfhFile", function(this, ...) {
  nbrOfBytes <- getFileSize(this)
  nbrOfSnps <- nbrOfSnps(this);
  dataOffset <- nbrOfBytes %% nbrOfSnps;
  bytesPerSnp <- nbrOfBytes %/% nbrOfSnps;
  dataOffset;
}, private=TRUE);


setMethodS3("nbrOfSnps", "CnagCfhFile", function(this, ...) {
  nbrOfSnps(getCdf(this));
})


setMethodS3("getTimestamp", "CnagCfhFile", function(this, format="%m/%d/%y %H:%M:%S", ...) {
  # Argument 'format':
  format <- Arguments$getCharacter(format, length=c(1,1));

  # Get the CEL v3 header of the CEL header
  header <- getHeader(this);

  # Get the DAT header
  timeStamp <- header$timeStamp;
  if (is.null(timeStamp))
    timeStamp <- NA;

  timeStamp;
}, private=TRUE)



setMethodS3("nbrOfCells", "CnagCfhFile", function(this, ...) {
  as.integer(NA);
})


###########################################################################/**
# @RdocMethod getChipType
#
# @title "Gets the chip type for this CEL file"
#
# \description{
#  @get "title" \emph{according} to the @see "AffymetrixCdfFile" object.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getChipType", "CnagCfhFile", function(this, ...) {
  getChipType(getCdf(this));
}, private=TRUE)



###########################################################################/**
# @RdocMethod readUnits
#
# @title "Reads CEL data unit by unit"
#
# \description{
#  @get "title" for all or a subset of units (probeset).
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be read. If @NULL, all units are read.}
#   \item{cdf}{An alternative CDF structure to be used.  This overrides
#     the \code{units} arguments.}
#   \item{...}{Arguments passed to \code{getUnits()} of the
#     @see "AffymetrixCdfFile" class (if \code{cdf} was not specified),
#     but also to the @see "affxparser::readCelUnits" methods.}
# }
#
# \value{
#  Returns the @list structure that @see "affxparser::readCelUnits" returns.
# }
#
# \section{Caching}{
#   CEL data is neither cached in memory nor on file by this method.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("readUnits", "CnagCfhFile", function(this, units=NULL, ..., verbose=FALSE) {
  # Argument 'units':
  cdf <- getCdf(this);
  if (!is.null(units)) {
    # A zero-offset index? /HB 2010-01-01
    units <- Arguments$getIndices(units, range=c(0, nbrOfUnits(cdf)));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  offset <- getHeaderLength(this);
  nbrOfSnps <- nbrOfSnps(this);
  bytesPerSnp <- 13;
  nbrOfBytes <- bytesPerSnp*nbrOfSnps;
  map <- offset + matrix(1:nbrOfBytes, nrow=bytesPerSnp);

  # Read subset of units
  if (!is.null(units)) {
    snps <- indexOf(cdf, "SNP_");
    # Map unit indices to CNAG SNP indices
    idxs <- match(units, snps);
    idxs <- na.omit(idxs);
    # Subset of SNPs to read
    map <- map[,idxs,drop=FALSE];
  }

  # Read from file
  pathname <- getPathname(this);
  raw <- readBin(pathname, what=raw(), n=nbrOfBytes);

  # Bytes 1:8 contains (thetaA,thetaB) as floats
  map <- map[1:8,,drop=FALSE];
  theta <- readBin(raw[map], what=double(), size=4, endian="little",
                                                           n=2*ncol(map));
  theta <- matrix(theta, ncol=2, byrow=TRUE);

  # Remap according to units
  if (!is.null(units)) {
    idxs <- match(snps[idxs], units);
    naValue <- as.double(NA);
    tmp <- matrix(naValue, ncol=2, nrow=length(units));
    tmp[idxs,] <- theta;
    theta <- tmp;
  }

  colnames(theta) <- c("A", "B");

  theta;
}, private=TRUE)


setMethodS3("range", "CnagCfhFile", function(this, ..., na.rm=TRUE) {
  x <- readDataFrame(this, ...);
  range(x, na.rm=na.rm);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2013-05-22
# o CLEANUP: Now getIdentifier() for CnagCfhFile utilizes getChecksum()
#   for the GenericDataFile class.
# 2012-11-20
# o CLEANUP: Deprecated "[" and "[[", because they should be used to
#   subset files and not units.
# 2011-02-24
# o BACKWARD COMPATILITY: getIdentifier() for CnagCfhFile generates
#   a checksum id based on the relative pathname.  For now, we simply
#   drop any tags from the root path such that it is compatible with
#   earlier version of aroma.*.
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: readUnits().
# 2008-05-09
# o Now CnagCfhFile inherits from AromaMicroarrayDataFile.
# 2007-06-11
# o BUG FIX: readUnits() of CnagCfhFile was broken because it used the
#   non-existing variable 'nbrOfBytes'.
# 2007-04-05
# o Created.
############################################################################
