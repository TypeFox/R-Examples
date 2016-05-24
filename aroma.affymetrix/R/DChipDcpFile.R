###########################################################################/**
# @RdocClass DChipDcpFile
#
# @title "The DChipDcpFile class"
#
# \description{
#  @classhierarchy
#
#  A DChipDcpFile object represents a DChip DCP file.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \seealso{
#   @see "DChipDcpSet".
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("DChipDcpFile", function(...) {
  this <- extend(AffymetrixFile(...), "DChipDcpFile",
    "cached:.header" = NULL
  );

#  setCdf(this, cdf);

  this;
})


setMethodS3("as.character", "DChipDcpFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("File format: %s", getFileFormat(this)));
  s <- c(s, sprintf("Number of cells: %s", nbrOfCells(this)));
  s <- c(s, sprintf("Number of units: %s", nbrOfUnits(this)));
  s <- c(s, sprintf("Has normalized data: %s", hasNormalizedData(this)));
  s <- c(s, sprintf("Has MBEI estimates: %s", hasMbeiData(this)));
  # Has CDF?
#  cdf <- getCdf(this);
#  if (!is.null(cdf)) {
#    s <- c(s, sprintf("Chip type: %s", getChipType(cdf)));
#  }
  s;
}, protected=TRUE)


setMethodS3("getFileFormat", "DChipDcpFile", function(this, ...) {
  hdr <- getHeader(this);

  ver <- sprintf("v%d", as.integer(hdr$Format));

  ver;
})


setMethodS3("getExtensionPattern", "DChipDcpFile", function(static, ...) {
  "[.](dcp|DCP)$";
}, static=TRUE, protected=TRUE)



setMethodS3("fromFile", "DChipDcpFile", function(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
  df <- newInstance(static, filename=filename, path=path, ...);

  # Try to read the header
  hdr <- getHeader(df);

  df;
}, protected=TRUE)



setMethodS3("getHeader", "DChipDcpFile", function(this, force=FALSE, ...) {
  hdr <- this$.header;
  if (force || is.null(hdr)) {
    hdr <- dChipIO::readDcpHeader(getPathname(this), ...);
    this$.header <- hdr;
  }
  hdr;
})

setMethodS3("nbrOfUnits", "DChipDcpFile", function(this, ...) {
  hdr <- getHeader(this);
  if (hdr$Format %in% c(3,4)) {
    nbrOfBytes <- file.info(getPathname(this))$size;
    fileHeaderSize <- 3028;
    nbrOfCells <- nbrOfCells(this);
    nbrOfUnitBytes <- nbrOfBytes - fileHeaderSize - 2*2*nbrOfCells;
    nbrOfUnits <- as.integer(nbrOfUnitBytes / 13);
  } else {
    throw("Cannot infer number of units for DCP file format v", hdr$Format);
  }

  nbrOfUnits;
})


setMethodS3("nbrOfCells", "DChipDcpFile", function(this, ...) {
  prod(dim(this, ...));
})

setMethodS3("dim", "DChipDcpFile", function(x) {
  # To please R CMD check
  this <- x;

  hdr <- getHeader(this);
  rep(hdr$CellDim, 2);
}, appendVarArgs=FALSE)


setMethodS3("hasNormalizedData", "DChipDcpFile", function(this, ...) {
  hdr <- getHeader(this, ...);
  as.logical(hdr$Normalized);
})


setMethodS3("hasMbeiData", "DChipDcpFile", function(this, ...) {
  hdr <- getHeader(this, ...);
  as.logical(hdr$ThetaValid);
})


setMethodS3("getRawIntensities", "DChipDcpFile", function(this, cells=NULL, force=FALSE, ...) {
  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, max=nbrOfCells(this));
  }

  nbrOfUnits <- nbrOfUnits(this);
  field <- "rawIntensities";
  data <- dChipIO::readDcp(getPathname(this), fields=field, cells=cells, ...);

  data <- data[[field]];

  data;
})


setMethodS3("getNormalizedIntensities", "DChipDcpFile", function(this, cells=NULL, force=FALSE, ...) {
  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, max=nbrOfCells(this));
  }

  nbrOfUnits <- nbrOfUnits(this);
  field <- "normalizedIntensities";
  data <- dChipIO::readDcp(getPathname(this), fields=field, cells=cells, ...);

  data <- data[[field]];

  data;
})


setMethodS3("getCalls", "DChipDcpFile", function(this, units=NULL, force=FALSE, ...) {
  nbrOfUnits <- nbrOfUnits(this);
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits);
  }

  field <- "calls";
  data <- dChipIO::readDcp(getPathname(this), fields=field, units=units, ...);

  data <- data[[field]];
  data <- as.integer(data) + as.integer(1);

  levels <- c("SNP_A", "SNP_B", "SNP_AB", "SNP_NO_CALL", "SNP_A_B", "SNP_B_A");
  levels <- c("AA", "BB", "AB", "NC", "ab", "ba");
  attr(data, "levels") <- levels;
  class(data) <- "factor";

  data;
})


setMethodS3("getThetas", "DChipDcpFile", function(this, units=NULL, force=FALSE, ...) {
  nbrOfUnits <- nbrOfUnits(this);
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits);
  }

  field <- "thetas";
  data <- dChipIO::readDcp(getPathname(this), fields=field, units=units, ...);

  data <- data[[field]];

  data;
})


setMethodS3("getThetaStds", "DChipDcpFile", function(this, units=NULL, force=FALSE, ...) {
  nbrOfUnits <- nbrOfUnits(this);
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits);
  }

  field <- "thetaStds";
  data <- dChipIO::readDcp(getPathname(this), fields=field, units=units, ...);

  data <- data[[field]];

  data;
})


setMethodS3("getExcludes", "DChipDcpFile", function(this, units=NULL, force=FALSE, ...) {
  nbrOfUnits <- nbrOfUnits(this);
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits);
  }

  field <- "excludes";
  data <- dChipIO::readDcp(getPathname(this), fields=field, units=units, ...);

  data <- data[[field]];

  data;
})


setMethodS3("getThetasAB", "DChipDcpFile", function(this, units=NULL, force=FALSE, ...) {
  nbrOfUnits <- nbrOfUnits(this);
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits);
  }

  fields <- c("thetas", "thetaStds");
  data <- dChipIO::readDcp(getPathname(this), fields=fields, units=units, ...);

  data <- data[fields];
  data <- unlist(data, use.names=FALSE);
  data <- matrix(data, ncol=2);
  colnames(data) <- c("thetaA", "thetaB");

  data;
})


setMethodS3("extractTheta", "DChipDcpFile", function(this, units=NULL, ..., drop=FALSE, nbrOfGroups=NULL, verbose=FALSE) {
  # Arguments 'nbrOfGroups':
  if (is.null(nbrOfGroups)) {
    nbrOfGroups <- this$.nbrOfGroups;
    if (is.null(nbrOfGroups))
      nbrOfGroups <- 1;
  }

  if (nbrOfGroups == 1) {
    data <- getThetas(this, units=units, ..., verbose=verbose);
    if (!drop) {
      data <- as.matrix(data);
    }
  } else if (nbrOfGroups == 2) {
    data <- getThetasAB(this, units=units, ..., verbose=verbose);
  }

  # Drop singleton dimensions?
  if (drop) {
    data <- drop(data);
  }

  data;
})



# setMethodS3("getCdf", "DChipDcpFile", function(this, ...) {
#   this$cdf;
# })
#
#
# setMethodS3("setCdf", "DChipDcpFile", function(this, cdf, ...) {
#   # Argument 'cdf':
#   if (!is.null(cdf)) {
#     cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");
#   }
#
#   this$cdf <- cdf;
# })
#
# setMethodS3("getUnitMap", "DChipDcpFile", function(this, ...) {
#   cdf <- getCdf(this);
#   if (is.null(cdf)) {
#     throw("Cannot infer the number of units. No CDF specified.");
#   }
#
#   unitNames <- getUnitNames(cdf, ...);
#   pattern <- "^(AFFX|Random)";
#   units <- which(regexpr(pattern, unitNames) == -1);
#   units;
# })
#
# setMethodS3("nbrOfUnits", "DChipDcpFile", function(this, ...) {
#   length(getUnitMap(this, ...));
# })



##############################################################################
# HISTORY:
# 2008-08-20
# o Updated DChipDcpFile to utilize the new dChipIO package.
# 2008-05-09
# o Now DChipDcpFile inherits from GenericDataFile.
# 2008-01-30
# o Created.
##############################################################################
