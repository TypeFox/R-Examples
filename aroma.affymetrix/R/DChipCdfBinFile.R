###########################################################################/**
# @RdocClass DChipCdfBinFile
#
# @title "The DChipCdfBinFile class"
#
# \description{
#  @classhierarchy
#
#  A DChipCdfBinFile object represents a DChip CDF.bin file.
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
# @author "HB"
#*/###########################################################################
setConstructorS3("DChipCdfBinFile", function(...) {
  this <- extend(AffymetrixFile(...), c("DChipCdfBinFile",
                                                       uses("UnitNamesFile")),
    "cached:.header" = NULL,
    "cached:.unitNames" = NULL
  );

  this;
})


setMethodS3("as.character", "DChipCdfBinFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  s <- c(s, sprintf("File format: %s", getFileFormat(this)));
  s <- c(s, sprintf("Number of cells: %s", nbrOfCells(this)));
  s <- c(s, sprintf("Number of units: %s", nbrOfUnits(this)));
  s;
}, protected=TRUE)


setMethodS3("getFileFormat", "DChipCdfBinFile", function(this, ...) {
  hdr <- getHeader(this);
  ver <- sprintf("v%d", as.integer(hdr$Format));
  ver;
})


setMethodS3("getChipType", "DChipCdfBinFile", function(this, fullname=TRUE, ...) {
  hdr <- getHeader(this);
  chipType <- hdr$ChipType;
  if (!fullname) {
    chipType <- gsub(",.*", "", chipType);
  }
  chipType;
})

setMethodS3("getPlatform", "DChipCdfBinFile", function(this, ...) {
  "affymetrix";
})


setMethodS3("nbrOfUnits", "DChipCdfBinFile", function(this, ...) {
  hdr <- getHeader(this);
  nbrOfUnits <- hdr$NumUnit;
  nbrOfUnits;
})

setMethodS3("getCellDim", "DChipCdfBinFile", function(this, ...) {
  hdr <- getHeader(this);
  hdr$CellDim;
}, protected=TRUE);

setMethodS3("nbrOfCells", "DChipCdfBinFile", function(this, ...) {
  nbrOfCells <- getCellDim(this)^2;
  nbrOfCells <- as.integer(nbrOfCells);
  nbrOfCells;
})


setMethodS3("findByChipType", "DChipCdfBinFile", function(this, chipType, tags=NULL, ...) {
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  fullname <- paste(c(chipType, tags), collapse=",");
  parts <- strsplit(fullname, split=",")[[1]];
  chipType <- parts[1];
  tags <- parts[-1];
  pattern <- sprintf("^%s.cdf.bin$", fullname);
  pathname <- findAnnotationDataByChipType(chipType, pattern=pattern);

  pathname;
}, protected=TRUE)


setMethodS3("byChipType", "DChipCdfBinFile", function(this, ...) {
  pathname <- findByChipType(this, ...);
  DChipCdfBinFile(pathname);
})


setMethodS3("fromFile", "DChipCdfBinFile", function(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
  df <- newInstance(static, filename=filename, path=path, ...);
  # Try to read the header
  hdr <- getHeader(df);
  df;
}, protected=TRUE)


setMethodS3("getHeader", "DChipCdfBinFile", function(this, force=FALSE, ...) {
  hdr <- this$.header;
  if (force || is.null(hdr)) {
    hdr <- dChipIO::readCdfBinHeader(getPathname(this), ...);
    this$.header <- hdr;
  }
  hdr;
})


setMethodS3("getUnitNames", "DChipCdfBinFile", function(this, units=NULL, ...) {
  # Arguments 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(this));
  }

  names <- this$.unitNames;

  if (is.null(names)) {
    names <- readDataFrame(this, ..., fields="unitName")[[1]];
    this$.unitNames <- names;
  }

  if (!is.null(units))
    names <- names[units];

  names;

})


setMethodS3("getUnitSizes", "DChipCdfBinFile", function(this, ...) {
  readDataFrame(this, ..., fields="unitSize")[[1]];
})

setMethodS3("readDataFrame", "DChipCdfBinFile", function(this, units=NULL, fields=c("unitName", "unitSize", "cellPos"), ...) {
  # Arguments 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(this));
  }

  data <- dChipIO::readCdfBin(getPathname(this), units=units, ...);

  # Keep only fields of interest
  names <- names(data);
  names <- gsub("unitNames", "unitName", names, fixed=TRUE);
  names <- gsub("numProbes", "unitSize", names, fixed=TRUE);
  names <- gsub("CellPos", "cellPos", names, fixed=TRUE);
  keep <- which(is.element(names, fields));
  data <- data[keep];
  names(data) <- fields;

  # Coerce into a data frame
  data <- as.data.frame(data, stringsAsFactors=FALSE);

  data;
})

##############################################################################
# HISTORY:
# 2009-02-13
# o Added argument fullname=TRUE to getChipType() of DChipCdfBinFile.
# 2008-08-20
# o Created.
##############################################################################
