# @author "HB"
setConstructorS3("AffymetrixCnChpFile", function(..., cdf=NULL) {
  this <- extend(AffymetrixFile(...), "AffymetrixCnChpFile",
    "cached:.header" = NULL,
    "cached:.unitReadMap" = NULL,
    .cdf = NULL
  );

  if (!is.null(cdf))
    setCdf(this, cdf);

  pathname <- getPathname(this)
  if (!is.null(pathname)) {
    # Make sure the name is non-empty
    name <- getName(this);
    if (!nzchar(name)) {
      throw("An ", class(this)[1], " must have a name of at least length one: ", pathname);
    }
  }

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})

setMethodS3("as.character", "AffymetrixCnChpFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("File format: %s", getFileFormat(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(this))));
  s <- c(s, sprintf("Timestamp: %s", as.character(getTimestamp(this))));
  s <- c(s, sprintf("Unit read map: %s", capture.output(str(getUnitReadMap(this)))));
  s;
}, protected=TRUE)


setMethodS3("clone", "AffymetrixCnChpFile", function(this, ..., verbose=TRUE) {
  # Clone itself (and clear the cached fields)
  object <- NextMethod("clone", clear=TRUE);

  # Clone the CDF here.
  if (!is.null(object$.cdf))
    object$.cdf <- clone(object$.cdf);

  object;
}, protected=TRUE)


setMethodS3("getExtensionPattern", "AffymetrixCnChpFile", function(static, ...) {
  "[.](cnchp|CNCHP)$";
}, static=TRUE, protected=TRUE)



setMethodS3("getFileFormat", "AffymetrixCnChpFile", function(this, ...) {
  pathname <- getPathname(this);

  # Read CEL header
  raw <- readBin(pathname, what=raw(), n=10);

  if (raw[1] == 59)
    return("v1 (binary; CC)");

  if (raw[1] == 64)
    return("v4 (binary; XDA)");

  if (rawToChar(raw[1:5]) == "[CEL]")
    return("v3 (text; ASCII)");

  naValue <- as.character(NA);
  return(naValue);
})




setMethodS3("fromFile", "AffymetrixCnChpFile", function(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
  if (.checkArgs) {
    # Argument 'filename' and 'path':
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
  } else {
    pathname <- filename;
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  # WORKAROUND: Currently the affxparser code crash R if the file is not
  # a valid CNCHP file.  The best we can do now is to test against the
  # filename.
  isCnChp <- (regexpr("[.](cnchp|CNCHP)$", pathname) != -1);
  if (!isCnChp) {
    throw("Could not read CNCHP file. Filename format error: ", pathname);
  }

  # Create a new instance of the same class
  newInstance(static, pathname, ...);
}, static=TRUE, protected=TRUE)


setMethodS3("getChipType", "AffymetrixCnChpFile", function(this, ...) {
  getChipType(getCdf(this));
}, private=TRUE)


setMethodS3("getUnitNamesFile", "AffymetrixCnChpFile", function(this, ...) {
  getCdf(this, ...);
})

setMethodS3("getUnitTypesFile", "AffymetrixCnChpFile", function(this, ...) {
  getCdf(this, ...);
})

setMethodS3("getCdf", "AffymetrixCnChpFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    # Look up the chip type in the header
    params <- getHeader(this)$dataHeader$parameters;
    chipType <- params[["affymetrix-algorithm-param-chip_type"]];

    # Retrieve the corresponding CDF
    cdf <- AffymetrixCdfFile$byChipType(chipType);

    this$.cdf <- cdf;
  }
  cdf;
})


setMethodS3("setCdf", "AffymetrixCnChpFile", function(this, cdf, ..., .checkArgs=TRUE) {
  if (.checkArgs) {
    # Argument 'cdf':
    cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

#    # Assure that the CDF is compatible with the CNCHP file
#    if (nbrOfCells(cdf) != nbrOfCells(this)) {
#      throw("Cannot set CDF. The specified CDF structure is not compatible with the CEL file. The number of cells do not match: ", nbrOfCells(cdf), " != ", nbrOfCells(this));
#    }
  }

  # Have to clear the cache
  clearCache(this);

  this$.cdf <- cdf;

  invisible(this);
})


setMethodS3("getHeader", "AffymetrixCnChpFile", function(this, ...) {
  header <- this$.header;
  if (is.null(header)) {
    pathname <- getPathname(this);
    suppressWarnings({
      header <- .readCcgHeader(pathname);
    });
    this$.header <- header;
  }
  header;
}, private=TRUE)


setMethodS3("getTimestamp", "AffymetrixCnChpFile", function(this, format="%a %b %d %H:%M:%S * %Y", ...) {
  # Argument 'format':
  format <- Arguments$getCharacter(format, length=c(1,1));

  # Look up the chip type in the header
  params <- getHeader(this)$dataHeader$parameters;
  timestamp <- params[["affymetrix-algorithm-param-create_date"]];
  timestamp <- trim(timestamp); # Unnecessary?

  # Ignore parts of the timestamp?
  if (regexpr("[*]", format) != -1) {
    timestamp <- unlist(strsplit(timestamp, split=" ", fixed=TRUE));
    format <- unlist(strsplit(format, split=" ", fixed=TRUE));
    keep <- (format != "*");
    timestamp <- timestamp[keep];
    format <- format[keep];
    timestamp <- paste(timestamp, collapse=" ");
    format <- paste(format, collapse=" ");
  }

  # Parse the identified timestamp into POSIXct
  res <- strptime(timestamp, format=format, ...);

  # If no valid timestamp was found, return NA.
  if (length(as.character(res)) == 0) {
    res <- as.POSIXct(NA);
  }

  # Keep the non-parsed timetstamp etc for debugging.
  attr(res, "text") <- timestamp;
  attr(res, "params") <- params;

  res;
}, protected=TRUE)


setMethodS3("hasUnitReadMap", "AffymetrixCnChpFile", function(this, ...) {
  !is.null(this$.unitReadMap);
}, private=TRUE);


setMethodS3("setUnitReadMap", "AffymetrixCnChpFile", function(this, readMap=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Setting unit read map");
  verbose && cat(verbose, "Argument 'readMap':");
  verbose && str(verbose, readMap);

  oldReadMap <- this$.unitReadMap;

  if (is.null(readMap)) {
  } else if (is.character(readMap)) {
    unitNamesOnFile <- readMap;
    unf <- getUnitNamesFile(this);
    unitNames <- getUnitNames(unf);
    verbose && cat(verbose, "Unit names (CDF):");
    verbose && str(verbose, unitNames);
    readMap <- match(unitNames, unitNamesOnFile);
  }

  verbose && cat(verbose, "readMap:");
  verbose && str(verbose, readMap);

  this$.unitReadMap <- readMap;

  verbose && exit(verbose);

  invisible(oldReadMap);
}, private=TRUE);


setMethodS3("getUnitReadMap", "AffymetrixCnChpFile", function(this, ...) {
  this$.unitReadMap;
}, private=TRUE)


setMethodS3("readRawData", "AffymetrixCnChpFile", function(this, fields=c("ProbeSetName", "Log2Ratio"), ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fields':
  fields <- intersect(fields, c("ProbeSetName", "Log2Ratio"));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting raw data from CNCHP file");
  pathname <- getPathname(this);
  verbose && cat(verbose, "Pathname: ", pathname);
  verbose && cat(verbose, "Fields: ", paste(fields, collapse=", "));

  verbose && enter(verbose, "Reading data from file");
#  res <- .readCcg(pathname, subset=list(CopyNumber=list(rows=..., fields=...));
  res <- .readCcg(pathname);
  verbose && exit(verbose);

  res <- res$dataGroups$MultiData$dataSets$CopyNumber$table;
  row.names(res) <- NULL;
  verbose && cat(verbose, "Raw data read (all fields):");
  verbose && str(verbose, res);

  # Trim all strings (should really be done by readCcg(). /HB 2008-08-23)
  for (kk in seq_len(ncol(res))) {
    values <- res[,kk];
    if (is.character(values)) {
      values <- trim(values);
      res[,kk] <- values;
    }
    # Not needed anymore
    values <- NULL;
  }

  if (!hasUnitReadMap(this)) {
    if ("ProbeSetName" %in% names(res)) {
      setUnitReadMap(this, res[,"ProbeSetName"], verbose=less(verbose, 1));
    }
  }

  verbose && enter(verbose, "Extracting fields of interest");
  if (!is.null(fields))
    res <- res[,fields, drop=FALSE];
  verbose && cat(verbose, "Raw data read:");
  verbose && str(verbose, res);
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}, protected=TRUE) # readRawData()



setMethodS3("getData", "AffymetrixCnChpFile", function(this, units=NULL, fields=c("ProbeSetName", "Log2Ratio"), ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  nbrOfUnits <- nbrOfUnits(getCdf(this));
  if (is.null(units)) {
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits);
    nbrOfUnits <- length(units);
  }

  # Argument 'fields':
  fields <- intersect(fields, c("ProbeSetName", "Log2Ratio"));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting data from CNCHP file");
  pathname <- getPathname(this);
  verbose && cat(verbose, "Pathname: ", pathname);
  verbose && cat(verbose, "Fields: ", paste(fields, collapse=", "));
  verbose && cat(verbose, "Units: ", capture.output(str(units)));

  # Check for cached results
  key <- list(method="getData", class="AffymetrixCnChpFile",
              header=unlist(getHeader(this)), fields=fields, units=units);
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="getData", header=unlist(getHeader(this)), fields=fields, units=units);
  }
  dirs <- c("aroma.affymetrix", "GTC");
  res <- loadCache(key=key, dirs=dirs);
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "Cached results found.");
    verbose && exit(verbose);
    return(res);
  }

  # Assert that we can retrieve the CDF
  cdf <- getCdf(this);

  # Get the raw data
  verbose && enter(verbose, "Reading raw data (takes time)");
  res <- readRawData(this, fields=fields, ..., verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Remapping rows (units) according to the CDF");
  # Extract units of interest
  readMap <- getUnitReadMap(this);  # Have been by readRawData()
  str(readMap);
  # Convert CDF unit indices to CNCHP file unit inidices
  if (is.null(units)) {
    units <- readMap;
  } else {
    units <- readMap[units];
  }
  verbose && cat(verbose, "Unit indices (on file):");
  verbose && str(verbose, units);
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting units of interest");
  res <- res[units,, drop=FALSE];
  row.names(res) <- NULL;
  verbose && str(verbose, res);
  verbose && exit(verbose);

  verbose && enter(verbose, "Caching to file");
  saveCache(res, key=key, dirs=dirs);
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # getData()


setMethodS3("extractLogRatios", "AffymetrixCnChpFile", function(this, units=NULL, ..., drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  data <- getData(this, units=units, ..., verbose=verbose);
  data <- data[,"Log2Ratio"];
  data;
})


############################################################################
# HISTORY:
# 2009-07-08
# o Added getUnitTypesFile() for AffymetrixCnChpFile.
# 2008-08-22
# o Added extractLogRatios().
# 2008-05-18
# o Added support for UnitNamesFile.
# 2008-01-16
# o Created.
############################################################################

