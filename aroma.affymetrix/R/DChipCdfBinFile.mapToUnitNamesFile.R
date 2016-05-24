# DChipCdfBinFile.mapToUnitNamesFile.R
setMethodS3("mapToUnitNamesFile", "DChipCdfBinFile", function(this, unf=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'unf':
  if (!is.null(unf)) {
    unf <- Arguments$getInstanceOf(unf, "UnitNamesFile");
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
 

  verbose && enter(verbose, "Generating index map of unit names from CDF.bin to UnitNamesFile (typically a CDF)");
  chipType <- getChipType(this);
  verbose && cat(verbose, "Chip type: ", chipType);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate Affymetrix CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(unf)) {
    verbose && enter(verbose, "Locating UnitNamesFile");
    cdf <- AffymetrixCdfFile$byChipType(chipType);
    unf <- cdf;
    verbose && exit(verbose);
  }

  verbose && print(verbose, "Target UnitNamesFile:");
  verbose && print(verbose, unf);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="mapToUnitNamesFile", class=class(this)[1], chipType=chipType, unf=list(platform=getPlatform(unf), chipType=getChipType(unf), nbrOfUnits=nbrOfUnits(unf), fullName=getFullName(unf), fileSize=getFileSize(unf)));
  dirs <- c("aroma.affymetrix", "dChip", chipType);
  if (!force) {
    units <- loadCache(key=key, dirs=dirs);
    if (!is.null(units)) {
      verbose && cat(verbose, "Found cached results.");
      verbose && exit(verbose);
      return(units);
    }
  }

  verbose && enter(verbose, "Reading unit names from CDF.bin");
  verbose && cat(verbose, "Number of units: ", nbrOfUnits(this));
  unitNames <- getUnitNames(this);
  verbose && str(verbose, unitNames);
  verbose && exit(verbose);

  verbose && enter(verbose, "Querying target UnitNamesFile");
  units <- indexOf(unf, names=unitNames);
  # Not needed anymore
  unitNames <- NULL;
  verbose && str(verbose, units);
  verbose && exit(verbose);

  saveCache(units, key=key, dirs=dirs);
 
  verbose && exit(verbose);

  units;
}, protected=TRUE)

############################################################################
# HISTORY:
# 2009-02-13
# o Created.
############################################################################
