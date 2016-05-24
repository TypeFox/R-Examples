setMethodS3("allocateFromCdf", "AromaUnitGcContentFile", function(static, cdf, ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  types <- "double";
  sizes <- 4L;

  res <- NextMethod("allocateFromCdf", types=types, sizes=sizes);
  res[,1] <- as.double(NA);

  res;
}, static=TRUE)



setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUnitGcContentFile", function(this, csv, rows=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csv':
  csv <- Arguments$getInstanceOf(csv, "AffymetrixNetAffxCsvFile");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Importing (unit name, GC content) data from ", class(csv)[1]);

  # Query CDF
  cdf <- getCdf(this);

  # Read data
  data <- readDataFrame(csv, colClasses=c("^(probeSetID|%GC)$"="character")); 
  unitNames <- data[,1];
  verbose && str(verbose, unitNames);
  values <- as.double(data[,2]);
  verbose && str(verbose, values);
  # Not needed anymore
  data <- NULL;
  
  # Keep the units in the CDF
  units <- indexOf(cdf, names=unitNames);
  verbose && cat(verbose, "Unit indices:");
  verbose && str(verbose, units);
  # Not needed anymore
  unitNames <- NULL;

  keep <- which(!is.na(units));
  verbose && printf(verbose, "Keeping %d of %d (%.2f%%)\n", 
       length(keep), length(units), 100*length(keep)/length(units));
  units <- units[keep];
  # Not needed anymore
  keep <- NULL;
  if (length(units) == 0) {
    warning("None of the unit names in the CSV match the ones in the CDF ('", getPathname(cdf), "'). Is the correct file ('", getPathname(csv), "'), being imported?");
  }

  # Update
  this[units,1] <- values;
  # Not needed anymore
  values <- NULL;

  verbose && exit(verbose);

  invisible(units);
})


############################################################################
# HISTORY:
# 2009-03-22
# o Created from AromaUflFile.R.
############################################################################
