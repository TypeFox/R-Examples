# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: AFFX
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("allocateFromCdf", "AromaUflFile", function(static, cdf, nbrOfEnzymes=1L, ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  # Argument 'nbrOfEnzymes':
  nbrOfEnzymes <- Arguments$getInteger(nbrOfEnzymes, range=c(1,10));

  types <- rep("integer", times=nbrOfEnzymes);
  sizes <- rep(2L, times=nbrOfEnzymes);

  NextMethod("allocateFromCdf", types=types, sizes=sizes);
}, static=TRUE)



setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUflFile", function(this, csv, enzymes=1:nbrOfEnzymes(this), enzymesToUpdate=1:nbrOfEnzymes(this), rows=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csv':
  csv <- Arguments$getInstanceOf(csv, "AffymetrixNetAffxCsvFile");

  # Argument 'enzymesToUpdate':
  enzymesToUpdate <- Arguments$getIndices(enzymesToUpdate, max=nbrOfEnzymes(this));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Importing (unit name, fragment length+) data from ", class(csv)[1]);

  # Query CDF
  cdf <- getCdf(this);

  # Read unit names
  unitNames <- getUnitNames(csv, rows=rows, ..., verbose=less(verbose));
  verbose && cat(verbose, "Unit names:");
  verbose && str(verbose, unitNames);
  
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
  if (length(units) == 0) {
    warning("None of the unit names in the CSV match the ones in the CDF ('", getPathname(cdf), "'). Is the correct file ('", getPathname(csv), "'), being imported?");
  }

#  verbose && summary(verbose, data);

  # Read data
  data <- readDataUnitFragmentLength(csv, enzymes=enzymes, rows=keep, ..., verbose=less(verbose));
  verbose && str(verbose, data);
#  verbose && summary(verbose, data);

  data <- data[,-1,drop=FALSE];
  verbose && str(verbose, data);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc, level=-10);

  # Update
  this[units,enzymesToUpdate] <- data;
  # Not needed anymore
  data <- NULL;

  gc <- gc();
  verbose && print(verbose, gc, level=-10);

  verbose && exit(verbose);

  invisible(units);
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: AFFX
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


############################################################################
# HISTORY:
# 2008-09-15
# o Added argument 'enzymesToUpdate' to importFromAffymetrixNetAffxCsvFile()
#   in order to make it possible to specify both which enzymes to read 
#   and to update.
# 2008-04-24
# o Updated importFromAffymetrixNetAffxCsvFile() to read the unit names,
#   then identify the ones that are in the CDF, and the only read the data
#   for the corresponding rows.  This was done so that one can see from the
#   verbose output that by excluding the extra units the remaining ones are
#   sound.
# 2008-04-14
# o Renamed readData() to readDataFrame() for AromaTabularBinaryFile.
# 2007-12-08
# o BUG FIX: importFromAffymetrixNetAffxCsvFile() of AromaUflFile failed
#   to import one of two enzymes.
# 2007-09-14
# o Added support for multiple fragment lengths, in case multiple enzymes
#   were used for the same assay.
# 2007-09-11
# o Created.
############################################################################
