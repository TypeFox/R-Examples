# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: Platform specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("importFromGenomeInformation", "AromaUgpFile", function(this, gi, ..., verbose=FALSE) {
  # Argument 'gi':
  gi <- Arguments$getInstanceOf(gi, "GenomeInformation");

  # AD HOC patch, since units==NULL does not work./HB 2007-03-03
  units <- seq_len(nbrOfUnits(gi));
  data <- getData(gi, units=units, fields=c("chromosome", "physicalPosition"));

  chr <- data[,"chromosome"];
  if (is.character(chr)) {
    chr[chr == "X"] <- 23;
    chr[chr == "Y"] <- 24;
    chr[chr %in% c("MT", "Z")] <- 25;
    suppressWarnings({
      chr <- as.integer(chr);
    })
  }
  
  pos <- data[,"physicalPosition"];
  suppressWarnings({
    pos <- as.integer(pos);
  })

  this[,1] <- chr;
  this[,2] <- pos;

  # A best guess of what was imported
  units <- units[!(is.na(chr) & is.na(pos))];

  invisible(units);
})



setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUgpFile", function(this, csv, shift="auto", onReplicates=c("median", "mean", "overwrite"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csv':
  csv <- Arguments$getInstanceOf(csv, "AffymetrixNetAffxCsvFile");

  # Argument 'shift':
  if (!identical(shift, "auto"))
    shift <- Arguments$getInteger(shift);

  # Argument 'onReplicates':
  onReplicates <- match.arg(onReplicates);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Importing (unit name, chromosome, position) data from ", class(csv)[1]);

  # Query CDF
  unf <- getUnitNamesFile(this);
  unfUnitNames <- getUnitNames(unf);

  # Read data
  data <- readDataUnitChromosomePosition(csv, ..., verbose=less(verbose));
  importNames <- attr(data, "importNames");

  # Map to unit names
  unfUnits <- match(data[,1], unfUnitNames);

  # Exclude units that are not in the annotation unit names file
  keep <- which(!is.na(unfUnits));
  unfUnits <- unfUnits[keep];
  if (length(unfUnits) == 0) {
    warning("None of the imported unit names match the ones in the annotation unit names file ('", getPathname(unf), "'). Is the correct file ('", getPathname(csv), "'), being imported?");
  }
  data <- data[keep,2:3,drop=FALSE];
  importNames <- importNames[2:3];

  # Garbage collect
  # Not needed anymore
  keep <- NULL;
  gc <- gc();


  # Shift positions?
  if (identical(shift, "auto")) {
    shift <- 0;
    if (any(regexpr("[sS]tart", importNames) != -1))
      shift <- 13;
  }
  if (shift != 0) {
    verbose && printf(verbose, "Shifting positions %d steps.\n", shift);
    data[,2] <- data[,2] + as.integer(shift);
  }
 
  # Replicated positions per unit?
  dups <- which(duplicated(unfUnits));
  if (length(dups) > 0) {
    verbose && enter(verbose, "Detected units with replicated positions");
    dupUnits <- unique(unfUnits[dups]);
    nDupUnits <- length(dupUnits);
    verbose && cat(verbose, "Number of units with replicated positions: ", 
                                                               nDupUnits);

    if (onReplicates %in% c("median", "mean")) {
      verbose && enter(verbose, "Calculate average positions for those (assuming they are on the same chromosome)");
      if (onReplicates == "mean") {
        avgFcn <- median;
      } else {
        avgFcn <- mean;
      }
      for (kk in seq_along(dupUnits)) {
        if (kk %% 500 == 0)
         verbose && printf(verbose, "%d, ", kk);
        dupUnit <- dupUnits[kk];
        # Identify position
        units <- which(unfUnits == dupUnit);
        # Average position
        avgPos <- median(data[units,2], na.rm=TRUE);
        avgPos <- round(avgPos);
        # Update (can we update just units[1]?)
        data[units,2] <- avgPos;
      }
      verbose && cat(verbose, kk);
      verbose && exit(verbose);
      verbose && enter(verbose, "Remove the extraneous cases");
      data <- data[-dups,,drop=FALSE];
      unfUnits <- unfUnits[-dups];
      verbose && exit(verbose);
      # Not needed anymore
      dupUnits <- units <- NULL;
      verbose && str(verbose, unfUnits);
      warning("The positions for ", nDupUnits, " units were calculated as the average of replicated positions, since that was what was available on file.");
    } else {
      verbose && cat(verbose, "Ignored replicated positions. The last one written will be the one available.");
    }
    verbose && exit(verbose);
  }
  # Not needed anymore
  dups <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc, level=-10);

  this[unfUnits,1] <- data[,1];
  this[unfUnits,2] <- data[,2];

  # Not needed anymore
  data <- NULL;
  gc <- gc();
  verbose && print(verbose, gc, level=-10);

  verbose && exit(verbose);

  invisible(unfUnits);
})



setMethodS3("importFromAffymetrixTabularFile", "AromaUgpFile", function(this, src, colClasses=c("*"="NULL", "^probeSetID$"="character", "^chromosome$"="character", "^(physicalPosition|position)$"="character"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'src':
  src <- Arguments$getInstanceOf(src, "AffymetrixTabularFile");

  units <- importFromGenericTabularFile(this, src=src, 
            colClasses=colClasses, camelCaseNames=TRUE, ...);

  invisible(units);
}, protected=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: Platform specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


############################################################################
# HISTORY:
# 2009-05-12
# o Removed getUnitNamesFile() from AromaUgpFile. This requires 
#   aroma.core v1.0.6 or newer.
# 2008-07-21
# o Removed allocateFromCdf() specific to AromaUgpFile since that is
#   now taken care of in the super classes.
# 2008-05-18
# o Made class a little bit less platform specific by utilizing new
#   UnitNamesFile interface.
# 2008-05-12
# o Added static allocate().
# 2008-04-29
# o BUG FIX: Name clash in getUnitsAt() after new argument 'chromosomes'.
# 2008-04-17
# o Renamed argument 'chromosome' of getUnitsAt() of AromaUgpFile to 
#   'chromosomes'.  This was done in order to make it consistent with 
#   getUnitsOnChromosome() of GenomeInformation. Thanks Tim Keighley at 
#   CSIRO for pointing this out.
# 2008-04-14
# o Renamed readData() to readDataFrame() for AromaTabularBinaryFile.
# 2007-09-16
# o Now importFromAffymetrixNetAffxCsvFile() averages positions if multiple
#   positions were available for a particular unit.
# o Now importFromGenomeInformation() tries to return units imported and
#   not all.  This is still a best guess, but still more informative than
#   before.
# 2007-09-14
# o Added getChromosomes(), which caches results in memory.
# o Added importFromAffymetrixTabularFile() to AromaUgpFile.
# 2007-09-13
# o Removed createFromGenomeInformation().
# o Updated AromaUgpFile according to changes in super class.
# 2007-09-10
# o Added importFromAffymetrixNetAffxCsvFile() to AromaUgpFile.
# 2007-03-04
# o Added findByChipType() and fromChipType().
# o Now the default path for createFromCdf() is the same as for the CDF.
# 2007-03-03
# o Now inherits from generic AromaGenomePositionFile.
# 2007-03-02
# o Created. Can import genome information data.
############################################################################
