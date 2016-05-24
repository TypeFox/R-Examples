setConstructorS3("AromaUgpFile", function(...) {
  this <- extend(AromaUnitChromosomeTabularBinaryFile(...), "AromaUgpFile");

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})

setMethodS3("getFilenameExtension", "AromaUgpFile", function(static, ...) {
  "ugp";
}, static=TRUE)

setMethodS3("getDefaultExtension", "AromaUgpFile", function(static, ...) {
  "ugp";
}, static=TRUE, protected=TRUE);

setMethodS3("getExtensionPattern", "AromaUgpFile", function(static, ...) {
  "[.](ugp)$";
}, static=TRUE, protected=TRUE)

setMethodS3("getDefaultColumnNames", "AromaUgpFile", function(this, ...) {
  c("chromosome", "position");
}, protected=TRUE)


setMethodS3("getPositions", "AromaUgpFile", function(this, units=NULL, ...) {
  # It is possible to pass NULL as an index here.
  cc <- indexOfColumn(this, name="position");
  this[units,cc,drop=TRUE];
}, protected=TRUE)



setMethodS3("getUnitsOnChromosome", "AromaUgpFile", function(this, ..., region=NULL) {
  # Arguments 'region':
  if (!is.null(region)) {
  }

  units <- NextMethod("getUnitsOnChromosome");

  # Subset by regions
  if (!is.null(region)) {
    cc <- indexOfColumn(this, "position");
    data <- this[units,cc,drop=TRUE];
    keep <- !is.na(data);
    keep <- keep & (region[1] <= data & data <= region[2]);
    units <- units[keep];
  }

  units;
}, protected=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: File I/O
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("allocate", "AromaUgpFile", function(static, ...) {
  NextMethod("allocate", types=rep("integer", times=2), sizes=c(1,4));
}, static=TRUE, protected=TRUE)



setMethodS3("importFromGenericTabularFile", "AromaUgpFile", function(this, src, colClasses=c("*"="NULL", "^Probe Set ID$"="character", "^Chromosome$"="character", "^Physical Position$"="character"), colOrder=NULL, shift=0, con=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'src':
  src <- Arguments$getInstanceOf(src, "GenericTabularFile");

  # Argument 'colOrder':
  if (!is.null(colOrder)) {
    colOrder <- Arguments$getIndices(colOrder, length=3);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Importing (unitName, chromosome, position) from ", class(src)[1], " file");

  data <- readDataFrame(src, colClasses=colClasses, ..., verbose=less(verbose));

  # Rearrange columns (optional)
  if (!is.null(colOrder))
    data <- data[,colOrder,drop=FALSE];

  # Map to unit names
  unf <- getUnitNamesFile(this);
  unfUnitNames <- getUnitNames(unf);
  unfUnits <- match(data[,1], unfUnitNames);

  # Exclude units that are not in the annotation unit names file
  keep <- which(!is.na(unfUnits));
  unfUnits <- unfUnits[keep];
  if (length(unfUnits) == 0) {
    warning("None of the imported unit names match the ones in the annotation unit names file ('", getPathname(unf), "'). Is the correct file ('", getPathname(src), "'), being imported?");
  }

  # Assume 'chromosome' is in 2nd column, and 'position' in 3rd.
  data <- data[keep,2:3,drop=FALSE];

  # Not needed anymore
  keep <- NULL;
  gc <- gc();

  # Convert chromosome strings to integers
  if (!is.integer(data[,1])) {
    map <- c(X=23, Y=24, Z=25);
    for (kk in seq_along(map)) {
      data[,1] <- gsub(names(map)[kk], map[kk], data[,1]);
    }
    suppressWarnings({
      data[,1] <- as.integer(data[,1]);
    })
    gc <- gc();
  }

  # Convert positions to integers
  if (!is.integer(data[,2])) {
    suppressWarnings({
      data[,2] <- as.integer(data[,2]);
    })
    gc <- gc();
  }

  # Shift positions?
  if (shift != 0) {
    verbose && printf(verbose, "Shifting positions %d steps.", shift);
    data[,2] <- data[,2] + as.integer(shift);
  }

  # Update to file
  this[unfUnits,1] <- data[,1];
  this[unfUnits,2] <- data[,2];

  verbose && exit(verbose);

  invisible(unfUnits);
}, protected=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: File I/O
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# OBSOLETE?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getUnitsAt", "AromaUgpFile", function(this, chromosomes, region=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument '...':
  args <- list(...);
  if ("chromosome" %in% names(args)) {
    chromosomes <- args[["chromosome"]];
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Cached in memory?
  key <- list(method="getUnitsAt", class="AromaUgpFile",
                  chromosomes=chromosomes, region=region);
  key <- getChecksum(key);
  if (!force) {
    res <- this$.memoryCache[["getUnitsAt"]][[key]];
    if (!is.null(res))
      return(res);
  }

  # Stratify by chromosome
  data <- this[,1,drop=TRUE];

  # Update known chromosomes, if not already done.
  allChromosomes <- getChromosomes(this, .chromosomes=data);

  keep <- !is.na(data) & (data %in% chromosomes);
  idxs <- which(keep);

  if (!is.null(region)) {
    data <- this[idxs,2,drop=TRUE];
    keep <- !is.na(data);
    keep <- keep & (region[1] <= data & data <= region[2]);
    idxs <- idxs[keep];
  }

  if (!is.list(this$.memoryCache))
    this$.memoryCache <- list();

  cache <- list(idxs); names(cache) <- key;
  this$.memoryCache[["getUnitsAt"]] <- cache;

  idxs;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2009-05-08
# o Can we drop getUnitsAt()?
# o Now AromaUgpFile identifies the position column by its name.
# o Now AromaUgpFile inherits from AromaUnitChromosomeTabularBinaryFile.
# o Now getUnitsOnChromosomes() of AromaUgpFile returns a list which one
#   element per requested chromosomes.
# 2009-05-03
# o Now readDataFrame() of AromaUgpFile accepts rows=integer(0).
# 2009-02-12
# o Added getUnitsOnChromosomes() to AromaUgpFile.
# 2008-05-24
# o Now allocate() of AromaUgpFile adds footer 'createdOn'.
# 2008-05-21
# o Added getUnitsOnChromosome() and getPositions().
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
