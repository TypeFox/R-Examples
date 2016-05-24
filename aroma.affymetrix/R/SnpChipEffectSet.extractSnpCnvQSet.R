setMethodS3("extractSnpCnvQSet", "SnpChipEffectSet", function(this, units=NULL, sortUnits=TRUE, transform=log2, ..., verbose=FALSE) {
  requireNamespace("Biobase") || throw("Package not loaded: Biobase")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this, fullname=FALSE);

  if (regexpr("^GenomeWideSNP_(5|6)$", chipType) == -1) {
    throw("Cannot extract SnpCnvQSet: Unsupported chip type: ", chipType);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'units':
  if (is.null(units)) {
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extracting SnpCnvQSet from ", class(this)[1]);

  if (is.null(units)) {
    verbose && enter(verbose, "Identifying all SNP_A-* units");
    # Identify all SNP_A-* units (as what is returned by oligo)
    units <- indexOf(cdf, pattern="^SNP_A-");
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Order units lexicographically by their names?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  unitNames <- getUnitNames(cdf, units=units);
  if (sortUnits) {
    verbose && enter(verbose, "Sorting units by their names");
    srt <- sort(unitNames, method="quick", index.return=TRUE);
    unitNames <- srt$x;
    units <- units[srt$ix];
    # Not needed anymore
    srt <- NULL;  # Not needed anymore
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting data");
  theta <- extractTheta(this, groups=1:2, units=units, verbose=verbose);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate and populate SnpCnvQSet
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocate and populate SnpCnvQSet");
  res <- new("SnpCnvQSet",
    thetaA = transform(theta[,1,,drop=TRUE]),
    thetaB = transform(theta[,2,,drop=TRUE])
  );

  # Not needed anymore
  # Not needed anymore
  theta <- NULL;

  # Assign feature data
  .featureNames(res) <- unitNames;
  # Not needed anymore
  unitNames <- NULL;

  # Assign annotation data
  pdPkgName <- .cleanPlatformName(chipType);
  .annotation(res) <- pdPkgName;

  # Assign sample names
  filenames <- sapply(this, getFilename);
  names(filenames) <- NULL;
  filenames <- gsub(",chipEffects", "", filenames);
  .sampleNames(res) <- filenames;

  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # extractSnpCnvQSet()


############################################################################
# HISTORY:
# 2012-09-01
# o ROBUSTNESS: extractSnpCnvQSet() for SnpChipEffectSet would throw an
#   exception if the 'Biobase' package was not loaded.
# 2008-12-05
# o Created.
############################################################################
