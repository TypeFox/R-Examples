setMethodS3("extractSnpQSet", "SnpChipEffectSet", function(this, units=NULL, sortUnits=TRUE, transform=log2, ..., verbose=FALSE) {
  requireNamespace("Biobase") || throw("Package not loaded: Biobase")


  # Assert oligo version
  pkg <- Package("oligo");
  if (!isOlderThan(pkg, "1.12.0")) {
    throw("extractSnpQSet() requires oligo v1.12.0 or older. Instead use extractAlleleSet(): ", getVersion(pkg));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this, fullname=FALSE);

  if (regexpr("^GenomeWideSNP_(5|6)$", chipType) != -1) {
    throw("Cannot extract SnpQSet: Unsupported chip type: ", chipType);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'units':
  if (is.null(units)) {
    # Identify all SNP_A-* units (as what is returned by oligo)
    units <- indexOf(cdf, pattern="^SNP_A-");
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Order units lexicographically by their names?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  unitNames <- getUnitNames(cdf, units=units);
  if (sortUnits) {
    verbose && enter(verbose, "Sort units by their names");
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
  theta <- extractTheta(this, groups=1:4, units=units, verbose=verbose);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Swap (antisense, sense) unit groups to (sense, antisense)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Ordering unit groups to be (sense, antisense)");
  # Make sure pairs are order as (sense, antisense)
  dirs <- getGroupDirections(cdf, units=units);
  names(dirs) <- NULL;
  # Not needed anymore
  units <- NULL;

  # Sanity check
  lens <- sapply(dirs, FUN=length);
  uLens <- unique(lens);
  if (any(!is.element(uLens, c(2,4)))) {
    throw("Internal error: Unexpected number of unit groups: ",
                                              paste(uLens, collapse=", "));
  }

  # Extract the direction/strand of the first group
#  dirs <- lapply(dirs, FUN=function(groups) groups[1]);
  dirs <- lapply(dirs, FUN=.subset, 1L);
  dirs <- unlist(dirs, use.names=FALSE);

  # Identify which to swap from (antisense,sense) to (sense,antisense)
  idxs <- which(dirs == 2);
  # Not needed anymore
  dirs <- NULL;

  verbose && cat(verbose, "Swapping elements:");
  verbose && str(verbose, idxs);
  theta[idxs,,] <- theta[idxs,c(3,4,1,2),];
  # Not needed anymore
  idxs <- NULL;   # Not needed anymore
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate and populate SnpQSet
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocate and populate SnpQSet");
  res <- new("SnpQSet",
    senseThetaA     = transform(theta[,1,,drop=TRUE]),
    senseThetaB     = transform(theta[,2,,drop=TRUE]),
    antisenseThetaA = transform(theta[,3,,drop=TRUE]),
    antisenseThetaB = transform(theta[,4,,drop=TRUE])
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

  res;
})


############################################################################
# HISTORY:
# 2012-09-01
# o ROBUSTNESS: extractSnpQSet() for SnpChipEffectSet would throw an
#   exception if the 'Biobase' package was not loaded.
# 2010-05-06
# o extractSnpQSet() now asserts that oligo v1.12.0 or older is installed.
# 2008-12-05
# o Created.
############################################################################
