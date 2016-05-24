# @author "HB"
setMethodS3("justSNPRMA", "character", function(...) {
  requireNamespace("oligo") || throw("Package not loaded: oligo")
  oligo::justSNPRMA(...);
})


# @author "HB"
setMethodS3("justSNPRMA", "AffymetrixCelSet", function(this, ..., normalizeToHapmap=TRUE, normalizeSNPsOnly="auto", returnESet=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'normalizeToHapmap':
  normalizeToHapmap <- Arguments$getLogical(normalizeToHapmap);

  # Argument 'normalizeSNPsOnly':
  if (normalizeSNPsOnly == "auto") {
  } else {
    normalizeSNPsOnly <- Arguments$getLogical(normalizeSNPsOnly);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Running SNPRMA on ", class(this)[1]);

  csR <- this;
  cdf <- getCdf(csR);
  chipType <- getChipType(cdf, fullname=FALSE);
  hasCNs <- (regexpr("^GenomeWideSNP_(5|6)$", chipType) != -1);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize SNPs only?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (identical(normalizeSNPsOnly, "auto")) {
    normalizeSNPsOnly <- hasCNs;
  }

  # Get the SNP only tag
  if (normalizeSNPsOnly) {
    snpOnlyTag <- "SNPs";
  } else {
    snpOnlyTag <- NULL;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Rank-based quantile normalization
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Rank-based quantile normalization");

  verbose && enter(verbose, "Setting up normalization model");
  if (normalizeToHapmap) {
    refTag <- "HapMapRef";
  } else {
    refTag <- NULL;
  }
  qn <- QuantileNormalization(csR, targetDistribution=NULL,
                                subsetToAvg=NULL, typesToUpdate="pm",
                                        tags=c("*", snpOnlyTag, refTag));

  if (!isDone(qn)) {
    if (normalizeToHapmap) {
      verbose && enter(verbose, "Loading HapMap reference target quantiles");

      pdPkgName <- .cleanPlatformName(chipType);
      verbose && cat(verbose, "Platform Design (PD) package: ", pdPkgName);

      # Load target from PD package
      path <- system.file(package=pdPkgName);
      if (path == "") {
        throw("Cannot load HapMap reference target quantiles. Package not installed: ", pdPkgName);
      }

      path <- file.path(path, "extdata");
      path <- Arguments$getReadablePath(path);

      verbose && enter(verbose, "Loading binary file");
      filename <- sprintf("%sRef.rda", pdPkgName);
      pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
      verbose && cat(verbose, "Pathname: ", pathname);
      target <- loadToEnv(pathname)$reference;
      verbose && str(verbose, target);
      verbose && exit(verbose);

      qn$.targetDistribution <- target;
      # Not needed anymore
      target <- NULL;
      verbose && exit(verbose);
    } # if (normalizeToHapMap)


    if (normalizeSNPsOnly) {
      verbose && enter(verbose, "Identifying cells of SNPs for fitting normalization function");

      # justSNPRMA() operates only on SNP* units (e.g. CN units ignored).
      # For this reason we here *estimate* the normalization function based
      # on these units only, but for convenience we will apply it to all
      # units (including CN units, if they exist).
      verbose && enter(verbose, "Identifying units");
      pattern <- "^SNP";
      verbose && cat(verbose, "Pattern: ", pattern);
      units <- indexOf(cdf, pattern=pattern);
      verbose && cat(verbose, "Units:");
      verbose && str(verbose, units);
      verbose && exit(verbose);

      verbose && enter(verbose, "Identifying cell indices of these units");
      cells <- getCellIndices(cdf, units=units, unlist=TRUE, useNames=FALSE);
      verbose && cat(verbose, "Cells:");
      verbose && str(verbose, cells);
      # Not needed anymore
      units <- NULL;
      verbose && exit(verbose);

      qn$.subsetToAvg <- cells;
      # Not needed anymore
      cells <- NULL;
      verbose && exit(verbose);
    } # if (normalizeSNPsOnly)
  } # if (!isDone(qn))

  verbose && print(verbose, qn);
  verbose && exit(verbose);

  verbose && enter(verbose, "Processing");
  csN <- process(qn, verbose=verbose);
  verbose && print(verbose, csN);
  verbose && exit(verbose);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Probe-level summarization
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Fitting probe-level (summarization) model");
  # We use the oligo estimator for fitting the log-additive model
  plm <- RmaSnpPlm(csN, mergeStrands=FALSE, flavor="oligo");
  verbose && print(verbose, plm);

  if (length(findUnitsTodo(plm)) > 0) {
    if (hasCNs) {
      verbose && enter(verbose, "Fitting CN probes");
      units <- fitCnProbes(plm, verbose=verbose);
      verbose && cat(verbose, "CN units fitted:");
      verbose && str(verbose, units);
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Fitting remaining units");
    units <- fit(plm, verbose=verbose);
    verbose && cat(verbose, "Units fitted:");
    verbose && str(verbose, units);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extracting chip effect set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting ChipEffectSet");
  ces <- getChipEffectSet(plm);
  verbose && print(verbose, ces);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extracting eSet
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (returnESet) {
    verbose && enter(verbose, "Extracting eSet");
    pkg <- Package("oligo");
    if (!isOlderThan(Package("oligo"), "1.12.0")) {
      # For oligo v1.12.0 and newer
      eSet <- extractAlleleSet(ces, verbose=verbose);
    } else {
      # For oligo v1.11.x and older
      if (hasCNs) {
        eSet <- extractSnpCnvQSet(ces, verbose=verbose);
      } else {
        eSet <- extractSnpQSet(ces, verbose=verbose);
      }
    }
    verbose && print(verbose, eSet);
    verbose && exit(verbose);

    res <- eSet;
  } else {
    res <- ces;
  }

  # Return result
  res;
})


############################################################################
# HISTORY:
# 2011-11-20
# o BUG FIX: Internally, justSNPRMA() would pass argument 'verbose=log'
#   instead of 'verbose=verbose', which would throw an error unless
#   'log' was assigned to be a logical value or a Verbose object.
# 2010-05-09
# o Made justSNPRMA(..., normalizeSNPsOnly="auto") for AffymetrixCelSet
#   the default.
# 2010-05-06
# o Now justSNPRMA(..., returnESet=TRUE) for AffymetrixCelSet returns an
#   AlleleSet due to updates of classes in oligo v1.12.0.
# 2009-05-10
# o BUG FIX: In the most recent version of oligo, its justSNPRMA() requires
#   that oligo is loaded.  Updated justSNPRMA() for character
#   to assert this.
# 2009-01-10
# o Added argument 'normalizeSNPsOnly'.
# 2008-12-05
# o Created.
############################################################################
