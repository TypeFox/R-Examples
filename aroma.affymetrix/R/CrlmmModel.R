setConstructorS3("CrlmmModel", function(dataSet=NULL, balance=1.5, minLLRforCalls=c(5, 1, 5), recalibrate=TRUE, flavor="v2", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  isMappingChipType <- FALSE;
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "SnpChipEffectSet");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Sanity checks
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    chipType <- getChipType(dataSet, fullname=FALSE);

    # For now, only allow know SNP chip types. /HB 2008-12-07
    if (regexpr("^Mapping(10|50|250)K_.*$", chipType) != -1) {
      isMappingChipType <- TRUE;
    } else if (regexpr("^GenomeWideSNP_(5|6)$", chipType) != -1) {
      throw("Cannot fit CRLMM model: Model fitting for this chip type is not supported/implemented: ", chipType);
    } else {
      throw("Cannot fit CRLMM model: Unsupported/unsafe chip type: ", chipType);
    }

    # CRLMM does not apply to chip effects for which strand-specific
    # estimates have been combined (summed together).
    if (getMergeStrands(dataSet)) {
      throw("Cannot fit CRLMM model: CRLMM requires that the probe-level model was fitted without merging the strands (mergeStrands=FALSE).");
    }
  }

  # Argument 'balance':
  balance <- Arguments$getDouble(balance, range=c(0.00001, 1e6));

  # Argument 'minLLRforCalls':
  minLLRforCalls <- Arguments$getDoubles(minLLRforCalls, length=3,
                                                         range=c(0, 1e6));

  # Argument 'recalibrate':
  recalibrate <- Arguments$getLogical(recalibrate);

  # Argument 'flavor':
  flavor <- match.arg(flavor);


  extend(Model(dataSet=dataSet, ...), "CrlmmModel",
    balance = balance,
    minLLRforCalls = minLLRforCalls,
    recalibrate = recalibrate,
    flavor = flavor,
    .isMappingChipType = isMappingChipType
  )
})


setMethodS3("getRootPath", "CrlmmModel", function(this, ...) {
  "crlmmData";
}, protected=TRUE)


setMethodS3("getAsteriskTags", "CrlmmModel", function(this, collapse=NULL, ...) {
  tags <- "CRLMM";
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
}, protected=TRUE)


setMethodS3("getParameters", "CrlmmModel", function(this, ...) {
  params <- NextMethod("getParameters");
  params$balance <- this$balance;
  params$minLLRforCalls <- this$minLLRforCalls;
  params$recalibrate <- this$recalibrate;
  params$flavor <- this$flavor;
  params;
}, protected=TRUE)




setMethodS3("getChipType", "CrlmmModel", function(this, ...) {
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  chipType <- getChipType(cdf, ...);
  chipType <- gsub(",monocell", "", chipType);
  chipType;
})



setMethodS3("getCallSet", "CrlmmModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  ces <- getDataSet(this);
  cdf <- getCdf(ces);
  chipType <- getChipType(this, fullname=FALSE);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocating parameter files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up output directory
  outPath <- getPath(this);

  # Setting up output pathnames
  fullnames <- getFullNames(ces);
  fullnames <- gsub(",chipEffects$", "", fullnames);
  filenames <- sprintf("%s,genotypes.acf", fullnames);

  nbrOfArrays <- length(ces);
  nbrOfUnits <- nbrOfUnits(cdf);
  platform <- getPlatform(cdf);

  verbose && enter(verbose, "Retrieving genotype call set");
  agcList <- list();
  for (kk in seq_along(filenames)) {
    filename <- filenames[kk];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, filename, nbrOfArrays));
    pathname <- filePath(outPath, filename);
    if (isFile(pathname)) {
      agc <- AromaUnitGenotypeCallFile(pathname);
    } else {
      verbose && enter(verbose, "Allocating new file");
     chipTypeF <- getChipType(this);
      agc <- AromaUnitGenotypeCallFile$allocate(filename=pathname, platform=platform, chipType=chipTypeF, nbrOfRows=nbrOfUnits, verbose=verbose);
      verbose && exit(verbose);
    }
    agcList[[kk]] <- agc;
    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  res <- AromaUnitGenotypeCallSet$byPath(outPath);

  res;
}) # getCallSet()



setMethodS3("getConfidenceScoreSet", "CrlmmModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  ces <- getDataSet(this);
  cdf <- getCdf(ces);
  chipType <- getChipType(this, fullname=FALSE);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocating parameter files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up output directory
  outPath <- getPath(this);

  # Setting up output pathnames
  fullnames <- getFullNames(ces);
  fullnames <- gsub(",chipEffects$", "", fullnames);
  filenames <- sprintf("%s,confidenceScores.acf", fullnames);

  nbrOfArrays <- length(ces);
  nbrOfUnits <- nbrOfUnits(cdf);
  platform <- getPlatform(cdf);

  verbose && enter(verbose, "Retrieving genotype call confidence set");
  agcList <- list();
  for (kk in seq_along(filenames)) {
    filename <- filenames[kk];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, filename, nbrOfArrays));
    pathname <- filePath(outPath, filename);
    if (isFile(pathname)) {
      agc <- AromaUnitSignalBinaryFile(pathname);
    } else {
      verbose && enter(verbose, "Allocating new file");
     chipTypeF <- getChipType(this);
      agc <- AromaUnitSignalBinaryFile$allocate(filename=pathname, platform=platform, chipType=chipTypeF, nbrOfRows=nbrOfUnits, verbose=verbose);
      naValue <- as.double(NA);
      agc[,1] <- naValue;
      verbose && exit(verbose);
    }
    agcList[[kk]] <- agc;
    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  res <- AromaUnitSignalBinarySet$byPath(outPath, pattern=",confidenceScores.acf$");

  res;
}) # getConfidenceScoreSet()


setMethodS3("getCrlmmParametersSet", "CrlmmModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  ces <- getDataSet(this);
  cdf <- getCdf(ces);
  chipType <- getChipType(this, fullname=FALSE);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocating parameter files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up output directory
  outPath <- getPath(this);

  # Setting up output pathnames
  fullnames <- getFullNames(ces);
  fullnames <- gsub(",chipEffects$", "", fullnames);
  filenames <- sprintf("%s,CRLMM.atb", fullnames);

  nbrOfArrays <- length(ces);
  nbrOfUnits <- nbrOfUnits(cdf);
  platform <- getPlatform(cdf);

  verbose && enter(verbose, "Retrieving CRLMM parameters set");
  atbList <- list();
  for (kk in seq_along(filenames)) {
    filename <- filenames[kk];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, filename, nbrOfArrays));
    pathname <- filePath(outPath, filename);
    if (isFile(pathname)) {
      atb <- CrlmmParametersFile(pathname);
    } else {
      verbose && enter(verbose, "Allocating new file");
      chipTypeF <- getChipType(this);
      atb <- CrlmmParametersFile$allocate(filename=pathname, platform=platform, chipType=chipTypeF, nbrOfRows=nbrOfUnits, verbose=verbose);
      verbose && exit(verbose);
    }
    atbList[[kk]] <- atb;
    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  res <- CrlmmParametersSet$byPath(outPath);

  res;
}) # getCrlmmParametersSet()



setMethodS3("getUnitsToFit", "CrlmmModel", function(this, ...) {
  # Identify which units CRLMM have fitted
  getCrlmmSNPs(this, ...);
})


setMethodS3("findUnitsTodo", "CrlmmModel", function(this, units=NULL, safe=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  unitsToFit <- getUnitsToFit(this, ...);
  verbose && str(verbose, unitsToFit);

  if (safe) {
    unitsTodo <- unitsToFit;
  } else {
    # Find all units with missing-value (NA) calls
    callSet <- getCallSet(this);
    unitsWithNAs <- findUnitsTodo(callSet, ...);
    unitsTodo <- intersect(unitsToFit, unitsWithNAs);
    verbose && str(verbose, unitsTodo);
  }

  if (!is.null(units)) {
    unitsTodo <- intersect(units, unitsTodo);
    verbose && str(verbose, unitsTodo);
  }

  unitsTodo;
})


setMethodS3("fit", "CrlmmModel", function(this, units="remaining", force=FALSE, ram=NULL, ..., verbose=FALSE) {
  # Early error
  requireNamespace("oligoClasses") || throw("Package not loaded: oligoClasses");
  .getM <- oligoClasses::getM

  requireNamespace("oligo") || throw("Package not loaded: oligo");

  # To please R CMD check
  ns <- loadNamespace("oligo");
  fitAffySnpMixture <- get("fitAffySnpMixture", mode="function", envir=ns);
  getInitialAffySnpCalls <- get("getInitialAffySnpCalls", mode="function", envir=ns);
  getAffySnpGenotypeRegionParams <- get("getAffySnpGenotypeRegionParams", mode="function", envir=ns);
  getGenotypeRegionParams <- get("getGenotypeRegionParams", mode="function", envir=ns);
  updateAffySnpParams <- get("updateAffySnpParams", mode="function", envir=ns);
  replaceAffySnpParams <- get("replaceAffySnpParams", mode="function", envir=ns);
  getAffySnpDistance <- get("getAffySnpDistance", mode="function", envir=ns);
  getAffySnpCalls <- get("getAffySnpCalls", mode="function", envir=ns);
  getAffySnpConfidence <- get("getAffySnpConfidence", mode="function", envir=ns);
  getAffySnpGenotypeRegionParams <- get("getAffySnpGenotypeRegionParams", mode="function", envir=ns);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pkg <- Package("oligo");
  if (!isOlderThan(Package("oligo"), "1.12.0")) {
    # For oligo v1.12.0 and newer

    extractESet <- function(..., hasQuartets=TRUE) {
      eSet <- extractAlleleSet(...);
      eSet;
    } # extractESet()

    extractLogRatios <- function(eSet, ...) {
      ad <- .assayData(eSet);
      M <- ad$senseAlleleA - ad$senseAlleleB;
      # Not needed anymore
      ad <- NULL;
      M;
    } # extractLogRatios()
  } else {
    # For oligo v1.11.x and older

    # To please R CMD check
    ns <- loadNamespace("oligo");
    thetaA <- get("thetaA", mode="function", envir=ns);
    thetaB <- get("thetaB", mode="function", envir=ns);

    extractESet <- function(..., hasQuartets=TRUE) {
      if (hasQuartets) {
        eSet <- extractSnpQSet(...);
      } else {
        eSet <- extractSnpCnvQSet(...);
      }
      eSet;
    } # extractESet()

    extractLogRatios <- function(eSet, ...) {
      M <- thetaA(eSet) - thetaB(eSet);
      M;
    } # extractLogRatios()
  }


  maleIndex <- c();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'units':
  doRemaining <- FALSE;
  if (is.null(units)) {
  } else if (is.numeric(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  } else if (identical(units, "remaining")) {
    doRemaining <- TRUE;
  } else {
    throw("Unknown mode of argument 'units': ", mode(units));
  }

  # Argument 'ram':
  if (identical(ram, "oligo")) {
  } else {
    # Argument 'ram':
    ram <- getRam(aromaSettings, ram);
  }

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling genotypes by CRLMM");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get algorithm parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- getParameters(this);
  balance <- params$balance;
  minLLRforCalls <- params$minLLRforCalls;
  recalibrate <- params$recalibrate;
  isMappingChipType <- this$.isMappingChipType;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identifying units to process
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying units to process");
  if (is.null(units)) {
    unitsToFit <- getUnitsToFit(this, verbose=less(verbose,1));
    units <- unitsToFit;
    # Not needed anymore
    unitsToFit <- NULL;
  } else if (doRemaining) {
    verbose && enter(verbose, "Identifying non-estimated units")
    units <- findUnitsTodo(this, safe=FALSE, verbose=less(verbose));
    verbose && str(verbose, units);
    verbose && exit(verbose);
  } else {
    unitsToFit <- getUnitsToFit(this, verbose=less(verbose,1));
    units <- intersect(units, unitsToFit);
    # Not needed anymore
    unitsToFit <- NULL;
    # Fit only unique units
    units <- unique(units);
  }
  nbrOfUnits <- length(units);
  verbose && str(verbose, units);
  verbose && printf(verbose, "Getting model fit for %d units.\n", nbrOfUnits);

  # Identify which of the requested units have *not* already been estimated
  if (!doRemaining) {
    if (force) {
      verbose && printf(verbose, "All of these are forced to be fitted.\n");
    } else {
      units <- findUnitsTodo(this, units=units, safe=FALSE, verbose=less(verbose));
      nbrOfUnits <- length(units);
      verbose && printf(verbose, "Out of these, %d units need to be fitted.\n", nbrOfUnits);
    }
  }
  verbose && exit(verbose);

  # Nothing to do?
  if (nbrOfUnits == 0) {
    verbose && exit(verbose);
    return(invisible(NULL));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setup");
  crlmm <- getCrlmmPriors(this, verbose=less(verbose,1));

  ces <- getDataSet(this);
  nbrOfArrays <- length(ces);
  data <- data.frame(gender=rep("female", times=nbrOfArrays));
  phenoData <- new("AnnotatedDataFrame", data=data);


  allSNPs <- getCrlmmSNPs(this, ..., verbose=less(verbose, 1));
  verbose && str(verbose, allSNPs);
  snpsOnChrX <- getCrlmmSNPsOnChrX(this, ..., verbose=less(verbose,1));
  verbose && str(verbose, snpsOnChrX);

  # Sanity check
  if (length(allSNPs) != length(crlmm$hapmapCallIndex)) {
    throw("Internal error: The number of identified SNPs and the number of prior HapMap calls does not match: ", length(allSNPs), " != ", length(crlmm$hapmapCallIndex));
  }

  dim <- dim(crlmm$params$centers);
  hasQuartets <- (length(dim) == 3);
  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get result set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  callSet <- getCallSet(this, verbose=less(verbose,1));
  paramSet <- getCrlmmParametersSet(this, verbose=less(verbose,1));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Processing units in chunk
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (identical(ram, "oligo")) {
    chunkSize <- 40000;
  } else {
    chunkSize <- ram * (500e3/nbrOfArrays);
  }
  unitList <- splitInChunks(units, chunkSize=chunkSize);
  nbrOfChunks <- length(unitList);

  # To keep a SNR estimate per array
  snrPerArray <- double(nbrOfArrays);

  chunk <- 1;
  while (length(unitList) > 0) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", chunk, nbrOfChunks));
    unitsChunk <- unitList[[1]];
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, unitsChunk);
    unitList <- unitList[-1];

    verbose && enter(verbose, "Extracting data");
    eSet <- extractESet(ces, units=unitsChunk, sortUnits=FALSE,
                             hasQuartets=hasQuartets, verbose=verbose);
    .phenoData(eSet) <- phenoData;
    verbose && exit(verbose);

    unitNames <- .featureNames(eSet);
    verbose && cat(verbose, "Unit names:");
    verbose && str(verbose, unitNames);

    verbose && enter(verbose, "Extract CRLMM priors");
    idxs <- match(unitNames, names(allSNPs));
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, idxs);
    hapmapCallIndex <- crlmm$hapmapCallIndex[idxs];
    params <- crlmm$params;
    if (hasQuartets) {
      params$centers <- params$centers[idxs,,,drop=FALSE];
      params$scales <- params$scales[idxs,,,drop=FALSE];
    } else {
      params$centers <- params$centers[idxs,,drop=FALSE];
      params$scales <- params$scales[idxs,,drop=FALSE];
    }
    params$N <- params$N[idxs,,drop=FALSE];
    # Not needed anymore
    idxs <- NULL;
    verbose && cat(verbose, "CRLMM prior parameters (estimated from HapMap):");
    verbose && str(verbose, params);
    verbose && exit(verbose);

    verbose && enter(verbose, "Fitting SNP mixtures");
    correction <- fitAffySnpMixture(eSet, verbose=as.logical(verbose));
    verbose && str(verbose, correction);
    verbose && cat(verbose, "SNR per array:");
    verbose && str(verbose, as.vector(correction$snr));
    verbose && exit(verbose); # "Fitting SNP mixtures"

    verbose && enter(verbose, "Genotype calling");
    naValue <- as.integer(NA);
    calls <- matrix(naValue, nrow=nrow(eSet), ncol=ncol(eSet));
    index <- which(!hapmapCallIndex);
    if (length(index) > 0) {
      verbose && enter(verbose, "Initial SNP calling");
      verbose && str(verbose, index);
      if (isMappingChipType) {
        # NOTE: Do not specify sqsClass=class(eSet). Instead it should
        # default to sqsClass="SnpQSet" regardless of the class of 'eSet'.
        # See oligo:::justCRLMMv3() of oligo v1.12.0. /HB 2010-05-06
        verbose && str(verbose, correction);
        calls[index,] <- getInitialAffySnpCalls(correction, subset=index, verbose=as.logical(verbose));
      } else {
        throw("Not implemented for GWS chip types");
      }
      verbose && exit(verbose);
    }
    verbose && exit(verbose);

    verbose && enter(verbose, "Estimate genotype regions");
    if (isMappingChipType) {
      # NOTE: Do not specify sqsClass=class(eSet). Instead it should
      # default to sqsClass="SnpQSet" regardless of the class of 'eSet'.
      # See oligo:::justCRLMMv3() of oligo v1.12.0. /HB 2010-05-06
      rparams <- getAffySnpGenotypeRegionParams(eSet, initialcalls=calls, f=correction$fs, subset=index, verbose=as.logical(verbose));
    } else {
      M <- extractLogRatios(eSet);
      rparams <- getGenotypeRegionParams(M, calls, correction$fs, verbose=as.logical(verbose));
    }
    verbose && exit(verbose);

    priors <- crlmm$priors;
    if (hasQuartets) {
      verbose && enter(verbose, "Updating sense & antisense genotype regions");
      eSet1 <- eSet[,1];
      M <- .getM(eSet1);  # From 'oligoClasses' (formely in 'oligo')
      dimnames(M) <- NULL;
      M <- M[,1,,drop=TRUE];
      oneStrand <- integer(nrow(M));
      for (ss in 1:2) {
        oneStrand[is.na(M[,ss])] <- ss;
      }
      rparams <- updateAffySnpParams(rparams, priors, oneStrand, verbose=as.logical(verbose));
      params  <- replaceAffySnpParams(params, rparams, index);
      dist <- getAffySnpDistance(eSet, params, correction$fs);
      dist[,,-2,] <- balance*dist[,,-2,];
      verbose && exit(verbose);
    } else {
      verbose && enter(verbose, "Updating genotype regions");
      rparams <- updateAffySnpParams(rparams, priors, verbose=as.logical(verbose));
      params  <- replaceAffySnpParams(params, rparams, index);
      dist <- getAffySnpDistance(eSet, params, correction$fs);
      dist[,,-2] <- balance*dist[,,-2];
      verbose && exit(verbose);
    }
    # Not needed anymore
    params <- index <- NULL;

    indexX <- which(is.element(unitNames, names(snpsOnChrX)));
    # NOTE: Do not specify sqsClass=class(eSet). Instead it should
    # default to sqsClass="SnpQSet" regardless of the class of 'eSet'.
    # See oligo:::justCRLMMv3() of oligo v1.12.0. /HB 2010-05-06
    calls <- getAffySnpCalls(dist, indexX, maleIndex, verbose=as.logical(verbose));
    llr <- getAffySnpConfidence(dist, calls, indexX, maleIndex, verbose=as.logical(verbose));

    if (recalibrate) {
      verbose && enter(verbose, "Recalibrating");
      naValue <- as.integer(NA);
      for (gg in 1:3) {
        calls[calls == gg & llr < minLLRforCalls[gg]] <- naValue;
      }
      # Not needed anymore
      llr <- NULL;
      # 'correction$snr' is a I vector, where I=#arrays.
      # Note that is has a dim() attribute making it look like a matrix.
      # Why exactly 3.675?!? /HB 2009-01-12
      calls[,(correction$snr < 3.675)] <- naValue;

      rparams <- getAffySnpGenotypeRegionParams(eSet, calls, correction$fs, verbose=as.logical(verbose));
      # Not needed anymore
      calls <- NULL;

      rparams <- updateAffySnpParams(rparams, priors, oneStrand);
      dist <- getAffySnpDistance(eSet, rparams, correction$fs, verbose=as.logical(verbose));
      if (hasQuartets) {
        dist[,,-2,] <- balance*dist[,,-2,];
      } else {
        dist[,,-2] <- balance*dist[,,-2];
      }
      # Not needed anymore
      oneStrand <- NULL;
      calls <- getAffySnpCalls(dist, indexX, maleIndex, verbose=as.logical(verbose));
      llr <- getAffySnpConfidence(dist, calls, indexX, maleIndex, verbose=as.logical(verbose));
      verbose && exit(verbose);
    } # if (recalibrate)

    snrPerArray <- snrPerArray + log(as.vector(correction$snr));
    verbose && cat(verbose, "Average SNR per array:");
    verbose && str(verbose, exp(snrPerArray/nbrOfChunks));

    # Clean up
    # Not needed anymore
    eSet <- indexX <- correction <- rparams <- priors <- NULL;

    verbose && enter(verbose, "Estimated genotype parameters");
    verbose && cat(verbose, "calls:");
    verbose && str(verbose, calls);
    verbose && cat(verbose, "llr:");
    verbose && str(verbose, llr);
    verbose && cat(verbose, "dist:");
    verbose && str(verbose, dist);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Storing results
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing CRLMM parameter estimates, confidence scores and genotypes");
    for (kk in seq_along(callSet)) {
      agc <- callSet[[kk]];
      atb <- paramSet[[kk]];
      verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, getName(agc), nbrOfArrays));

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (i) CRLMM specific parameter estimates
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "CRLMM specific parameter estimates");
      # LLR
      atb[unitsChunk,1] <- llr[,kk,drop=TRUE];

      # Distances
      distKK <- dist[,kk,,,drop=TRUE];
      dim(distKK) <- c(dim(distKK)[1], prod(dim(distKK)[2:3]));
      for (cc in 1:ncol(distKK)) {
        atb[unitsChunk,cc+1] <- distKK[,cc];
      }
      # Not needed anymore
      distKK <- atb <- NULL;
      verbose && exit(verbose); # "CRLMM specific parameter estimates"

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (ii) Genotype calls
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Genotype calls");

      verbose && enter(verbose, "Tranlating oligo calls to {NC,AA,AB,BB}");
      callsKK <- calls[,kk,drop=TRUE];
      callsT <- character(nrow(calls));
      callsT[(callsKK == 0)] <- "NC";
      callsT[(callsKK == 1)] <- "AA";
      callsT[(callsKK == 2)] <- "AB";
      callsT[(callsKK == 3)] <- "BB";
      # Not needed anymore
      callsKK <- NULL;
      verbose && exit(verbose); # "Tranlating oligo calls to {NC,AA,AB,BB}"

      updateGenotypes(agc, units=unitsChunk, calls=callsT,
                                              verbose=less(verbose,5));
      # Not needed anymore
      callsT <- NULL;

      verbose && cat(verbose, "Genotypes stored (as on file):");
      verbose && str(verbose, extractGenotypes(agc, units=unitsChunk));
      # Not needed anymore
      agc <- NULL;

      verbose && exit(verbose); # "Genotype calls"

      verbose && exit(verbose); # "Array #%d ('%s') of %d"
    } # for (kk ...)
    verbose && exit(verbose);  # "Storing CRLMM parameter estimates, confidence scores and genotypes"


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Next chunk
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    chunk <- chunk + 1;
    # Not needed anymore
    unitsChunk <- calls <- llr <- dist <- NULL;
    verbose && exit(verbose);  # "Chunk #%d of %d"
  } # while(length(unitsList) > 0)
  # Not needed anymore
  callSet <- NULL;

  verbose && enter(verbose, "Storing average SNR per arrays");

  # Calculate average SNR per array (on the non-log scale)
  snrPerArray <- exp(snrPerArray / nbrOfChunks);
  verbose && cat(verbose, "Average SNR per array (over all chunks):");
  verbose && str(verbose, snrPerArray);

  for (kk in seq_along(paramSet)) {
    pf <- paramSet[[kk]];
    updateParameter(pf, "snr", snrPerArray[kk], verbose=less(verbose, -20));
  }

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Confidence scores
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating confidence scores for each call");
  callSet <- calculateConfidenceScores(this, verbose=verbose);
  verbose && exit(verbose); # "Calculating confidence scores for each call"

  verbose && exit(verbose);

  # Return fitted units
  invisible(units);
}) # fit()


setMethodS3("calculateConfidenceScores", "CrlmmModel", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calculating confidence scores for each call");

  callSet <- getCallSet(this, verbose=less(verbose,1));
  confSet <- getConfidenceScoreSet(this, verbose=less(verbose,1));
  paramSet <- getCrlmmParametersSet(this, verbose=less(verbose,1));

  nbrOfArrays <- length(callSet);
  verbose && cat(verbose, "Number of arrays: ", nbrOfArrays);
  cf <- getOneFile(callSet);
  nbrOfUnits <- nbrOfUnits(cf);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);

  verbose && enter(verbose, "Retrieving SNR estimates");
  snrs <- sapply(paramSet, FUN=readParameter, name="snr", mode="double");
  snrs <- unlist(snrs, use.names=FALSE);
  # Sanity check
  if (length(snrs) != nbrOfArrays) {
    throw("Number of SNRs read from CRLMM parameter set does not match the number of arrays modelled: ", length(snrs), " != ", nbrOfArrays);
  }
  names(snrs) <- NULL;
  verbose && cat(verbose, "SNRs:");
  verbose && str(verbose, snrs);
  verbose && exit(verbose); # "Retrieving SNR estimates"

  # Ported from oligo:::LLR2conf()

  verbose && enter(verbose, "Retrieving prior spline parameters");
  splineParams <- getCrlmmSplineParameters(this, verbose=less(verbose,20));
  verbose && print(verbose, ll(envir=splineParams));
  verbose && exit(verbose); # "Retrieving prior spline parameters"

  verbose && enter(verbose, "Thresholding SNRs according to prior estimates");
  truncSnrs <- pmin(log(snrs), splineParams$SNRK);
  verbose && str(verbose, truncSnrs);
  verbose && exit(verbose); # "Thresholding SNRs according to prior estimates"

  verbose && enter(verbose, "Calculating transformed SNRs using prior lm fit");
  coef <- splineParams$SNRlm$coef;
  verbose && cat(verbose, "Prior parameters used:");
  verbose && str(verbose, coef);
  snrs2 <- coef[1] + coef[2]*truncSnrs;
  verbose && cat(verbose, "Transformed SNRs:");
  verbose && str(verbose, snrs2);
  # Not needed anymore
  coef <- NULL;
  verbose && exit(verbose); # "Calculating transformed SNRs using prior lm fit"

  # Allocate results
  naValue <- as.double(NA);
  conf <- matrix(naValue, nrow=nbrOfUnits, ncol=nbrOfArrays);

  params <- list(
    homozygotes=list(
      coefs = splineParams$lm1$coef,
      k2 = splineParams$HmzK2,
      k3 = splineParams$HmzK3
    ),
    heterozygotes=list(
      coefs = splineParams$lm2$coef,
      k2 = splineParams$HtzK2,
      k3 = splineParams$HtzK3
    )
  );
  verbose && str(verbose, "Prior parameters:");
  verbose && str(verbose, params);

  # If you call the genotype by chance, the probability that
  # you are correct is 1/3. Since CRLMM always call the genotype
  # and picks the one with greatest confidence, the smallest
  # possible confidence score is 1/3. /HB 2009-01-12
  minConf <- 1/3;

  for (kk in seq_len(nbrOfArrays)) {
    cf <- callSet[[kk]];
    pf <- paramSet[[kk]];
    sf <- confSet[[kk]];

    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, getName(cf), nbrOfArrays));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identifying heterozygote units
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Identifying heterozygote units");
    isHet <- isHeterozygous(cf, drop=TRUE, verbose=less(verbose, 25));
    verbose && exit(verbose); # "Identifying heterozygote units"

    # Identified units called by CRLMM.  This will for instance also
    # skip CN units.
    unitsKK <- which(!is.na(isHet));
    isHet <- isHet[unitsKK];
    nbrOfCalled <- length(unitsKK);
    verbose && printf(verbose, "Number of called units: %d (%.2f%%)\n",
                                      nbrOfCalled, 100*nbrOfCalled/nbrOfUnits);
    verbose && str(verbose, unitsKK);

    # Sanity check
    if (nbrOfCalled == 0) {
      throw("Cannot calculate confidence scores. Genotypes are not called: ",
                                                              getPathname(cf));
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Allocating/retrieving confidence scores
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    conf <- sf[unitsKK,1,drop=TRUE];

    # Identifying subset to be calculated
    if (force) {
    } else {
      verbose && enter(verbose, "Identifying subset of units not yet calculated units");
      idxs <- which(is.na(conf));
      unitsKK <- unitsKK[idxs];
      conf <- conf[idxs];
      # Not needed anymore
      idxs <- NULL;
      verbose && str(verbose, unitsKK);
      verbose && exit(verbose); # "Identifying subset of units not yet calculated units"
    }

    if (length(unitsKK) == 0) {
      verbose && cat(verbose, "Nothing to do.");
      verbose && exit(verbose);
      next;
    }

    # Special case
    if (snrs[kk] <= 3) {
      conf[unitsKK] <- minConf;
    } else {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Retrieving LLRs
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Retrieving LLRs");
      llr <- sqrt(pf[unitsKK,1,drop=TRUE]);
      verbose && str(verbose, llr);
      verbose && exit(verbose); # "Retrieving LLRs"


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Calculating confidence scores stratified by homozygotes & heterozygotes
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Calculating confidence scores stratified by homozygotes and heterozygotes");
      for (what in names(params)) {
        verbose && enter(verbose, sprintf("Stratify by '%s'", what));

        paramsT <- params[[what]];
        verbose && cat(verbose, "Parameters:");
        verbose && str(verbose, paramsT);

        if (what == "homozygotes") {
          idxs <- which(!isHet);
        } else {
          idxs <- which(isHet);
        }

        verbose && cat(verbose, "Number of ", what, ": ", length(idxs));
        if (length(idxs) > 0) {
          coefs <- paramsT$coefs;
          k2 <- paramsT$k2;
          k3 <- paramsT$k3;

          llrT <- pmin(llr[idxs], k3);
          confT <- coefs[1] + coefs[2]*llrT;

          delta <- llrT - k2;
          idxsT <- which(delta > 0);
          if (length(idxsT) > 0) {
            confT[idxsT] <- confT[idxsT] + coefs[3]*delta[idxsT];
          }

          conf[idxs] <- confT;
          # Not needed anymore
          idxsT <- delta <- llrT <- confT <- NULL;
        }
        # Not needed anymore
        idxs <- NULL;

        verbose && exit(verbose); # "Stratify by '%s'"
      } # for (what ...)
      verbose && exit(verbose); # "Calculating confidence scores stratified by homozygotes and heterozygotes"

      # Clean up
      # Not needed anymore
      isHet <- llr <- NULL;

      verbose && cat(verbose, "Confidence \"scores\":");
      verbose && str(verbose, conf);
      verbose && summary(verbose, conf);

      verbose && cat(verbose, "Corrected confidence \"scores\":");
      conf <- conf + snrs2[kk];
      verbose && summary(verbose, conf);

      verbose && cat(verbose, "Final confidence scores:");
      conf <- 1/(1+exp(-conf));

      conf[(conf < minConf)] <- minConf;
    } # if (snrs[kk] <= 3)

    verbose && cat(verbose, "Final confidence scores:");
    verbose && str(verbose, conf);
    verbose && summary(verbose, conf);

    # Sanity check
    if (any((conf < 0 | conf > 1), na.rm=TRUE)) {
      throw(verbose, "Internal error: Identified ", sum((conf < 0 | conf > 1), na.rm=TRUE), " confidence scores that were out of range");
    }

    verbose && enter(verbose, "Storing confidence scores");
    sf[unitsKK,1] <- conf;

    # Updating footer
    footer <- readFooter(sf);
    key <- "sourceFiles";
    data <- footer[[key]];
    if (is.null(data)) {
      data <- list();
    }
    srcFiles <- list(callFile=cf, paramFile=pf);
    for (name in names(srcFiles)) {
      srcFile <- srcFiles[[name]];
      attr <- list(
        nbrOfArrays = length(callSet),
        filename = getFilename(srcFile),
        filesize = getFileSize(srcFile),
        checksum = getChecksum(srcFile)
      );
      data[[name]] <- attr;
    }
    footer[[key]] <- data;
    writeFooter(sf, footer);
    # Not needed anymore
    footer <- srcFiles <- srcFile <- attr <- data <- key <- NULL;
    verbose && exit(verbose); # "Storing confidence scores"

    # Not needed anymore
    conf <- unitsKK <- NULL;

    # Next array
    verbose && exit(verbose); "Array #%d ('%s') of %d"
  } # for (kk ...)

  verbose && exit(verbose); # "Calculating confidence scores for each call"

  invisible(confSet);
}, protected=TRUE);



############################################################################
# HISTORY:
# 2013-06-02
# o BUG FIX: calculateConfidenceScores() for CrlmmModel used defunct
#   method isHeterozygote() instead of isHeterozygous().
# 2012-09-05
# o ROBUSTNESS: Now the CrlmmModel constructor asserts that the chip
#   effects were estimated without merging the strands (mergeStrands=FALSE).
#   If not, an informative exception is thrown.
# 2011-11-05
# o FIX: The verbose enter/exit statements of fit() of CrlmmModel()
#   did not match up, resulting deeper and deeper indentations.
# 2011-04-25
# o Clarified the error message that CRLMM is not supported for GWS6.
# 2010-05-06
# o Now CrlmmModel(..., recalibrate=TRUE) is the default.
# o BUG FIX: fit() of CrlmmModel would not work with oligo v1.12.0
#   and newer.
# o BUG FIX: getCallSet() and getCrlmmParametersSet() of CrlmmModel used
#   non-existing verbose object 'log' instead of 'verbose'.
# 2010-01-06
# o CLEAN UP: No need for assign NAs when allocating new files; this is now
#   always the default way (in aroma.core v1.4.1).
# 2009-05-17
# o BUG FIX: fit() for CrlmmModel was calling oligo::getM(), but that
#   method was later moved to oligoClasses.  Now we just do getM().
# 2009-01-12
# o Added calculateConfidenceScores().
# 2009-01-10
# o Updated to work with latest aroma.core and aroma.affymetrix.
# 2009-01-07
# o Overriding getAsteriskTags() and getRootPath().
# 2008-12-08
# o Now setup is much more like fit() for ProbeLevelModel.
# 2008-12-07
# o Starting to make use of AromaCrlmmBinarySet.
# o Created CrlmmModel from justCRLMMv2().
# 2008-12-05
# o Created from justCRLMMv2() of oligo v1.7.3.
############################################################################
