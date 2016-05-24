setMethodS3("doSegmentByPairedPSCBS", "AromaUnitPscnBinarySet", function(pscnT, pscnN, ..., ascn=c("classic", "paired"), avgDH=c("median", "mean"), tbn=TRUE, B=1000L, cache=TRUE, subset=NULL, verbose=FALSE) {
  # Assert packages
  use("R.filesets (>= 2.6.0)")

  pairedAlleleSpecificCopyNumbers <- .pairedAlleleSpecificCopyNumbers

  use("PSCBS (>= 0.43.0)")
  segmentByPairedPSCBS <- PSCBS::segmentByPairedPSCBS
  dropSegmentationOutliers <- PSCBS::dropSegmentationOutliers
  bootstrapSegmentsAndChangepoints <- PSCBS::bootstrapSegmentsAndChangepoints
  bootstrapTCNandDHByRegion <- PSCBS::bootstrapTCNandDHByRegion


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pscnT':
  pscnT <- Arguments$getInstanceOf(pscnT, "AromaUnitPscnBinarySet")
  ugpT <- getAromaUgpFile(pscnT)

  # Argument 'pscnN':
  if (inherits(pscnN, "AromaUnitPscnBinaryFile")) {
    pscnN <- newInstance(pscnT, pscnN);
  }
  pscnN <- Arguments$getInstanceOf(pscnN, "AromaUnitPscnBinarySet")
  ugpN <- getAromaUgpFile(pscnN)
  if (!equals(ugpN, ugpT)) {
    throw("Argument 'pscnT' and 'pscnN' appears to be for different chip types: ", getFullName(ugpT), " != ", , getFullName(ugpT))
  }
  # Use a single match normal for all tumors?
  if (length(pscnN) == 1L) {
    pscnN <- pscnN[rep(1L, times=length(pscnT))]
  }
  stopifnot(length(pscnN) == length(pscnT))

  # Argument 'ascn':
  ascn <- match.arg(ascn)
  if (ascn == "paired") {
    use("aroma.light (>= 1.34.0)")
  }

  # Argument 'avgDH':
  avgDH <- match.arg(avgDH)

  # Argument 'tbn':
  tbn <- Arguments$getLogical(tbn)

  # Argument 'B':
  B <- Arguments$getInteger(B, range=c(1,Inf))

  # Argument 'cache':
  cache <- Arguments$getLogical(cache)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Paired tumor-normal PSCBS segmentation")
  verbose && print(verbose, "Tumor data set:")
  verbose && print(verbose, pscnT)
  verbose && print(verbose, "Matched normal data set:")
  verbose && print(verbose, pscnN)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify sample pairs to process
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying samples to be processed")
  pairNames <- mapply(pscnT, pscnN, FUN=function(dfT, dfN) {
    sprintf("%s_vs_%s", getName(dfT), getName(dfN))
  })
  sampleNames <- pairNames
  ## FIXME
  ## sampleNames <- sprintf("HC1143,TvsN,%s", sampleNames)

  # Output pathnames
  ascnTag <- switch(ascn, classic=NULL, paired="pASCN");
  subsetTag <- sprintf("subset=%g", subset)
  tags <- c(ascnTag, subsetTag)
  datasetF <- paste(c(getFullName(pscnT), tags), collapse=",")
  path <- file.path("pscbsData", datasetF, getChipType(pscnN))
  path <- Arguments$getWritablePath(path)
  filenames <- sprintf("%s,PairedPSCBS.xdr", sampleNames)
  pathnames <- file.path(path, filenames)

  # Non-processed samples
  todo <- which(!sapply(pathnames, FUN=isFile))

  # Nothing to do?
  if (length(todo) == 0L) {
    verbose && cat(verbose, "Already processed. Skipping.")
    verbose && exit(verbose)
    verbose && exit(verbose)
    res <- lapply(pathnames, FUN=PairedPSCBSFile, mustExist=FALSE)
    res <- PairedPSCBSFileSet(res)
    return(res)
  }

  # Tumor-normal pairs to process
  pathnames <- pathnames[todo]
  pairNames <- pairNames[todo]
  pscnT <- pscnT[todo]
  pscnN <- pscnN[todo]
  stopifnot(length(pscnT) == length(pscnN))

  verbose && print(verbose, "Tumor data set:")
  verbose && print(verbose, pscnT)
  verbose && print(verbose, "Matched normal data set:")
  verbose && print(verbose, pscnN)
  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving annotation data")
  ugp <- getAromaUgpFile(pscnT)
  dataA <- ugp[,1:2]
  colnames(dataA)[2] <- "x"
  verbose && str(verbose, dataA)
  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Process via dsApplyInPairs()
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- dsApplyInPairs(pscnT, pscnN, FUN=function(dsPair, dataA, ..., ascn=c("classic", "paired"), tbn=TRUE, B=NULL, cache=FALSE, subset=NULL, seed=NULL, verbose=FALSE) {
    use("R.utils (>= 1.34.0)", verbose=TRUE)
    use("aroma.light", verbose=TRUE)
    use("PSCBS (>= 0.43.0)", verbose=TRUE)
    use("aroma.cn (>= 1.5.5)", verbose=TRUE)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'dsPair':
    dsPair <- AromaUnitPscnBinarySet(dsPair)

    # Argument 'dataA':
    dataA <- Arguments$getInstanceOf(dataA, "data.frame");

    # Argument 'ascn':
    ascn <- match.arg(ascn)
    if (ascn == "paired") {
      use("aroma.light (>= 1.34.0)", verbose=TRUE)
    }

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose)
    if (verbose) {
      pushState(verbose)
      on.exit(popState(verbose))
    }


    verbose && enter(verbose, "Paired tumor-normal PSCBS segmentation")

    pairName <- paste(sapply(dsPair, FUN=getName), collapse="_vs_")
    sampleName <- pairName
    verbose && cat(verbose, "Sample name: ", sampleName)

    ascnTag <- switch(ascn, classic=NULL, paired="pASCN")
    subsetTag <- sprintf("subset=%g", subset)
    tags <- c(ascnTag, subsetTag)
    datasetF <- paste(c(getFullName(pscnT), tags), collapse=",")
    path <- file.path("pscbsData", datasetF, getChipType(dsPair))
    path <- Arguments$getWritablePath(path)

    filename <- sprintf("%s,PairedPSCBS.xdr", sampleName)
    pathname <- file.path(path, filename)
    verbose && cat(verbose, "Output pathname: ", pathname)

    # Nothing todo?
    if (isFile(pathname)) {
      verbose && cat(verbose, "Already processed. Skipping.")
      verbose && exit(verbose)
      res <- PairedPSCBSFile(pathname)
      return(res)
    }


    # Extract tumor-normal pair data
    dfT <- dsPair[[1]]
    dfN <- dsPair[[2]]
    dataII <- cbind(dfT[,1:2], dfN[,1:2])
    colnames(dataII) <- c("thetaT", "betaT", "thetaN", "betaN")
    verbose && cat(verbose, "Loaded tumor-normal PSCN data:")
    verbose && str(verbose, dataII)

    # Calculate paired ASCNs?
    if (ascn == "paired") {
      verbose && enter(verbose, "Calculating paired ASCNs")
      dataII <- do.call(pairedAlleleSpecificCopyNumbers, args=dataII)
      verbose && cat(verbose, "Tumor-normal PSCN data:")
      verbose && str(verbose, dataII)
      verbose && exit(verbose)

      # Disable TumorBoost
      tbn <- FALSE
    }

    data <- cbind(dataA, dataII)
    dataA <- dataII <- NULL # Not needed anymore
    verbose && cat(verbose, "Annotated tumor-normal PSCN data:")
    verbose && str(verbose, data)

    # Drop outliers
    verbose && enter(verbose, "Dropping TCN outliers")
    data <- dropSegmentationOutliers(data)
#    verbose && print(verbose, head(data))
    verbose && exit(verbose)

    if (is.numeric(subset)) {
      verbose && enter(verbose, "Subsetting by sampling")
      verbose && cat(verbose, "Subset: ", subset)
      set.seed(seed)
      keep <- sample(1:nrow(data), size=subset*nrow(data))
      data <- data[keep,]
      keep <- NULL; # Not needed anymore
 #     verbose && print(verbose, head(data))
      verbose && exit(verbose)
    }

    verbose && cat(verbose, "PSCN signals:")
#    verbose && print(verbose, head(data))

    # Segment
    verbose && enter(verbose, "Segmenting")
    fit <- segmentByPairedPSCBS(data, ..., avgDH=avgDH, tbn=tbn, preserveScale=FALSE, seed=seed, verbose=verbose)
    verbose && str(verbose, fit)
    verbose && exit(verbose)

    data <- NULL # Not needed anymore

    # Bootstrap and cache?
    if (!is.null(B)) {
      verbose && enter(verbose, "Bootstrapping")
      boot <- bootstrapSegmentsAndChangepoints(fit, B=B, seed=seed, cache=cache, verbose=less(verbose, 10))
      verbose && exit(verbose)

      # Call segments
      verbose && enter(verbose, "Calling")
      fit <- bootstrapTCNandDHByRegion(fit, boot=boot)
      fit <- callROH(fit, verbose=less(verbose, 10))
      fit <- callAB(fit, verbose=less(verbose, 10))
      fit <- callLOH(fit, verbose=less(verbose, 10))
      verbose && exit(verbose)
    }

    verbose && enter(verbose, "Saving results")
    saveObject(fit, file=pathname)
    verbose && exit(verbose)

    verbose && exit(verbose)

    PairedPSCBSFile(pathname)
  }, dataA=dataA, ascn=ascn, tbn=tbn, B=B, cache=cache, subset=subset, verbose=verbose, ...)

  # Load available ones
  res <- lapply(pathnames, FUN=PairedPSCBSFile, mustExist=FALSE)
  res <- PairedPSCBSFileSet(res)

  verbose && exit(verbose)

  res
}) # doSegmentByPairedPSCBS()


setMethodS3("findLargeGaps", "AromaUgpFile", function(this, ...) {
  findLargeGaps <- PSCBS::findLargeGaps

  data <- this[,1:2];
  colnames(data)[2L] <- "x";
  findLargeGaps(data, ...);
}, protected=TRUE)


setMethodS3("findLargeGaps", "AromaUnitPscnBinarySet", function(this, ...) {
  findLargeGaps <- PSCBS::findLargeGaps

  ugp <- getAromaUgpFile(this);
  findLargeGaps(ugp, ...);
}, protected=TRUE)

##############################################################################
# HISTORY
# 2014-03-31 [HB]
# o Added findLargeGaps() for AromaUgpFile and AromaUnitPscnBinarySet.
# 2014-03-30 [HB]
# o Added doSegmentByPairedPSCBS().
# o Created.
##############################################################################
