setConstructorS3("SegmentationDataSet", function(files=NULL, ...) {
  extend(GenericDataFileSet(files=files, ...), "SegmentationDataSet");
})



setMethodS3("getDefaultFullName", "SegmentationDataSet", function(this, ...) {
  path <- getPath(this);
  path <- getParent(path);
  basename(path);
}, protected=TRUE)


setMethodS3("getChipType", "SegmentationDataSet", function(this, ...) {
  path <- getPath(this);
  basename(path);
})


setMethodS3("as.character", "SegmentationDataSet", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("Chip type(s): %s", getChipType(this)));

  nf <- length(this);

  snames <- getSampleNames(this);
  ns <- length(snames);
  s <- c(s, sprintf("Sample names: %s [%d]", hpaste(snames), ns));

  chrs <- getChromosomes(this);
  nc <- length(chrs);
  s <- c(s, sprintf("Chromosomes: %s [%d]", hpaste(chrs), nc));

  s <- c(s, sprintf("Number of \"missing\" files: %d [%d*%d-%d=%d]",
                     (nc*ns-nf), nc, ns, (nc*ns-nf), nf));

  rnames <- getReferenceNames(this);
  nr <- length(rnames);
  s <- c(s, sprintf("Reference names: %s [%d]", hpaste(rnames), nr));
  s;
}, protected=TRUE)



setMethodS3("byPath", "SegmentationDataSet", abstract=TRUE, static=TRUE, protected=TRUE);



setMethodS3("getReferenceNames", "SegmentationDataSet", function(this, ..., force=FALSE) {
  referenceNames <- this$.referenceNames;
  if (force || is.null(referenceNames)) {
    referenceNames <- c();
    for (kk in seq_along(this)) {
      df <- this[[kk]];
      referenceNames <- c(referenceNames, getReferenceName(df));
    }
    referenceNames <- unique(referenceNames);
    referenceNames <- sort(referenceNames);

    this$.referenceNames <- referenceNames;
  }

  referenceNames;
})


setMethodS3("getSampleNames", "SegmentationDataSet", function(this, ..., force=FALSE) {
  sampleNames <- this$.sampleNames;
  if (force || is.null(sampleNames)) {
    sampleNames <- c();
    for (kk in seq_along(this)) {
      df <- this[[kk]];
      sampleNames <- c(sampleNames, getSampleName(df));
    }
    sampleNames <- unique(sampleNames);
    sampleNames <- sort(sampleNames);

    this$.sampleNames <- sampleNames;
  }

  sampleNames;
})


setMethodS3("getChromosomes", "SegmentationDataSet", function(this, ..., force=FALSE) {
  chromosomes <- this$.chromosomes;

  if (force || is.null(chromosomes)) {
    chromosomes <- logical(100);
    for (kk in seq_along(this)) {
      df <- this[[kk]];
      chromosome <- getChromosome(df);
      chromosomes[chromosome] <- TRUE;
    }
    chromosomes <- which(chromosomes);
    this$.chromosomes <- chromosomes;
  }

  chromosomes;
})


setMethodS3("extractByReferenceName", "SegmentationDataSet", function(this, referenceName, ...) {
  # Argument 'referenceName':
  referenceName <- Arguments$getCharacter(referenceName);

  names <- sapply(this, getReferenceName);
  keep <- is.element(names, referenceName);
  idxs <- which(keep);
  extract(this, idxs, ...);
}, protected=TRUE)


setMethodS3("extractBySampleName", "SegmentationDataSet", function(this, sampleName, ...) {
  # Argument 'sampleName':
  sampleName <- Arguments$getCharacter(sampleName);

  names <- sapply(this, getSampleName);
  keep <- is.element(names, sampleName);
  idxs <- which(keep);
  extract(this, idxs, ...);
}, protected=TRUE)




setMethodS3("extractCopyNumberRegions", "SegmentationDataSet", function(this, ..., mergeBySample=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extracting segmentations");
  res <- list();
  for (kk in seq_along(this)) {
    df <- this[[kk]];
    verbose && enter(verbose, sprintf("Segmentation #%d of %d", kk, length(this)));
    resKK <- extractCopyNumberRegions(df, ...);
    res[[kk]] <- resKK;
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  if (mergeBySample) {
    verbose && enter(verbose, "Merge segmentations by sample");
    sampleNames <- sapply(res, FUN=function(cnr) cnr$sampleName);
    verbose && str(verbose, sampleNames);
    uSampleNames <- unique(sampleNames);
    resT <- list();
    for (ii in seq_along(uSampleNames)) {
      sampleName <- uSampleNames[ii];
      verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", ii, sampleName, length(uSampleNames)));
      idxs <- which(sampleNames == sampleName);
      verbose && str(verbose, idxs);

      resII <- res[idxs];
      resII <- Reduce(append, resII);
      resT[[sampleName]] <- resII;
      verbose && exit(verbose);
    }
    res <- resT;
    verbose && exit(verbose);
  }

  res;
}, protected=TRUE);




############################################################################
# HISTORY:
# 2011-04-03
# o Updated as.character() to utilize hpaste().
# 2010-08-05
# o Created.
############################################################################
