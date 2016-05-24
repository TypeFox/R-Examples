###########################################################################/**
# @RdocClass TotalCnSmoothing
#
# @title "The abstract TotalCnSmoothing class"
#
# \description{
#  @classhierarchy
#
# }
#
# @synopsis
#
# \arguments{
#  \item{dataSet}{An @see "aroma.core::AromaUnitTotalCnBinarySet".}
#  \item{...}{Arguments passed to @see "aroma.core::AromaTransform".}
#  \item{targetUgp}{An @see "aroma.core::AromaUgpFile" specifying the
#    target loci for which smoothed copy-number are generated.}
#  \item{.reqSetClass}{(internal only)}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("TotalCnSmoothing", function(dataSet=NULL, ..., targetUgp=NULL, .reqSetClass="AromaUnitTotalCnBinarySet") {
  if (!is.null(dataSet)) {
    # Argument 'targetUgp':
    targetUgp <- Arguments$getInstanceOf(targetUgp, "AromaUgpFile");
  }

  extend(AromaTransform(dataSet=dataSet, ..., .reqSetClass=.reqSetClass), "TotalCnSmoothing",
    .targetUgp = targetUgp
  );
}, abstract=TRUE)


setMethodS3("getParameters", "TotalCnSmoothing", function(this, ...) {
  params <- NextMethod("getParameters");
  params$targetUgp <- this$.targetUgp;
  params;
}, protected=TRUE);


setMethodS3("getAsteriskTags", "TotalCnSmoothing", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Add class-specific tags

  params <- getParameters(this);
  # "Parameter" 'by'
  byTag <- grep("(b|kb|Mb)$", getTags(params$targetUgp), value=TRUE);
  if (length(byTag) > 0) {
    byTag <- sprintf("by=%s", byTag[1]);
  }
  tags <- c(tags, byTag);

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
}, protected=TRUE)



setMethodS3("getRootPath", "TotalCnSmoothing", function(this, ...) {
  "smoothCnData";
}, protected=TRUE)



setMethodS3("getTargetUgpFile", "TotalCnSmoothing", function(this, ...) {
  this$.targetUgp;
})

setMethodS3("getPath", "TotalCnSmoothing", function(this, create=TRUE, ...) {
  path <- NextMethod("getPath", create=FALSE);
  path <- dirname(path);
  targetUgp <- getTargetUgpFile(this);
  chipType <- getChipType(targetUgp, fullname=FALSE);

  # The full path
  path <- filePath(path, chipType);

  # Create path?
  if (create) {
    path <- Arguments$getWritablePath(path);
  } else {
    path <- Arguments$getReadablePath(path, mustExist=FALSE);
  }

  # Verify that it is not the same as the input path
  inPath <- getPath(getInputDataSet(this));
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  path;
}, protected=TRUE)



setMethodS3("getTargetPositions", "TotalCnSmoothing", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- getParameters(this);
  targetUgp <- params$targetUgp;

  # For now, always use all chromosomes
  chromosomes <- NULL;

  # Argument 'chromosomes':
  if (is.null(chromosomes)) {
    chromosomes <- getChromosomes(targetUgp);
  } else {
    chromosomes <- Arguments$getIndices(chromosomes);
  }
  nbrOfChromosomes <- length(chromosomes);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  res <- this$.targetPositions;
  if (!force && !is.null(res)) {
    return(res);
  }

  verbose && enter(verbose, "Identifying all target positions");
  verbose && cat(verbose, "Chromosomes:");
  verbose && print(verbose, chromosomes);

  verbose && print(verbose, targetUgp);

  res <- list();
  for (cc in chromosomes) {
    chrTag <- sprintf("Chr%02d", cc);
    verbose && enter(verbose, sprintf("Chromosome %d ('%s') of %d",
                                             cc, chrTag, nbrOfChromosomes));
    verbose && cat(verbose, "Target positions:");
    units <- getUnitsOnChromosome(targetUgp, chromosome=cc);
    xOut <- getPositions(targetUgp, units=units);
    verbose && str(verbose, xOut);
    res[[chrTag]] <- list(chromosome=cc, units=units, xOut=xOut);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  this$.targetPositions <- res;

  res;
}, protected=TRUE)



setMethodS3("smoothRawCopyNumbers", "TotalCnSmoothing", abstract=TRUE, protected=TRUE);


setMethodS3("getOutputFileExtension", "TotalCnSmoothing", function(this, ...) {
  ds <- getInputDataSet(this);
  df <- getFile(ds, 1);
  ext <- getFilenameExtension(df);
  sprintf(".%s", ext);
}, protected=TRUE)


setMethodS3("getOutputFileSetClass", "TotalCnSmoothing", function(this, ...) {
  ds <- getInputDataSet(this);
  className <- class(ds)[1];
  Class$forName(className);
}, protected=TRUE)


setMethodS3("getOutputFileClass", "TotalCnSmoothing", function(this, ...) {
  setClass <- getOutputFileSetClass(this, ...);
  setInstance <- newInstance(setClass);
  className <- getFileClass(setInstance);
  Class$forName(className);
}, protected=TRUE)


setMethodS3("getOutputDataSet0", "TotalCnSmoothing", function(this, pattern=NULL, className=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'pattern':
  if (!is.null(pattern)) {
    pattern <- Arguments$getRegularExpression(pattern=pattern);
  }


  verbose && enter(verbose, "Retrieving existing set of output files");
  ds <- getInputDataSet(this);
  outPath <- getPath(this);
  if (is.null(className)) {
    clazz <- getOutputFileSetClass(this);
    className <- getName(clazz);
  }
  verbose && cat(verbose, "Class: ", className);

  path <- getPath(this);
  verbose && cat(verbose, "Path: ", path);

  if (is.null(pattern)) {
    # Default filename pattern find non-private (no dot prefix) files with
    # the same file name extension as the input data set.
    fileExt <- getOutputFileExtension(this);
    fileExt <- c(fileExt, tolower(fileExt), toupper(fileExt));
    fileExt <- sprintf("(%s)", paste(unique(fileExt), collapse="|"));
    verbose && cat(verbose, "Expected file extensions: ", fileExt);
    pattern <- sprintf("^[^.].*%s$", fileExt);
  }
  verbose && cat(verbose, "Pattern: ", pattern);

  verbose && enter(verbose, sprintf("Calling %s$forName()", className));
  clazz <- Class$forName(className);
  args <- list(path=path, pattern=pattern, ...);
  verbose && str(verbose, args);
  args$verbose <- less(verbose);
  staticMethod <- clazz$byPath;
  dsOut <- do.call("staticMethod", args=args);
  # Not needed anymore
  staticMethod <- args <- NULL;
  verbose && exit(verbose);

  verbose && exit(verbose);

  dsOut;
}, protected=TRUE)



setMethodS3("process", "TotalCnSmoothing", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  if (isDone(this)) {
    dsOut <- getOutputDataSet(this);
    return(invisible(dsOut));
  }

  verbose && enter(verbose, "Smoothing copy-number towards set of target loci");

  params <- getParameters(this);

  verbose && print(verbose, "Input data set:");
  ds <- getInputDataSet(this);
  verbose && print(verbose, ds);

  verbose && enter(verbose, "Identifying all target positions");
  targetList <- getTargetPositions(this, ...);
  nbrOfChromosomes <- length(targetList);
  verbose && str(verbose, targetList);

  targetUgp <- params$targetUgp;
  platform <- getPlatform(targetUgp);
  chipType <- getChipType(targetUgp);
  nbrOfUnits <- nbrOfUnits(targetUgp);
  # Not needed anymore
  targetUgp <- NULL;
  verbose && cat(verbose, "Total number of target units:", nbrOfUnits);
  verbose && exit(verbose);

  # Get Class object for the output files
  clazz <- getOutputFileClass(this);

  # Get the filename extension for output files
  ext <- getOutputFileExtension(this);


  nbrOfArrays <- length(ds);
  for (kk in seq_along(ds)) {
    df <- getFile(ds, kk);
    verbose && enter(verbose, sprintf("Array %d ('%s') of %d",
                                            kk, getName(df), nbrOfArrays));

    path <- getPath(this);
    fullname <- getFullName(df);
    filename <- sprintf("%s%s", fullname, ext);
    pathname <- Arguments$getReadablePathname(filename, path=path,
                                                         mustExist=FALSE);
    verbose && cat(verbose, "Output pathname: ", pathname);

    if (isFile(pathname)) {
      dfOut <- newInstance(clazz, filename=pathname);
      if (nbrOfUnits != nbrOfUnits(dfOut)) {
        throw("The number of units in existing output file does not match the number of units in the output file: ", nbrOfUnits, " != ", nbrOfUnits(dfOut));
      }
      verbose && cat(verbose, "Skipping already existing output file.");
      verbose && exit(verbose);
      next;
    }

    verbose && print(verbose, df);

    # Preallocate vector
    M <- rep(as.double(NA), times=nbrOfUnits);

    verbose && enter(verbose, "Reading and smoothing input data");
    for (cc in seq_along(targetList)) {
      target <- targetList[[cc]];
      chromosome <- target$chromosome;
      chrTag <- sprintf("Chr%02d", chromosome);

      verbose && enter(verbose, sprintf("Chromosome %d ('%s') of %d",
                                               cc, chrTag, nbrOfChromosomes));
      verbose && cat(verbose, "Extracting raw CNs:");
      rawCNs <- extractRawCopyNumbers(df, chromosome=chromosome,
                                                  verbose=less(verbose, 10));
      verbose && print(verbose, rawCNs);
      verbose && summary(verbose, rawCNs);

      verbose && cat(verbose, "Smoothing CNs:");
      verbose && cat(verbose, "Target positions:");
      verbose && str(verbose, target$xOut);

      smoothCNs <- smoothRawCopyNumbers(this, rawCNs=rawCNs,
                                        target=target, verbose=verbose);

      verbose && print(verbose, smoothCNs);
      verbose && summary(verbose, smoothCNs);

      M[target$units] <- getSignals(smoothCNs);
      verbose && exit(verbose);
    } # for (cc ...)

    verbose && cat(verbose, "Smoothed CNs across all chromosomes:");
    verbose && str(verbose, M);
    verbose && summary(verbose, M);
    verbose && printf(verbose, "Missing values: %d (%.1f%%) out of %d\n",
                   sum(is.na(M)), 100*sum(is.na(M))/nbrOfUnits, nbrOfUnits);
    verbose && exit(verbose);

    verbose && enter(verbose, "Storing smoothed data");
    verbose && cat(verbose, "Pathname: ", pathname);

    params2 <- params;
    params2[["targetUgp"]] <- NULL;
    footer <- list(
      sourceDataFile=list(
        fullname=getFullName(df),
        platform=getPlatform(df),
        chipType=getChipType(df),
        checksum=getChecksum(df)
      ), parameters=list(
        targetUgp=list(
          fullname=getFullName(params$targetUgp),
          platform=getPlatform(params$targetUgp),
          chipType=getChipType(params$targetUgp),
          checksum=getChecksum(params$targetUgp)
        ),
        params=params2
      )
    );

    # Write to a temporary file
    pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

    dfOut <- clazz$allocate(filename=pathnameT, nbrOfRows=nbrOfUnits,
                            platform=platform, chipType=chipType,
                            footer=footer, verbose=less(verbose, 50));

    dfOut[,1] <- M;
    # Not needed anymore
    M <- NULL;

    # Renaming temporary file
    pathname <- popTemporaryFile(pathnameT, verbose=verbose);

    verbose && exit(verbose); # Storing

    verbose && exit(verbose);
  } # for (kk ...)

  verbose && exit(verbose);

  dsOut <- getOutputDataSet(this);
  invisible(dsOut);
})



setMethodS3("getOutputFiles", "TotalCnSmoothing", function(this, ...) {
  NextMethod("getOutputFiles", pattern=".*[.]asb$");
}, protected=TRUE)



############################################################################
# HISTORY:
# 2012-10-11
# o Added custom getOutputDataSet0() for TotalCnSmoothing that utilizes
#   getOutputFileExtension().
# o Added getOutputFileClass() and getOutputFileExtension() for
#   TotalCnSmoothing.
# 2011-12-15
# o Moved argument 'bandwidth' to TotalCnKernelSmoothing.
# o ROBUSTNESS: Now process() of TotalCnSmoothing write output atomically.
# 2009-05-05
# o BUG FIX: process() of TotalCnSmoothing would not "recognize" fullname
#   translators, that is, the output filenames were always identical to
#   the input ones.
# 2009-05-04
# o BUG FIX: Added missing argument 'verbose' in getTargetPositions() of
#   TotalCnSmoothing.  This caused unwanted verbose output in some cases.
# 2009-02-08
# o Now the root path is smoothCnData/ and no longer cnData/.
# o Any subclass must implement smoothRawCopyNumbers().
# o Made TotalCnSmoothing an abstract class, cf. TotalCnKernelSmoothing.
# 2009-01-26
# o Adopted to the new AromaUnitTotalCnBinarySet.
# o Added Rdoc comments.
# 2008-05-23
# o Created.
############################################################################
