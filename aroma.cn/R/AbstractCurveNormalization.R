###########################################################################/**
# @RdocClass AbstractCurveNormalization
#
# @title "The AbstractCurveNormalization class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#  \item{dataSet}{An @see "aroma.core::AromaUnitTotalCnBinarySet" of
#     "test" samples to be normalized.}
#  \item{targetSet}{An @see "aroma.core::AromaUnitTotalCnBinarySet" of
#     paired target samples.}
#  \item{subsetToFit}{The subset of loci to be used to fit the
#    normalization functions.
#    If @NULL, loci on chromosomes 1-22 are used, but not on ChrX and ChrY.
#  }
#  \item{tags}{(Optional) Sets the tags for the output data sets.}
#  \item{copyTarget}{If @TRUE, target arrays are copied to the output
#     data set, otherwise not.}
#  \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("AbstractCurveNormalization", function(dataSet=NULL, targetSet=NULL, subsetToFit=NULL, tags="*", copyTarget=TRUE, ...) {
  # Validate arguments
  if (!is.null(dataSet)) {
    # Arguments 'dataSet' and 'targetSet'
    dsList <- list(dataSet=dataSet, targetSet=targetSet);
    className <- "AromaUnitTotalCnBinarySet";
    for (kk in seq_along(dsList)) {
      key <- names(dsList)[kk];
      ds <- dsList[[kk]];
      ds <- Arguments$getInstanceOf(ds, className, .name=key);
    } # for (kk ...)

    # Assert that each data set contains the same number of files
    for (jj in 1:(length(dsList)-1)) {
      keyJJ <- names(dsList)[jj];
      dsJJ <- dsList[[jj]];
      nJJ <- length(dsJJ);
      chipTypeJJ <- getChipType(dsJJ);
      for (kk in (jj+1):length(dsList)) {
        keyKK <- names(dsList)[kk];
        dsKK <- dsList[[kk]];
        nKK <- length(dsKK);
        chipTypeKK <- getChipType(dsKK);

        # Assert that each data set contains the same number of files
        if (nKK != nJJ) {
          throw(sprintf("The number of files in '%s' and '%s' does not match: %s != %s", keyKK, keyJJ, nKK, nJJ));
        }

        # Assert that each data set is for the same chip type
        if (chipTypeKK != chipTypeJJ) {
          throw(sprintf("The chip types for '%s' and '%s' does not match: %s != %s", keyKK, keyJJ, chipTypeKK, chipTypeJJ));
        }
      } # for (kk ...)
    } # for (jj ...)

    # Assert that the UGP file exists
    ugp <- getAromaUgpFile(dataSet);

    # Argument 'subsetToFit':
    if (is.null(subsetToFit)) {
    } else if (is.character(subsetToFit)) {
      throw("Yet not implemented: Argument 'subsetToFit' is of type character.");
    } else {
      subsetToFit <- Arguments$getIndices(subsetToFit, max=nbrOfUnits(ugp));
    }
  } # if (!is.null(dataSet))

  # Argument 'copyTarget':
  copyTarget <- Arguments$getLogical(copyTarget);

  # Arguments '...':
  args <- list(...);
  if (length(args) > 0) {
    argsStr <- paste(names(args), collapse=", ");
    throw("Unknown arguments: ", argsStr);
  }

  this <- extend(Object(...), "AbstractCurveNormalization",
    .dataSet = dataSet,
    .targetSet = targetSet,
    .subsetToFit = subsetToFit,
    .copyTarget = copyTarget
  );

  if (!is.null(dataSet)) {
    setTags(this, tags);
  }

  this;
})


setMethodS3("as.character", "AbstractCurveNormalization", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  dsList <- getDataSets(this);
  s <- c(s, sprintf("Data sets (%d):", length(dsList)));
  for (kk in seq_along(dsList)) {
    ds <- dsList[[kk]];
    s <- c(s, sprintf("<%s>:", capitalize(names(dsList)[kk])));
    s <- c(s, as.character(ds));
  }
  GenericSummary(s);
}, protected=TRUE)


setMethodS3("getAsteriskTags", "AbstractCurveNormalization", function(this, ...) {
  # Create a default asterisk tags for any class by extracting all
  # capital letters and pasting them together, e.g. AbcDefGhi => ADG.
  name <- class(this)[1];

  # Remove any 'Model' suffixes
  name <- gsub("Model$", "", name);

  name <- capitalize(name);

  # Vectorize
  name <- strsplit(name, split="")[[1]];

  # Identify upper case
  name <- name[(toupper(name) == name)];

  # Paste
  name <- paste(name, collapse="");

  tags <- name;

  tags;
}, protected=TRUE)



setMethodS3("getName", "AbstractCurveNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  getName(ds);
})

setMethodS3("getSubsetToFit", "AbstractCurveNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  units <- this$.subsetToFit;
  if (is.null(units)) {
    verbose && enter(verbose, "Identify subset of units for fitting the normalization function");

    verbose && enter(verbose, "Retrieving the UGP file");
    ds <- getInputDataSet(this);
    ugp <- getAromaUgpFile(ds);
    verbose && print(verbose, ugp);
    verbose && exit(verbose);

    verbose && enter(verbose, "Querying UGP for units on chromosomes of interest");
    chromosomes <- 1:22;
    verbose && cat(verbose, "Chromosomes to fit: ",
                                             seqToHumanReadable(chromosomes));
    units <- sapply(chromosomes, FUN=function(cc) {
      getUnitsOnChromosome(ugp, cc);
    });
    units <- unlist(units, use.names=FALSE);
    units <- unique(units);
    units <- sort(units);
    verbose && str(verbose, units);
    verbose && exit(verbose);

    this$.subsetToFit <- units;

    verbose && exit(verbose);
  }

  units;
}, protected=TRUE)



setMethodS3("getTags", "AbstractCurveNormalization", function(this, collapse=NULL, ...) {
  # "Pass down" tags from input data set
  ds <- getInputDataSet(this);
  tags <- getTags(ds, collapse=collapse);

  # Get class-specific tags
  tags <- c(tags, this$.tags);

  # Update default tags
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }

  if (length(tags) == 0)
    tags <- NULL;

  tags;
})


setMethodS3("setTags", "AbstractCurveNormalization", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }

  this$.tags <- tags;

  invisible(this);
})


setMethodS3("getFullName", "AbstractCurveNormalization", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getDataSets", "AbstractCurveNormalization", function(this, ...) {
  list(test=this$.dataSet, target=this$.targetSet);
}, protected=TRUE)

setMethodS3("getInputDataSet", "AbstractCurveNormalization", function(this, ...) {
  this$.dataSet;
})

setMethodS3("getTargetDataSet", "AbstractCurveNormalization", function(this, ...) {
  this$.targetSet;
})

setMethodS3("getRootPath", "AbstractCurveNormalization", function(this, ...) {
  "totalAndFracBData";
}, protected=TRUE)


setMethodS3("getPath", "AbstractCurveNormalization", function(this, create=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Chip type
  ds <- getInputDataSet(this);
  chipType <- getChipType(ds, fullname=FALSE);

  # The full path
  path <- filePath(rootPath, fullname, chipType);

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


setMethodS3("nbrOfFiles", "AbstractCurveNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  length(ds);
}, protected=TRUE)


setMethodS3("getOutputDataSet", "AbstractCurveNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting output data set");
  ds <- getInputDataSet(this);
  path <- getPath(this);
  verbose && cat(verbose, "Output path: ", path);
  res <- byPath(ds, path=path, ..., verbose=less(verbose, 20));

  verbose && enter(verbose, "Keeping output data files matching input data files");
  # Identify output data files that match the input data files
  fullnames <- getFullNames(ds);
  df <- getFile(ds, 1);
  translator <- getFullNameTranslator(df);
  setFullNamesTranslator(res, translator);
  fullnamesOut <- getFullNames(res);
  idxs <- match(fullnames, fullnamesOut);
  verbose && str(verbose, idxs);
  if (anyMissing(idxs)) {
    throw("Should not happen.");
  }
  verbose && cat(verbose, "Number of files dropped: ", length(res) - length(idxs));
  verbose && cat(verbose, "Number of files kept: ", length(idxs));
  res <- extract(res, idxs);
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
})


setMethodS3("getPairedDataSet", "AbstractCurveNormalization", function(this, array, ..., verbose=FALSE) {
  ds <- getInputDataSet(this);
  nbrOfArrays <- length(ds);

  # Argument 'array':
  array <- Arguments$getIndex(array, max=nbrOfArrays);

  df <- getFile(ds, array);
  name <- getName(df);

  verbose && enter(verbose, sprintf("Extracting paired data set for array %d ('%s') of %d", array, name, nbrOfArrays));
  dsT <- getTargetDataSet(this);
  dfT <- getFile(dsT, array);

  dsPair <- newInstance(ds, list(dfT, df));
  verbose && cat(verbose, "Pair:");
  verbose && print(verbose, dsPair);

  verbose && exit(verbose);

  dsPair;
}, protected=TRUE)



setMethodS3("fitOne", "AbstractCurveNormalization", abstract=TRUE, protected=TRUE);
setMethodS3("backtransformOne", "AbstractCurveNormalization", abstract=TRUE, protected=TRUE);


setMethodS3("process", "AbstractCurveNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Paired (x,y)-curve normalization");
  nbrOfFiles <- nbrOfFiles(this);
  verbose && cat(verbose, "Number of arrays: ", nbrOfFiles);

  copyTarget <- this$.copyTarget;
  verbose && cat(verbose, "Copying target arrays: ", copyTarget);

  ds <- getInputDataSet(this);
  chipType <- getChipType(ds, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);

  outPath <- getPath(this);
  verbose && cat(verbose, "Output path: ", outPath);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Optional function h() and g() for transforming and backtransform
  # signals.  Typically x = g(h(x)), although maybe only for positive
  # values, e.g. h(x) = log2(x) and g(y) = 2^x.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hFcn <- this$.transformFcn;
  gFcn <- this$.untransformFcn;

  # Sanity check (none or both functions must be specified)
  if (is.function(hFcn) || is.function(gFcn)) {
    if (!is.function(hFcn) || !is.function(gFcn)) {
      throw("Either none or both h(.) and g(.) functions must be given.");
    }
  }


  for (kk in seq_len(nbrOfFiles)) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Extract array pair
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dsPair <- getPairedDataSet(this, array=kk, verbose=less(verbose,5));
    verbose && cat(verbose, "Pair:");
    verbose && print(verbose, dsPair);
    verbose && cat(verbose, "Fullnames of pair:");
    verbose && print(verbose, getFullNames(dsPair));

    df <- getFile(dsPair, 2);
    fullname <- getFullName(df);
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d",
                                                    kk, fullname, nbrOfFiles));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Copy target file?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (copyTarget) {
      verbose && enter(verbose, "Copying target data file");
      dfT <- getFile(dsPair, 1);
      # Output file
      fullname <- getFullName(dfT);
      ext <- getFilenameExtension(dfT);
      filename <- sprintf("%s.%s", fullname, ext);
      pathname <- Arguments$getReadablePathname(filename, path=outPath, mustExist=FALSE);
      # Nothing to do?
      if (isFile(pathname)) {
        verbose && cat(verbose, "Already copied: ", pathname);
      } else {
        verbose && cat(verbose, "Output pathname: ", pathname);
        copyFile(getPathname(dfT), pathname, copy.mode=FALSE, verbose=less(verbose,50));
      }
      # Not needed anymore
      filename <- pathname <- NULL;
      verbose && exit(verbose);
    } # if (copyTarget)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Output file
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fullname <- getFullName(df);
    ext <- getFilenameExtension(df);
    filename <- sprintf("%s.%s", fullname, ext);
    pathname <- Arguments$getReadablePathname(filename, path=outPath, mustExist=FALSE);

    # Nothing to do?
    if (isFile(pathname)) {
      verbose && cat(verbose, "Already normalized: ", pathname);
      verbose && exit(verbose);
      next;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Reading data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Reading all data");
    theta <- extractMatrix(dsPair, verbose=less(verbose,5));
    nbrOfUnits <- nrow(theta);
    verbose && str(verbose, theta);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Transforming data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.function(hFcn)) {
      verbose && enter(verbose, "Transforming");
      verbose && cat(verbose, "Function y <- h(x):");
      verbose && str(verbose, hFcn);
      theta <- hFcn(theta);
      verbose && str(verbose, theta);
      verbose && exit(verbose);
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fitting
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Fitting");

    verbose && enter(verbose, "Identifying subset to fit");
    subsetToFit <- getSubsetToFit(this);
    verbose && cat(verbose, "Subset to fit:");
    verbose && str(verbose, subsetToFit);
    verbose && exit(verbose);

    verbose && enter(verbose, "Extracting subset used for fitting normalization function");
    thetaFit <- theta[subsetToFit,,drop=FALSE];
    # Not needed anymore
    subsetToFit <- NULL;
    verbose && str(verbose, thetaFit);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calling fit function");
    fit <- fitOne(this, theta=thetaFit, ..., verbose=verbose);
    # Not needed anymore
    thetaFit <- NULL;
    verbose && str(verbose, fit);
    verbose && exit(verbose);

#    verbose && enter(verbose, "Saving fit");
#    filename <- sprintf("%s,fit.Rbin", fullname);
#    pathnameF <- Arguments$getWritablePathname(filename, path=outPath, mustNotExist=FALSE);
#    verbose && cat(verbose, "Fit pathname: ", pathnameF);
#    saveObject(fit, file=pathnameF);
#    verbose && exit(verbose);

    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Normalizing
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Normalizing all data (using fitted normalization function)");
    verbose && cat(verbose, "theta:");
    verbose && str(verbose, theta);
    verbose && cat(verbose, "fit:");
    verbose && str(verbose, fit);

    # Allocate normalized signals
    thetaN <- backtransformOne(this, theta=theta, fit=fit, targetDimension=1);
    # Not needed anymore
    theta <- fit <- NULL;
    thetaN <- thetaN[,2,drop=TRUE];
    verbose && str(verbose, thetaN);
    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Back-transforming data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.function(gFcn)) {
      verbose && enter(verbose, "Backtransforming");
      verbose && cat(verbose, "Function x* <- g(y*) where g(h(x)) = x (should be):");
      verbose && str(verbose, gFcn);
      thetaN <- gFcn(thetaN);
      verbose && str(verbose, thetaN);
      verbose && exit(verbose);
    }

    # Sanity check
    stopifnot(length(thetaN) == nbrOfUnits);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Storing
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing normalized data");
    verbose && cat(verbose, "Output pathname: ", pathname);

    verbose && enter(verbose, "Allocating to temporary file");
    pathnameT <- sprintf("%s.tmp", pathname);
    file.copy(getPathname(df), pathnameT);
    dfN <- newInstance(df, pathnameT);
    srcFiles <- lapply(dsPair, function(df) {
      list(
        filename = getFilename(df),
        filesize = getFileSize(df),
        checksum = getChecksum(df)
      )
    });
    footer <- readFooter(dfN);
    footer$srcFiles <- srcFiles;
    writeFooter(dfN, footer);
    # Not needed anymore
    srcFiles <- footer <- NULL;
    verbose && exit(verbose);

    verbose && enter(verbose, "Writing to temporary file");
    dfN[,1] <- thetaN;
    # Not needed anymore
    thetaN <- NULL;
    verbose && exit(verbose);

    # Renaming
    verbose && enter(verbose, "Renaming temporary file");
    file.rename(pathnameT, pathname);
    if (!isFile(pathname)) {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    verbose && exit(verbose);
    # Not needed anymore
    dfN <- pathnameT <- pathname <- NULL;
    verbose && exit(verbose);

    # Not needed anymore
    df <- NULL;

    verbose && exit(verbose);
  } # for (kk ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return output data set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- getOutputDataSet(this, verbose=less(verbose, 1));

  verbose && exit(verbose);

  invisible(res);
})

############################################################################
# HISTORY:
# 2014-09-04
# o ROBUSTNESS: It could be that process() for AbstractCurveNormalization
#   would generate an error due to read-only permissions introduced
#   by copying the target file without resetting the file permissions.
# 2012-04-16
# o DOCUMENTATION: Removed reference to aroma.light::fitPrincipalCurve().
# 2010-01-05
# o BUG FIX: getOutputDataSet() of AbstractCurveNormalization returned all
#   files, not just the ones matching the input set.
# o Added support for transform/untransform functions h(.) and g(.) to
#   AbstractCurveNormalization, which allows us to fit say on the log
#   scale, e.g. h(x)=log2(x), g(y)=2^y.
# 2009-07-15
# o Created.
############################################################################
