setMethodS3("calculateAverageColumnAcrossFiles", "GenericTabularFileSet", function(this, method=c("mean", "median"), na.rm=TRUE, ..., ram=NULL, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'this':
  nbrOfFiles <- length(this);
  if (nbrOfFiles == 0L) {
    throw("Cannot calculate average across data files. No data files in data set: ", getFullName(this));
  }

  # Argument 'method':
  method <- match.arg(method);

  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Estimating average column across arrays");
  verbose && cat(verbose, "Averaging method: ", method);
  verbose && cat(verbose, "na.rm: ", na.rm);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dataSet <- getFullName(this);
  chipType <- getChipType(this);
  key <- list(dataSet=dataSet, chipType=chipType, fullnames=getFullNames(this), method=method, na.rm=na.rm);
  dirs <- c("aroma.affymetrix", dataSet, chipType);
  if (!force) {
    res <- loadCache(key=key, dirs=dirs);
    if (!is.null(res)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(res);
    }
  }


  if (method == "mean") {
    rowFcn <- rowMeans;
  } else if (method == "median") {
    rowFcn <- rowMedians;
  }

  verbose && cat(verbose, "Number of files: ", nbrOfFiles);

  df <- getOneFile(this);
  nbrOfRows <- nbrOfRows(df);
  units <- seq_len(nbrOfRows);
  nbrOfRows <- length(units);
  verbose && cat(verbose, "Number of rows: ", nbrOfRows);

  nbrOfRowsPerChunk <- ram*50e6;
  chunkSize <- ceiling(nbrOfRowsPerChunk/nbrOfFiles);

  verbose && cat(verbose, "Number of rows per chunk: ", chunkSize);

  unitChunks <- splitInChunks(units, chunkSize=chunkSize);
  # Not needed anymore
  units <- NULL;

  res <- lapply(unitChunks, FUN=function(units) {
    data <- extractMatrix(this, units=units, ...);
    est <- rowFcn(data, na.rm=na.rm);
    # Sanity check
    stopifnot(length(est) == length(units));
    est;
  });
  # Not needed anymore
  unitChunks <- NULL;

  res <- unlist(res, use.names=FALSE);
  verbose && str(verbose, res);

  # Sanity check
  stopifnot(length(res) == nbrOfRows);

  saveCache(res, key=key, dirs=dirs);

  verbose && exit(verbose);

  res;
}, protected=TRUE)

############################################################################
# HISTORY:
# 2009-02-13
# o Added support for memoization.
# o Created.
############################################################################
