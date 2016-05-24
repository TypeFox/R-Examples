setMethodS3("getMaxLengthRepeats", "AromaCellSequenceFile", function(this, cells, positions=1:getProbeLength(this), bases=c("A","C","G","T"), ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, max=nbrOfCells(this));
    nbrOfCells <- length(cells);
  } else {
    nbrOfCells <- nbrOfCells(this);
  }

  # Argument 'bases':
  bases <- Arguments$getCharacters(bases);
  bases <- unique(bases);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Identifying cells with repeats");
  verbose && cat(verbose, "Nucleotide positions:");
  verbose && str(verbose, positions);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this);
  key <- list(method="getMaxRepeatLength", class=class(this)[1],
              chipType=chipType, tags=getTags(this),
              cells=cells, ...);
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="getMaxRepeatLength", chipType=chipType,
                                     tags=getTags(this), cells=cells, ...);
  }
  dirs <- c("aroma.affymetrix", chipType);
  if (!force) {
    verbose && enter(verbose, "Checking for cached results");
    res <- loadCache(key=key, dirs=dirs);
    if (!is.null(res)) {
      verbose && cat(verbose, "Found cached results");
      verbose && exit(verbose);
      verbose && exit(verbose);
      return(res);
    }
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify repeats
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Searching for repeats");
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);

  # Allocate
  startPosition <- rep(positions[1], times=nbrOfCells);
  maxRepeatLength <- rep(as.integer(1), times=nbrOfCells);
  verbose && str(verbose, startPosition);
  verbose && str(verbose, maxRepeatLength);

  verbose && enter(verbose, "Missing probes sequences");
  nok <- isMissing(this, cells=cells);
  naValue <- as.integer(NA);
  startPosition[nok] <- naValue;
  maxRepeatLength[nok] <- naValue;

  if (is.null(cells)) {
    cells <- which(!nok);
  } else {
    cells <- cells[!nok];
  }
  verbose && str(verbose, startPosition);
  verbose && str(verbose, maxRepeatLength);
  verbose && exit(verbose);

  verbose && enter(verbose, "Remaining probe sequences");
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);

  if (length(cells) > 0) {
    counts <- rep(as.integer(0), times=length(cells));
    startPositionT <- rep(as.integer(0), times=length(cells));
    maxRepeatLengthT <- rep(as.integer(1), times=length(cells));

    for (pp in seq_along(positions)) {
      verbose && enter(verbose, sprintf("Position %d of %d", pp, length(positions)));
      pos <- positions[pp];
      b1 <- readSequenceMatrix(this, cells=cells, position=pos, what="raw", drop=TRUE);

      if (pp == 1) {
        map <- attr(b1, "map");
        basesToKeep <- map[bases];
        isRepeat <- rep(TRUE, times=length(cells));
      } else {
        # Check for repeats
        isRepeat <- (b1 == b0);
      }
      skip <- !is.element(b1, basesToKeep);
      isRepeat[skip] <- FALSE;

      # Increment lengths of current repeats
      counts[isRepeat] <- counts[isRepeat] + as.integer(1);
      counts[!isRepeat] <- as.integer(1);

      # Check which ones are greater
      isGreater <- rep(FALSE, times=length(cells));
      keep <- is.finite(counts);
      isGreater[keep] <- (counts[keep] > maxRepeatLengthT[keep]);

      maxRepeatLengthT[isGreater] <- counts[isGreater];
      startPositionT[isGreater] <- pos - maxRepeatLengthT[isGreater] + as.integer(1);

      # Next position
      b0 <- b1;
      verbose && exit(verbose);
    } # for (pp ...)

    startPositionT[startPositionT == 0] <- as.integer(-1);
    maxRepeatLength[!nok] <- maxRepeatLengthT;
    startPosition[!nok] <- startPositionT;
    # Not needed anymore
    nok <- startPositionT <- maxRepeatLengthT <- NULL;
  }

  res <- cbind(startPosition, maxRepeatLength);
  # Not needed anymore
  startPosition <- maxRepeatLength <- NULL;

  verbose && cat(verbose, "Results:");
  verbose && str(verbose, res);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Caching result");
  saveCache(res, key=key, dirs=dirs);
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # getMaxLengthRepeats()



############################################################################
# HISTORY:
# 2008-12-30
# o Created.
############################################################################
