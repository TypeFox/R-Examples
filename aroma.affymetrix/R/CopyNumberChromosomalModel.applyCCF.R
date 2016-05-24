setMethodS3("applyAnyOrder", "CopyNumberChromosomalModel", function(this, chipTypes=NULL, arrays=NULL, chromosomes=NULL, FUN, order=c("cca", "cac"), ..., verbose=FALSE) {
  allChipTypes <- getChipTypes(this);
  allChromosomes <- getChromosomes(this);
  allArrays <- getNames(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipTypes':
  if (is.null(chipTypes)) {
    chipTypes <- allChipTypes;
  } else if (is.numeric(chipTypes)) {
    chipTypes <- Arguments$getIndices(chipTypes, max=length(allChipTypes));
    chipTypes <- allChipTypes[chipTypes];
  } else {
    if (!all(chipTypes %in% allChipTypes)) {
      throw("Argument 'chipTypes' contains an unknown value: ",
                                          paste(chipTypes, collapse=", "));
    }
  }

  # Argument 'chromosomes':
  if (is.null(chromosomes)) {
    chromosomes <- allChromosomes;
  } else if (is.numeric(chromosomes)) {
    chromosomes <- Arguments$getIndices(chromosomes, max=length(allChromosomes));
    chromosomes <- allChromosomes[chromosomes];
  } else {
    if (!all(chromosomes %in% allChromosomes)) {
      throw("Argument 'chromosomes' contains an unknown value: ",
                                        paste(chromosomes, collapse=", "));
    }
  }

  # Argument 'arrays':
  if (is.null(arrays)) {
    arrayNames <- allArrays;
    arrays <- seq_along(allArrays);
  } else if (is.character(arrays)) {
    if (!all(arrays %in% allArrays)) {
      throw("Argument 'arrays' contains an unknown value: ",
                                        paste(arrays, collapse=", "));
    }
    arrayNames <- arrays;
    arrays <- match(arrayNames, allArrays);
  } else {
    arrays <- Arguments$getIndices(arrays, max=length(allArrays));
    arrayNames <- allArrays[arrays];
  }

  # Argument 'FUN':
  if (!is.function(FUN)) {
    throw("Argument 'FUN' is not a function: ", mode(FUN));
  }

  # Argument 'order':
  order <- match.arg(order);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "applyAnyOrder()");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the test and the reference set tuples
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Get set tuples");
  cesTuple <- getSetTuple(this);
  verbose && exit(verbose);
  verbose && enter(verbose, "Get reference set tuples");
  refTuple <- getReferenceSetTuple(this); # Different from getRefSetTuple().
  verbose && exit(verbose);

  verbose && enter(verbose, "Extract chip types");
  # Extract chip types of interest
  if (identical(chipTypes, getChipTypes(cesTuple))) {
    cesTuple <- newInstance(cesTuple, getSets(cesTuple)[chipTypes]);
  }
  if (identical(chipTypes, getChipTypes(refTuple))) {
    refTuple <- newInstance(refTuple, getSets(refTuple)[chipTypes]);
  }
  verbose && exit(verbose);

  # Extract the arrays of interest
  verbose && enter(verbose, "Extract arrays of interest");
  cesTuple <- extract(cesTuple, arrays=arrays, onDuplicates="error");
  refTuple <- extract(refTuple, arrays=arrays, onDuplicates="error");
  verbose && exit(verbose);

  nbrOfArrays <- length(cesTuple);

  # Sanity check
  stopifnot(length(refTuple) == nbrOfArrays);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the (UGC, chromosome, position) maps
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (1) Identify one ChipEffectFile from each chip type
  arrayTable <- getTableOfArrays(cesTuple);
  pair <- apply(arrayTable, MARGIN=2, FUN=function(idx) {
    which.min(is.finite(idx));
  });
  if (any(is.na(pair))) {
    throw("Internal error: Could not identify a pair of ChipEffectFile:s");
  }
  # Not needed anymore
  arrayTable <- NULL;

  cesSets <- getSets(cesTuple);
  nbrOfChipTypes <- nbrOfChipTypes(cesTuple);
  ceFiles <- vector("list", nbrOfChipTypes);
  names(ceFiles) <- chipTypes;
  for (kk in seq_len(nbrOfChipTypes)) {
    ces <- cesSets[[kk]];
    ceFiles[[kk]] <- ces[[pair[kk]]];
    # Not needed anymore
    ces <- NULL;
  }

  # (2) For each of them, extract one map per chromsome
  nbrOfChromosomes <- length(chromosomes);
  maps <- vector("list", nbrOfChipTypes);
  names(maps) <- chipTypes;
  for (kk in seq_len(nbrOfChipTypes)) {
    chipType <- chipTypes[[kk]];
    verbose && enter(verbose, sprintf("Chip type #%d ('%s') of %d",
                                          kk, chipType, nbrOfChipTypes));

    mapsT <- vector("list", nbrOfChromosomes);
    names(mapsT) <- sprintf("Chr%02d", as.integer(chromosomes));
    maps[[chipType]] <- mapsT;
    # Not needed anymore
    mapsT <- NULL;

    ce <- ceFiles[[kk]];
    for (ll in seq_len(nbrOfChromosomes)) {
      chromosome <- chromosomes[ll];
      verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d",
                                      ll, chromosome, nbrOfChromosomes));

      map <- getUnitGroupCellChromosomePositionMap(ce, chromosome=chromosome, verbose=less(verbose, 5));
      maps[[kk]][[ll]] <- map;
      # Not needed anymore
      map <- NULL;

      verbose && exit(verbose);
    }
    # Not needed anymore
    ce <- NULL;
    verbose && exit(verbose);
  }
  # Not needed anymore
  ceFiles <- NULL;
  verbose && cat(verbose, "(UGC, chromosome, position) maps:");
  verbose && str(verbose, maps);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculating raw CNs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cesSets <- getSets(cesTuple);
  refSets <- getSets(refTuple);

  if (order == "cca") {
    values <- applyCCF0(this, cesSets=cesSets, refSets=refSets,
                      maps=maps, FUN=FUN, ..., verbose=less(verbose, 5));
  } else if (order == "cac") {
    values <- applyCFC0(this, cesSets=cesSets, refSets=refSets,
                      maps=maps, FUN=FUN, ..., verbose=less(verbose, 5));
  } else {
    throw("Unknown order: ", paste(order, collapse=", "));
  }

  verbose && exit(verbose);

  invisible(values);
}, private=TRUE);



setMethodS3("applyCCF0", "CopyNumberChromosomalModel", function(this, cesSets, refSets, maps, FUN, ..., verbose=FALSE) {
  nbrOfChipTypes <- length(maps);
  chipTypes <- names(maps);

  nbrOfChromosomes <- length(maps[[1]]);
  chromosomes <- names(maps[[1]]);

  nbrOfArrays <- length(cesSets[[1]]);
  arrays <- seq_len(nbrOfArrays);


  valuesCT <- vector("list", nbrOfChipTypes);
  names(valuesCT) <- chipTypes;
  for (kk in seq_len(nbrOfChipTypes)) {
    chipType <- chipTypes[[kk]];
    verbose && enter(verbose, sprintf("Chip type #%d ('%s') of %d",
                                          kk, chipType, nbrOfChipTypes));

    ces <- cesSets[[kk]];
    res <- refSets[[kk]];

    valuesC <- vector("list", nbrOfChromosomes);
    names(valuesC) <- chromosomes;
    valuesCT[[kk]] <- valuesC;

    for (ll in seq_len(nbrOfChromosomes)) {
      chromosome <- chromosomes[ll];
      verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d",
                                      ll, chromosome, nbrOfChromosomes));

      map <- maps[[kk]][[ll]];
      verbose && cat(verbose, "(UGC, chromosome, position) map:");
      verbose && str(verbose, map);

      valuesA <- vector("list", nbrOfArrays);
      names(valuesA) <- getNames(ces);
      valuesCT[[kk]][[ll]] <- valuesA;

      for (mm in seq_len(nbrOfArrays)) {
        ce <- ces[[mm]];
        re <- res[[mm]];
        verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                          mm, getName(ce), nbrOfArrays));

        # Work around... /HB 2008-03-11
        class(map) <- c("UnitGroupCellMap", class(map));

        value <- FUN(chipType=kk, chromosome=ll, array=mm, map=map,
                                     ce=ce, re=re, ..., verbose=verbose);

        valuesCT[[kk]][[ll]][[mm]] <- value;

        verbose && exit(verbose);
      }

      # Not needed anymore
      map <- NULL;
      verbose && exit(verbose);
    }

    verbose && exit(verbose);
  }

  invisible(valuesCT);
}, private=TRUE)


setMethodS3("applyCFC0", "CopyNumberChromosomalModel", function(this, cesSets, refSets, maps, FUN, ..., verbose=FALSE) {
  nbrOfChipTypes <- length(maps);
  chipTypes <- names(maps);

  nbrOfChromosomes <- length(maps[[1]]);
  chromosomes <- names(maps[[1]]);

  nbrOfArrays <- length(cesSets[[1]]);
  arrays <- seq_len(nbrOfArrays);
  arrayNames <- getNames(cesSets[[1]]);

  valuesC <- vector("list", nbrOfChromosomes);
  names(valuesC) <- chromosomes;

  for (ll in seq_len(nbrOfChromosomes)) {
    chromosome <- chromosomes[ll];
    verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d",
                                    ll, chromosome, nbrOfChromosomes));

    valuesA <- vector("list", nbrOfArrays);
    names(valuesA) <- arrayNames;
    valuesC[[ll]] <- valuesA;

    for (mm in seq_len(nbrOfArrays)) {
      ces <- cesSets[[1]];
      ce <- ces[[mm]];
      verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                        mm, getName(ce), nbrOfArrays));
      valuesCT <- vector("list", nbrOfChipTypes);
      names(valuesCT) <- chipTypes;
      valuesC[[ll]][[mm]] <- valuesCT;

      for (kk in seq_len(nbrOfChipTypes)) {
        chipType <- chipTypes[[kk]];
        verbose && enter(verbose, sprintf("Chip type #%d ('%s') of %d",
                                            kk, chipType, nbrOfChipTypes));

        ces <- cesSets[[kk]];
        res <- refSets[[kk]];

        map <- maps[[kk]][[ll]];
        verbose && cat(verbose, "(UGC, chromosome, position) map:");
        verbose && str(verbose, map);

        # Work around... /HB 2008-03-11
        class(map) <- c("UnitGroupCellMap", class(map));

        ce <- ces[[mm]];
        re <- res[[mm]];

        value <- FUN(chipType=kk, chromosome=ll, array=mm, map=map,
                                     ce=ce, re=re, ..., verbose=verbose);

        valuesC[[ll]][[mm]][[kk]] <- value;

        verbose && exit(verbose);
      }

      # Not needed anymore
      map <- NULL;
      verbose && exit(verbose);
    }

    verbose && exit(verbose);
  }

  invisible(valuesC);
}, private=TRUE)



setMethodS3("calcRawCnStats", "default", function(M, ...) {
  # Number of data points
  nAll <- length(M);

  # Ignore NAs
  ok <- !is.na(M);
  M <- M[ok];
  # Not needed anymore
  ok <- NULL;

  # Number of non-missing data points
  n <- length(M);

  # Mean log-ratio
  mu <- median(M);

  # Standard deviation
  dM <- diff(M);
  sigma <- mad(dM)/sqrt(2);

  list(mu=mu, sigma=sigma, n=n, nAll=nAll);
}, private=TRUE) # calcRawCnStats()


setMethodS3("calcRawCnStats", "CopyNumberChromosomalModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  myFUN <- function(chipType=NULL, chromosome=NULL, array=NULL, map=NULL, ce=NULL, re=NULL, ..., cumulative=FALSE, verbose=FALSE) {
    theta0  <- extractMatrix(ce, units=map, verbose=less(verbose, 25));
    thetaR0 <- extractMatrix(re, units=map, verbose=less(verbose, 25));

    # Cumulative and same chip type?
    if (cumulative) {
      if (!exists("lastArray", mode="numeric")) {
        lastArray <- NA;
      }
      if (identical(array, lastArray)) {
        theta <<- c(theta, theta0);
        thetaR <<- c(thetaR, thetaR0);
        pos <<- c(pos, map[,"physicalPosition"]);
        # Reorder
        o <- order(pos);
        pos <<- pos[o];
        theta <<- theta[o];
        thetaR <<- thetaR[o];
        # Not needed anymore
        o <- NULL;
      } else {
        theta <<- theta0;
        thetaR <<- thetaR0;
        pos <<- map[,"physicalPosition"];
        lastArray <<- array;
      }
    } else {
      theta <- theta0;
      thetaR <- thetaR0;
    }

    # The log-ratios are already ordered along the chromosome
    calcRawCnStats(log2(theta/thetaR));
  } # myFUN()


  verbose && enter(verbose, "calcRawCnStats()");

  theta <- thetaR <- pos <- NULL;  # Dummies for '<<-' above.
  values <- applyAnyOrder(this, ..., order="cac",
                          FUN=myFUN, cumulative=TRUE, verbose=verbose);

  # Extract the cumulative statistics (stored in the last chip type)
  nbrOfChipTypes <- nbrOfChipTypes(this);
  values <- lapply(values, FUN=function(files) {
    lapply(files, FUN=.subset2, nbrOfChipTypes);
  });
#  verbose && str(verbose, values);

  mu <- t(sapply(values, FUN=function(files) {
    sapply(files, FUN=.subset2, "mu");
  }))
  verbose && str(verbose, mu);

  sigma <- t(sapply(values, FUN=function(files) {
    sapply(files, FUN=.subset2, "sigma");
  }))
  verbose && str(verbose, sigma);

  verbose && exit(verbose);

  list(mu=mu, sigma=sigma);
}, private=TRUE);


##############################################################################
# HISTORY:
# 2008-03-11
# o Created.
##############################################################################
