# (thetaA,thetaB) -> (theta, freqB)
#  theta = thetaA+thetaB
#  freqB = thetaB/theta
#
# (thetaA,thetaB) -> (theta, freqB)
#  thetaB = freqB*theta
#  thetaA = theta - thetaB = (1 - freqB)*theta
#
# (thetaA, thetaB) => (theta, freqB)
# ----------------------------------
#   (0,0) => (0, NA) => (NA,NA)***
#   (0,eps) => (eps, 1) => (0,eps)
#   (1,0) => (1,  0) => (1,0)
#   (0,1) => (1,  1) => (0,1)
#   (1,1) => (2,1/2) => (1,1)
#   (2,0) => (2,  0) => (2,0)
#   (1,1) => (2,1/2) => (1,1)
#   (0,2) => (2,  1) => (0,2)
#   (3,0) => (3,  0) => (3,0)
#   (2,1) => (3,1/3) => (2,1)
#   (1,2) => (3,2/3) => (1,2)
#   (0,3) => (3,  1) => (0,3)
#
#  (-1, 0) => (-1,  0) => (-1,  0)
#  ( 0,-1) => (-1,  1) => ( 0, -1)
#  (-2, 0) => (-2,  0) => (-2,  0)
#  (-1,-1) => (-2,1/2) => (-1, -1)
#  ( 0,-2) => (-2,  1) => ( 0, -2)
#  (-3, 0) => (-3,  0) => (-3,  0)
#  (-2,-1) => (-3,1/3) => (-2, -1)
#  (-1,-2) => (-3,2/3) => (-1, -2)
#  ( 0,-3) => (-3,  1) => ( 0, -3)
#
#  (-1, 1) => ( 0, NA)  => (NA, NA)***
#  (-1, 2) => ( 1,   2) => (-1,  2)
#  (-2, 1) => (-1,  -1) => (-2,  1)
#  (-1, 3) => ( 2, 3/2) => (-1,  3)
#  (-2, 4) => ( 2,   2) => (-2,  4)
#  (-3, 1) => (-2,-1/2) => (-3,  1)
#

setMethodS3("extractTotalAndFreqB", "CnChipEffectFile", function(this, units=NULL, ..., drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'units':
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    ugcMap <- NULL;
  } else if (isUnitGroupCellMap(units)) {
    ugcMap <- units;
    units <- unique(ugcMap[,"unit"]);
    nbrOfUnits <- length(units);
    # Not needed anymore
    units <- NULL;
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
    nbrOfUnits <- length(units);
    ugcMap <- NULL;
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Extracting (total, freqB)");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify possible groups
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (this$combineAlleles && this$mergeStrands) {
    # theta == (theta)
    groups <- 1;
  } else if (this$combineAlleles && !this$mergeStrands) {
    # theta == (theta+, theta-)
    groups <- c(1,3);
  } else if (!this$combineAlleles && this$mergeStrands) {
    # theta == (thetaA, thetaB)
    groups <- c(1,2);
  } else if (!this$combineAlleles && !this$mergeStrands) {
    # theta == (thetaA+, thetaB+, thetaA-, thetaB-)
    groups <- c(1,2,3,4);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the UGC map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(ugcMap)) {
    verbose && enter(verbose, "Getting (unit, group, cell) map");
    ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }

  if (!is.null(groups)) {
    idxs <- which(ugcMap$group %in% groups);
    ugcMap <- ugcMap[idxs,,drop=FALSE];
  } else {
    groups <- sort(unique(ugcMap$group));
  }

  verbose && cat(verbose, "Using (unit,group,cell) map:");
  verbose && str(verbose, ugcMap);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  theta <- extractTheta(this, units=ugcMap, groups=groups, ...,
                                                 verbose=less(verbose, 5));
  nbrOfUnits <- nrow(theta);

  # Calculating total chip effect
  thetaTotal <- rowSums(theta, na.rm=TRUE);

  # Calculating Allele B frequency
  if (this$combineAlleles) {
    freqB <- rep(NA, nbrOfUnits);
    # Not needed anymore
    theta <- NULL;
  } else {
    if (ncol(theta) == 2) {
      thetaB <- theta[,2];
    } else if (ncol(theta) == 4) {
      thetaB <- rowSums(theta[,c(2,4)], na.rm=TRUE);
    }
    # Not needed anymore
    theta <- NULL;
    freqB <- thetaB/thetaTotal;
    # Not needed anymore
    thetaB <- NULL;
  }

  data <- matrix(c(thetaTotal, freqB), nrow=nbrOfUnits, ncol=2);
  colnames(data) <- c("total", "freqB");

  # Drop singleton dimensions?
  if (drop) {
    data <- drop(data);
  }

  verbose && cat(verbose, "Results:");
  verbose && str(verbose, data);

  verbose && exit(verbose);

  data;
})



setMethodS3("extractTotalAndFracB", "SnpChipEffectFile", function(this, units=NULL, ..., drop=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'units':
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    ugcMap <- NULL;
  } else if (isUnitGroupCellMap(units)) {
    ugcMap <- units;
    units <- unique(ugcMap[,"unit"]);
    nbrOfUnits <- length(units);
    # Not needed anymore
    units <- NULL;
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
    nbrOfUnits <- length(units);
    ugcMap <- NULL;
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Extracting (total, freqB)");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify possible groups
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (this$mergeStrands) {
    # theta == (thetaA, thetaB)
    groups <- c(1,2);
  } else {
    # theta == (thetaA+, thetaB+, thetaA-, thetaB-)
    groups <- c(1,2,3,4);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the UGC map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(ugcMap)) {
    verbose && enter(verbose, "Getting (unit, group, cell) map");
    ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }

  if (!is.null(groups)) {
    idxs <- which(ugcMap$group %in% groups);
    ugcMap <- ugcMap[idxs,,drop=FALSE];
  } else {
    groups <- sort(unique(ugcMap$group));
  }

  verbose && cat(verbose, "Using (unit,group,cell) map:");
  verbose && str(verbose, ugcMap);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  theta <- extractTheta(this, units=ugcMap, groups=groups, ...,
                                                 verbose=less(verbose, 5));
  nbrOfUnits <- nrow(theta);

  # Calculating total chip effect
  thetaTotal <- rowSums(theta, na.rm=TRUE);

  # Calculating Allele B frequencies
  if (ncol(theta) == 2) {
    thetaB <- theta[,2];
  } else if (ncol(theta) == 4) {
    thetaB <- rowSums(theta[,c(2,4)], na.rm=TRUE);
  }
  # Not needed anymore
  theta <- NULL;
  freqB <- thetaB/thetaTotal;
  # Not needed anymore
  thetaB <- NULL;

  data <- matrix(c(thetaTotal, freqB), nrow=nbrOfUnits, ncol=2);
  colnames(data) <- c("total", "freqB");

  # Drop singleton dimensions?
  if (drop) {
    data <- drop(data);
  }

  verbose && cat(verbose, "Results:");
  verbose && str(verbose, data);

  verbose && exit(verbose);

  data;
})


setMethodS3("extractTotalAndFreqB", "default", function(this, ...) {
  extractTotalAndFracB(this, ...);
})



############################################################################
# HISTORY:
# 2009-01-10
# o Added extractTotalAndFracB().
# 2008-07-16
# o Added argument 'drop=FALSE' to all extractTotalAndFreqB().
# 2008-05-11
# o The (thetaA,thetaB) -> (theta, freqB) is bijective *except* when
#   theta = thetaA+thetaB = 0.  With the assumption that thetaA,thetaB > 0
#   then it is truly bijective.
# 2008-05-10
# o Now extractTotalAndFreqB() takes and UGC map via argument 'units'.
# 2008-05-09
# o Created.
############################################################################
