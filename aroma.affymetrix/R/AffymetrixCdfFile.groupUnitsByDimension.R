###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod groupUnitsByDimension
#
# @title "Groups units by dimensions"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{An optional @integer @vector specifying the units to be
#     queried.}
#   \item{...}{Not used.}
#   \item{sort}{If @TRUE, the sets are ordered by number of groups per
#     units and then by the number of cells per group.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a named @list.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("groupUnitsByDimension", "AffymetrixCdfFile", function(this, units=NULL, ..., sort=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(this));
    allUnits <- units;
  } else {
    allUnits <- seq_len(nbrOfUnits(this));
  }

  # Argument 'sort':
  sort <- Arguments$getLogical(sort);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Grouping units by dimensions");

  verbose && cat(verbose, "Units: ");
  verbose && str(verbose, allUnits);

  # Identify SNPs with equal number of cells
  sizes <- nbrOfCellsPerUnitGroup(this, units=units, verbose=less(verbose,1));

  # Identify sets of units of equal dimensions
  nbrOfGroupsPerUnit <- sapply(sizes, FUN=length);
  t <- table(nbrOfGroupsPerUnit);
  verbose && print(verbose, t);
  ## nbrOfGroupsPerUnit
  ##      2      4
  ## 906874   2748

  uSizes <- as.integer(names(t));

  # Sort?
  if (sort) {
    uSizes <- sort(uSizes);
  }

  res <- list();

  for (uu in seq_along(uSizes)) {
    sizeUU <- uSizes[uu];
    verbose && enter(verbose, sprintf("Unit size #%d (nbrOfGroups=%d) of %d", uu, sizeUU, length(uSizes)));
    verbose && cat(verbose, "Number of groups per unit: ", sizeUU);

    # Pull out all units 'size' number of groups
    idxsUU <- which(nbrOfGroupsPerUnit == sizeUU);
    unitsUU <- allUnits[idxsUU];
    verbose && cat(verbose, "Units: ");
    verbose && str(verbose, unitsUU);

    # Group all those units into those with equal number cells per group
    sizesUU <- sizes[idxsUU];
    sizesUU <- unlist(sizesUU, use.names=FALSE);
    sizesUU <- matrix(sizesUU, ncol=sizeUU, byrow=TRUE);
    uDimsUU <- unique(sizesUU);

    # Sort?
    if (sort) {
      o <- order(uDimsUU[,1]);
      uDimsUU <- uDimsUU[o,,drop=FALSE];
    }

    verbose && cat(verbose, "Unique dimensions: ");
    verbose && print(verbose, uDimsUU);

    resUU <- list();
    resUU$nbrOfGroups <- sizeUU;
    resUU$units <- unitsUU;
    resUU$sets <- list();

    # For each group of equal dimension
    for (kk in seq_len(nrow(uDimsUU))) {
      dimKK <- uDimsUU[kk, ,drop=TRUE];
      verbose && enter(verbose, sprintf("Dimension #%d of %d", kk, nrow(uDimsUU)));
      verbose && cat(verbose, "Number of cells per group: ", paste(dimKK, collapse=", "));

      keep <- rep(TRUE, times=length(unitsUU));
      for (cc in seq_len(ncol(uDimsUU))) {
        keep <- keep & (dimKK[cc] == sizesUU[,cc, drop=TRUE]);
      }
      idxsKK <- which(keep);
      unitsKK <- unitsUU[idxsKK];
#      verbose && cat(verbose, "Units: ");
#      verbose && str(verbose, unitsKK);

      resKK <- list();
      resKK$nbrOfGroups <- length(dimKK);
      resKK$nbrOfUnits <- length(unitsKK);
      resKK$nbrOfCellsPerGroup <- dimKK;
      resKK$units <- unitsKK;
      verbose && str(verbose, resKK);

      resUU$sets[[kk]] <- resKK;
      # Not needed anymore
      idxsKK <- unitsKK <- resKK <- NULL;

      verbose && exit(verbose);
    } # for (kk ...)

    res[[uu]] <- resUU;
    # Not needed anymore
    keep <- idxsUU <- unitsUU <- resUU <- NULL;

    verbose && exit(verbose);
  } # for (uu ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Asserting correctness");
  # (a) Check units
  unitsT <- lapply(res, FUN=function(x) x$units);
  unitsT <- unlist(unitsT, use.names=FALSE);
  unitsT <- sort(unitsT);
  stopifnot(identical(unitsT, allUnits));

  # (b) Check units in subelements
  setsT <- lapply(res, FUN=function(x) x$sets);
  unitsT <- lapply(setsT, FUN=function(x) {
    lapply(x, FUN=function(y) y$units)
  });
  unitsT <- unlist(unitsT, use.names=FALSE);
  unitsT2 <- sort(unitsT);
  stopifnot(identical(unitsT2, allUnits));

  # Flatten sets
  sets <- lapply(res, FUN=function(x) { x$sets });
  # Merge sets
  sets <- Reduce(append, sets);

  # Tabulate dimensions per set
  dims <- sapply(sets, FUN=function(set) c(set$nbrOfGroups, set$nbrOfUnits));
  dims <- t(dims);
  colnames(dims) <- c("nbrOfGroups", "nbrOfUnits");
  dims <- as.data.frame(dims);
  nbrOfCellsPerGroup <- lapply(sets, FUN=function(x) { x$nbrOfCellsPerGroup });
  dims$nbrOfCellsPerGroup <- nbrOfCellsPerGroup;
  dims <- dims[,c(1L,3L,2L)];

  # Sanity check
  stopifnot(sum(dims[,"nbrOfUnits"]) == length(allUnits));
  verbose && exit(verbose);

  res <- list(nestedSets=res, sets=sets, setDimensions=dims);

  verbose && exit(verbose);

  res;
}, protected=TRUE) # groupUnitsByDimension()


############################################################################
# HISTORY:
# 2011-11-18
# o Added column 'nbrOfCellsPerGroup' containing a list of ragged vectors
#   to element 'setDimensions' returned by groupUnitsByDimension().
# 2011-11-10
# o ROBUSTNESS: Added a sanity check to groupUnitsByDimension() for
#   AffymetrixCdfFile.
# 2010-04-21
# o Created.
############################################################################
