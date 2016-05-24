###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod isResequenceChip
#
# @title "Static method to check if a CDF is for a resequencing (CustomSeq) chip"
#
# \description{
#   @get "title".  Note, this method is not bullet proof.  Several
#   resequencing CDF does not carry that information.  For such, we add
#   tests based on their chip type, as we become aware of them.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE if the chip type refers to a resequence array,
#   otherwise @FALSE.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("isResequenceChip", "AffymetrixCdfFile", function(this, ...) {
  chipType <- getChipType(this);

  # First some hardwired return values
  if (regexpr("^Mitochip_2.*$", chipType) != -1) {
    return(TRUE);
  }
  if (regexpr("^MitoP-.*$", chipType) != -1) {
    return(TRUE);
  }

  # Then, check for resequencing units
  types <- getUnitTypes(this, ...);
  hasReseqUnits <- any(types == 3);

  hasReseqUnits;
}, private=TRUE)





###########################################################################/**
# @RdocMethod readUnitsByQuartets
#
# @title "Gets the cell quartets for each base position"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{units}{Subset of units to be queried. If @NULL, all units are used.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @vector of @factors.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("readUnitsByQuartets", "AffymetrixCdfFile", function(this, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (is.null(units)) {
    units <- seq_len(nbrOfUnits(this));
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(this));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Read (pbase, tbase, cell index)
  verbose && enter(verbose, "Reading (pbase, tbase, cell index)");
  cdfUnits <- readUnits(this, units=units, readIndices=TRUE, readBases=TRUE, readXY=FALSE, readExpos=FALSE, readType=FALSE, readDirection=FALSE);
  verbose && cat(verbose, "Number of units read: ", length(cdfUnits));
  verbose && exit(verbose);

  verbose && enter(verbose, "Restructuring (and validating assumptions about) fields 'pbase', 'tbase', and 'indices'");
  verbose && cat(verbose, "Number of units: ", length(cdfUnits));
  # Restructure and validate
  for (uu in seq_along(cdfUnits)) {
    verbose && enter(verbose, sprintf("Unit #%d ('%s') of %d", uu, names(cdfUnits)[uu], length(cdfUnits)));

    cdfUnit <- cdfUnits[[uu]];
    cdfGroups <- cdfUnit$groups;
    verbose && cat(verbose, "Number of groups: ", length(cdfGroups));

    for (gg in seq_along(cdfGroups)) {
      verbose && enter(verbose, sprintf("Group #%d ('%s') of %d", gg, names(cdfGroups)[gg], length(cdfGroups)));

      cdfGroup <- cdfGroups[[gg]];
      # Sanity check of assumption of ordering of cells
      pbase <- cdfGroup$pbase;
      pbase <- matrix(pbase, nrow=4, byrow=FALSE);
      pbase <- t(pbase);
      pbase <- unique(pbase);
      if (nrow(pbase) != 1) {
        throw("Assumption exception: The probe bases ('pbase') are not ordered consistently for this unit: ", units[uu]);
      }
      pbase <- as.vector(pbase);
      pbase <- toupper(pbase);

      tbase <- cdfGroup$tbase;
      tbase <- matrix(tbase, nrow=4, byrow=FALSE);
      tbase1 <- tbase[1,,drop=TRUE];
      # Sanity check
      if (!all(apply(tbase, MARGIN=1, FUN=identical, tbase1))) {
        throw("Assumption exception: The target bases ('tbase') are not ordered consitently for this unit: ", units[uu]);
      }
      tbase <- tbase1;
      tbase <- toupper(tbase);

      cells <- cdfGroup$indices;
      cells <- matrix(cells, nrow=4, byrow=FALSE);
      rownames(cells) <- pbase;
#      colnames(cells) <- tbase;

      cdfGroup <- list(indices=cells, tbase=tbase);

      cdfGroups[[gg]] <- cdfGroup;
      # Not needed anymore
      pbase <- cells <- cdfGroup <- NULL;
      verbose && exit(verbose);
    } # for (gg ...)

    cdfUnit$groups <- cdfGroups;
    cdfUnits[[uu]] <- cdfUnit;
    # Not needed anymore
    cdfUnit <- NULL;

    verbose && exit(verbose);
  } # for (uu ...)
  verbose && exit(verbose);

  cdfUnits;
}, private=TRUE)



setMethodS3("getCellQuartets", "AffymetrixCdfFile", function(this, units=NULL, mergeGroups=TRUE, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  cdfGetFields <- affxparser::cdfGetFields


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (is.null(units)) {
  } else {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(this)));
  }

  # Argument 'mergeGroups':
  mergeGroups <- Arguments$getLogical(mergeGroups);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting cell-index matrices");

  key <- list(method="getCellQuartets", class=class(this)[1],
              units=units, mergeGroups=mergeGroups);
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="getCellQuartets", units=units, mergeGroups=mergeGroups);
  }
  dirs <- c("aroma.affymetrix", getChipType(this));
  if (!force) {
    cells <- loadCache(key=key, dirs=dirs);
    if (!is.null(cells)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(cells);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read raw CDF data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdfUnits <- readUnitsByQuartets(this, units=units, verbose=verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract cell quartets and annotate with (pbase, tbase)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pbase <- rownames(cdfUnits[[1]]$groups[[1]]$indices);
  # Sanity check
  if (is.null(pbase)) {
    throw("No resequencing cell indices available.");
  }
  if (!all(is.element(c("A", "C", "G", "T"), pbase))) {
    throw("No resequencing cell indices available.");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract cell quartets by CDF groups
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Restructuring into a matrices");

  # Extract field 'tbase'
  tbase <- .applyCdfGroups(cdfUnits, cdfGetFields, "tbase");
  tbase <- .applyCdfGroups(tbase, function(groups) {
    lapply(groups, FUN=.subset2, 1);
  });
  tbase <- lapply(tbase, FUN=.subset2, "groups");

  # Extract field 'cells'
  cells <- .applyCdfGroups(cdfUnits, cdfGetFields, "indices");
  cells <- .applyCdfGroups(cells, function(groups) {
    lapply(groups, FUN=.subset2, 1);
  });
  cells <- lapply(cells, FUN=.subset2, "groups");
  # Not needed anymore
  cdfUnits <- NULL;  # Not needed anymore

  # Attach 'tbase' as column names to 'cells'
  for (uu in seq_along(cells)) {
    cellsUU <- cells[[uu]];
    tbaseUU <- tbase[[uu]];
    for (gg in seq_along(cellsUU)) {
      cellsGG <- cellsUU[[gg]];
      tbaseGG <- tbaseUU[[gg]];
      colnames(cellsGG) <- tbaseGG;
      cellsUU[[gg]] <- cellsGG;
    } # for (gg ...)
    cells[[uu]] <- cellsUU;
    # Not needed anymore
    cellsUU <- tbaseUU <- NULL;
  } # for (uu ...)
  # Not needed anymore
  tbase <- cellsGG <- tbaseGG <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Merging groups
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (mergeGroups) {
    for (uu in seq_along(cells)) {
      cellsUU <- cells[[uu]];
      cells[[uu]] <- cellsUU;
      cellsUUMerged <- NULL;
      for (gg in seq_along(cellsUU)) {
        cellsGG <- cellsUU[[gg]];
        cellsUUMerged <- cbind(cellsUUMerged, cellsGG);
      } # for (gg ...)
      cells[[uu]] <- list(all=cellsUUMerged);
      # Not needed anymore
      cellsUU <- cellsUUMerged <- NULL;
    } # for (uu ...)
  }

  # Save to cache?
  if (cache) {
    saveCache(cells, key=key, dirs=dirs);
  }

  verbose && exit(verbose);

  cells;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2008-09-29
# o ROBUSTNESS: Now readUnitsByQuartets() for AffymetrixCdfFile translates
#   lower-case pbase and tbase letter to upper case.  It also asserts that
#   order of the 'tbase' is consistent with the expectations.
#   TODO: This method should be renamed to indicate that it is intended
#   for the resequencing arrays.
# 2008-08-29
# o Added argument 'mergeGroups' to getCellQuartets().
# 2008-08-18
# BUG FIX: Used non-existing 'cdf' instead of 'this'.
# 2008-08-10
# o Created.
############################################################################
