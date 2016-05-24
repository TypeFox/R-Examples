###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod getAlleleCellPairs
#
# @title "Gets the cell indices of allele pairs"
#
# \description{
#   @get "title"
#   in units of type "genotyping" with 2 or 4 unit groups.
# }
#
# @synopsis
#
# \arguments{
#  \item{units}{A @integer @vector of units to query.
#    If @NULL, all units are considered.}
#  \item{stratifyBy}{A @character string specifying what type of probes
#    to return.}
#  \item{...}{Not used.}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a Nx2 @integer @matrix of cell indices, where each row is
#   a (PMA, PMB) probe pair.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAlleleCellPairs", "AffymetrixCdfFile", function(this, units=NULL, stratifyBy=c("pm", "pmmm", "mm"), force=FALSE, ..., verbose=FALSE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser")
  cdfMergeAlleles <- affxparser::cdfMergeAlleles


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(this));
  }

  # Argument 'stratifyBy':
  stratifyBy <- match.arg(stratifyBy);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Identifying the probe pairs");
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this);
  key <- list(method="getAlleleCellPairs", class=class(this)[1],
                    chipType=chipType, units=units, stratifyBy=stratifyBy);
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="getAlleleCellPairs", chipType=chipType, units=units, stratifyBy=stratifyBy);
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify genotype units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying genotyping units (with 2 or 4 groups)");
  # Use only units that are SNPs...
  verbose && enter(verbose, "Reading unit types");
  unitTypes <- getUnitTypes(this, units=units, verbose=less(verbose, 1));
  verbose && exit(verbose);

  keep <- which(unitTypes == 2);
  if (is.null(units)) {
    units <- keep;
  } else {
    units <- units[keep];
  }
  # Not needed anymore
  unitTypes <- keep <- NULL;

  # ...and with either 2 or 4 groups
  verbose && enter(verbose, "Reading number of groups per SNP unit");
  unitSizes <- nbrOfGroupsPerUnit(this, units=units);
  verbose && cat(verbose, "Detected unit sizes:");
  verbose && print(verbose, table(unitSizes));
  verbose && exit(verbose);

  keep <- which(is.element(unitSizes, c(2,4)));
  units <- units[keep];
  # Not needed anymore
  unitSizes <- keep <- NULL;

  nbrOfUnits <- length(units);
  verbose && cat(verbose, "Number of SNP units to query: ", nbrOfUnits);

  verbose && exit(verbose);

  if (nbrOfUnits == 0) {
    verbose && exit(verbose);
    return(NULL);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading cell indices");
  verbose && cat(verbose, "Stratify by: ", stratifyBy);
  cells <- getCellIndices(this, units=units, stratifyBy=stratifyBy,
                                          useNames=FALSE, verbose=verbose);
  # Not needed anymore
  units <- NULL;

  verbose && enter(verbose, "Merging groups by allele pair");
  verbose && printf(verbose, "Units left: ");
  for (uu in seq_along(cells)) {
    if (uu %% 5000 == 0)
      verbose && writeRaw(verbose, length(cells)-uu, ", ");
    unit <- cells[[uu]];
    groups <- unit$groups;
    groups <- cdfMergeAlleles(groups);
    unit$groups <- groups;
    cells[[uu]] <- unit;
    if (uu %% 100000 == 0) {
      gc <- gc();
      verbose && print(verbose, gc);
    }
  }
  verbose && writeRaw(verbose, "0.\n");
  # Not needed anymore
  unit <- groups <- uu <- NULL;
  verbose && exit(verbose);

  cells <- unlist(cells, use.names=FALSE);
  cells <- matrix(cells, nrow=2);
  rownames(cells) <- c("A", "B");

  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Caching result");
  saveCache(cells, key=key, dirs=dirs);
  verbose && exit(verbose);

  verbose && exit(verbose);

  cells;
}, protected=TRUE) # getAlleleCellPairs()



############################################################################
# HISTORY:
# 2008-09-02
# o Added getAlleleCellPairs() for AffymetrixCdfFile.
# o Created.
############################################################################

