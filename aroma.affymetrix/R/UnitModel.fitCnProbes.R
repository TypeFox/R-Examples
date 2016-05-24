# SPEEDUP TRICK: Update CN units (which are all single probe units) by copying
setMethodS3("fitCnProbes", "UnitModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Constants
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  .Machine$float.eps <- sqrt(.Machine$double.eps);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Estimating single-probe CN units");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting chip-effect set");
  # Get chip-effect set
  ces <- getChipEffectSet(this, verbose=verbose);
  verbose && exit(verbose);

  ds <- getDataSet(this);

  verbose && enter(verbose, "Identifying CN units");
  cdf <- getCdf(ds);
  units <- which(getUnitTypes(cdf) == 5);
  verbose && str(verbose, units);
  verbose && exit(verbose);

  # Nothing to do?
  if (length(units) == 0) {
    verbose && cat(verbose, "Nothing to do: This chip type does not have copy-number units");
    verbose && exit(verbose);
    return(invisible(units));
  }

  verbose && enter(verbose, "Identifying subset to fit");
  units2 <- findUnitsTodo(ces, verbose=verbose);
  units <- intersect(units, units2);
  # Not needed anymore
  units2 <- NULL;
  # Nothing to do?
  if (length(units) == 0) {
    verbose && cat(verbose, "Nothing to do: All copy-number units are fitted");
    verbose && exit(verbose);
    return(invisible(units));
  }
  verbose && str(verbose, units);
  verbose && exit(verbose);

  params <- getParameters(this);
  shift <- params$shift;

  verbose && enter(verbose, "Identifying cell indices for CN units");
  cells <- getCellIndices(cdf, units=units, useNames=FALSE);
  cells <- sapply(cells, FUN=unlist, use.names=FALSE);
  verbose && exit(verbose);

  verbose && enter(verbose, "Keeping only single-cell units");
  unitCellCounts <- sapply(cells, FUN=length);
  keep <- which(unitCellCounts == 1);
  # Not needed anymore
  unitCellCounts <- NULL;
  units <- units[keep];
  verbose && cat(verbose, "Single-cell units:");
  verbose && str(verbose, units);

  ## Nothing todo?
  if (length(units) == 0L) {
    verbose && cat(verbose, "Nothing to do: All copy-number units are fitted")
    verbose && exit(verbose)
    return(invisible(units))
  }

  cells <- cells[keep];
  cells <- unlist(cells, use.names=FALSE);
  verbose && cat(verbose, "Cell indices:");
  verbose && str(verbose, cells);
  # Sanity check
  if (length(cells) != length(units)) {
    throw("Detected units with more than one cell: ", getChipType(cdf));
  }
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying cell indices for estimates");
  cdfM <- getCdf(ces);
  cellsM <- getCellIndices(cdfM, units=units, useNames=FALSE, unlist=TRUE);
  verbose && str(verbose, cellsM);
  # Sanity check
  if (length(cellsM) != length(units))
    throw("Detected units with more than one cell: ", getChipType(cdfM));
  if (length(cellsM) != length(cells)) {
    throw("Input 'cells' and output 'cellsM' are of different lengths: ",
                                  length(cellsM), " != ", length(cells));
  }
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # "Fitting"
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Fitting ", length(ds), " arrays");
  res <- listenv()

  for (kk in seq_along(ds)) {
    df <- ds[[kk]];
    cef <- ces[[kk]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, getName(df), length(ds)));

    res[[kk]] %<=% {
      verbose && enter(verbose, sprintf("fitCnProbes(): Array '%s'", getName(df)));

      verbose && enter(verbose, "Reading signals");
      y <- extractMatrix(df, cells=cells, drop=TRUE);
      stopifnot(length(y) == length(cells));
      verbose && str(verbose, y);
      verbose && exit(verbose);

      # Shift?
      if (shift != 0) {
        verbose && enter(verbose, "Shifting signals");
        y <- y + shift;
        verbose && str(verbose, y);
        verbose && exit(verbose);
      }

      verbose && enter(verbose, "Transforming signals to estimates");
      sdTheta <- .Machine$float.eps;  # Smallest float > 0.
      data <- data.frame(cell=cellsM, theta=y, sdTheta=sdTheta, outliers=FALSE);
      # Not needed anymore
      y <- NULL;
      verbose && str(verbose, data);
      verbose && exit(verbose);

      verbose && enter(verbose, "Writing estimates");
      updateDataFlat(cef, data=data, verbose=verbose);
      # Not needed anymore
      data <- NULL;
      verbose && exit(verbose);

      verbose && exit(verbose);
    } # %<=%

    verbose && exit(verbose);
  } # for (kk ...)
  gc <- gc();
  verbose && print(verbose, gc);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(units);
}, private=TRUE) # fitCnProbes()


############################################################################
# HISTORY:
# 2009-05-20
# o Updated fitCnProbes() of UnitModel to identify single-cell CN units,
#   and ignore multi-cell CN units, which will be process like the other
#   units.  By not assuming single-cell CN units, this methods should also
#   apply to other CDFs, e.g. the new Cytogenetics_Array.
# 2008-12-01
# o Now fitCnProbes() of UnitModel only fits non-fitted CN units.
# 2008-09-05
# o Created.
############################################################################
