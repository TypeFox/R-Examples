###########################################################################/**
# @RdocClass ProbeLevelTransform3
#
# @title "The ProbeLevelTransform3 class"
#
# \description{
#  @classhierarchy
#
#  This abstract class is specialized from @see "ProbeLevelTransform" and
#  provides methods to identify subsets and types of probes that are used
#  for fitting and/or updating the signals.
# }
#
# @synopsis
#
# \arguments{
#   \item{dataSet}{A @see "AffymetrixCelSet".}
#   \item{...}{Arguments passed to the constructor of
#     @see "ProbeLevelTransform".}
#   \item{unitsToFit}{The units from which the normalization curve should
#     be estimated.  If @NULL, all are considered.}
#   \item{typesToFit}{Types of probes to be used when fitting the model.}
#   \item{unitsToUpdate}{The units to be updated.
#     If @NULL, all are considered.}
#   \item{typesToUpdate}{Types of probes to be updated.}
#   \item{shift}{An optional amount to shift data before fitting and updating.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("ProbeLevelTransform3", function(dataSet=NULL, ..., unitsToFit="-XY", typesToFit=typesToUpdate, unitsToUpdate=NULL, typesToUpdate="pm", shift=0) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  extraTags <- NULL;

  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "AffymetrixCelSet");

    # Argument 'typesToUpdate':
    if (!is.null(typesToUpdate)) {
      typesToUpdate <- match.arg(typesToUpdate, choices=c("pm", "mm", "pmmm"));
    }

    # Argument 'unitsToUpdate':
    if (is.null(unitsToUpdate)) {
    } else if (is.character(unitsToUpdate)) {
      if (unitsToUpdate %in% c("-X", "-Y", "-XY")) {
      } else {
        throw("Unknown value of argument 'unitsToUpdate': ", unitsToUpdate);
      }
      extraTags <- c(extraTags, unitsToUpdate=unitsToUpdate);
    } else {
      df <- getOneFile(dataSet);
      nbrOfCells <- nbrOfCells(df);
      unitsToUpdate <- Arguments$getIndices(unitsToUpdate, max=nbrOfCells);
      unitsToUpdate <- unique(unitsToUpdate);
      unitsToUpdate <- sort(unitsToUpdate);
    }

    # Argument 'typesToFit':
    if (!is.null(typesToFit)) {
      typesToFit <- match.arg(typesToFit, choices=c("pm", "mm", "pmmm"));
    }

    # Argument 'unitsToFit':
    if (is.null(unitsToFit)) {
    } else if (is.character(unitsToFit)) {
      if (unitsToFit %in% c("-X", "-Y", "-XY")) {
      } else {
        throw("Unknown value of argument 'unitsToFit': ", unitsToFit);
      }
      extraTags <- c(extraTags, unitsToFit=unitsToFit);
    } else {
      df <- getOneFile(dataSet);
      nbrOfCells <- nbrOfCells(df);
      unitsToFit <- Arguments$getIndices(unitsToFit, max=nbrOfCells);
      unitsToFit <- unique(unitsToFit);
      unitsToFit <- sort(unitsToFit);
    }
  }

  # Argument 'shift':
  shift <- Arguments$getDouble(shift, disallow=c("NA", "NaN", "Inf"));


  extend(ProbeLevelTransform(dataSet=dataSet, ...), "ProbeLevelTransform3",
    "cached:.cellsToUpdate" = NULL,
    "cached:.cellsToFit" = NULL,
    .typesToUpdate = typesToUpdate,
    .unitsToUpdate = unitsToUpdate,
    .typesToFit = typesToFit,
    .unitsToFit = unitsToFit,
    .extraTags = extraTags,
    shift = shift
  )
})


setMethodS3("getAsteriskTags", "ProbeLevelTransform3", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Extra tags?
  tags <- c(tags, this$.extraTags);

  # Add class-specific tags
  shift <- as.integer(round(this$shift));
  if (shift != 0) {
    tags <- c(tags, sprintf("%+d", shift));
  }

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)



setMethodS3("getUnitsTo", "ProbeLevelTransform3", function(this, what, ..., force=FALSE) {
  what <- capitalize(what);
  idxField <- sprintf(".unitsTo%s", what);
  units <- this[[idxField]];
  dataSet <- getInputDataSet(this);
  cdf <- getCdf(dataSet);
  units <- getSubsetOfUnits(cdf, units=units, ..., force=force);
  units;
}, protected=TRUE)


setMethodS3("getUnitsToFit", "ProbeLevelTransform3", function(this, ...) {
  getUnitsTo(this, what="fit", ...);
}, protected=TRUE)


setMethodS3("getUnitsToUpdate", "ProbeLevelTransform3", function(this, ...) {
  getUnitsTo(this, what="update", ...);
}, protected=TRUE)




setMethodS3("getCellsToInternal", "ProbeLevelTransform3", function(this, units, stratifyBy, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (is.null(units)) {
  } else if (is.numeric(units)) {
  } else {
    throw("Internal error. Unknown data type for argument 'units': ", mode(units));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting cell indices for subset of units");

  dataSet <- getInputDataSet(this);
  cdf <- getCdf(dataSet);
  chipType <- getChipType(cdf);

  verbose && cat(verbose, "Dataset class:");
  verbose && cat(verbose, class(dataSet)[1]);
  verbose && cat(verbose, "Chip type: ", chipType);
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Treat ChipEffectSet specially
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (inherits(dataSet, "ChipEffectSet")) {
    verbose && enter(verbose, "Treating ", class(dataSet)[1], " specially");
    verbose && cat(verbose, "Requested 'units': ");
    verbose && str(verbose, units);
    stratifyBy <- NULL;
    extractKey <- list(class(dataSet)[1], params=getParameters(dataSet));
    verbose && exit(verbose);
  } else {
    extractKey <- NULL;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="getCellsToInternal", class=class(this)[1],
    chipType=chipType, units=units, stratifyBy=stratifyBy,
    extractKey=extractKey
  );
  dirs <- c("aroma.affymetrix", chipType);
  res <- loadCache(key=key, dirs=dirs);
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "Found cached results!");
    verbose && exit(verbose);
    return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the cell indices for the requested units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the subset of cells according to the CDF
  cells <- getSubsetOfCellIndices(cdf, units=units, stratifyBy=stratifyBy,
                                                 verbose=less(verbose, 5));
  verbose && cat(verbose, "Cells: ");
  verbose && str(verbose, cells);


  verbose && cat(verbose, "Identified cells:");
  verbose && str(verbose, cells);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Treat ChipEffectSet specially
  # If a ChipEffectSet, then so called 'restructuring' might occur,
  # which will exclude additional cells.  For this reason, we have
  # check what possible cells the ChipEffectSet uses and take the
  # intersection between those and the above found cells.
  # /HB 2008-07-25.  This was borrowed from the same idea done in
  # ScaleNormalization on 2007-04-11.  Ideally we will transfer the
  # whole restructuring to the CDF with a more intuitive API.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (inherits(dataSet, "ChipEffectSet")) {
    df <- getOneFile(dataSet);

    # Cannot use 'unlist=TRUE' next, because restructuring might occur.
    possibleCells <- getCellIndices(df, verbose=less(verbose));
    possibleCells <- unlist(possibleCells, use.names=FALSE);
    possibleCells <- sort(possibleCells);
    verbose && cat(verbose, "Cells used by this ", class(df)[1], ":");
    verbose && str(verbose, possibleCells);

    # Take the intersection
    if (is.null(cells)) {
      cells <- possibleCells;
    } else {
      cells <- intersect(cells, possibleCells);
    }
    # Not needed anymore
    possibleCells <- NULL;

    verbose && cat(verbose, "Cells: ");
    verbose && str(verbose, cells);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Cache results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache)
    saveCache(cells, key=key, dirs=dirs);

  cells;
}, private=TRUE)




setMethodS3("getCellsTo", "ProbeLevelTransform3", function(this, what, ..., force=FALSE) {
  what <- capitalize(what);
  idxField <- sprintf(".cellsTo%s", what);

  # Check for cached results
  cells <- this[[idxField]];
  if (!force && !is.null(cells))
    return(cells);

  # Get the units
  units <- getUnitsTo(this, what=what, ..., force=force);

  # Construct subset
  typesField <- sprintf(".typesTo%s", what);
  stratifyBy <- this[[typesField]];

  cells <- getCellsToInternal(this, units=units, stratifyBy=stratifyBy, ..., force=force);

  # Cache results
  this[[idxField]] <- cells;

  cells;
}, protected=TRUE)


setMethodS3("getCellsToFit", "ProbeLevelTransform3", function(this, ...) {
  getCellsTo(this, what="fit", ...);
}, protected=TRUE)


setMethodS3("getCellsToUpdate", "ProbeLevelTransform3", function(this, ...) {
  getCellsTo(this, what="update", ...);
}, protected=TRUE)





setMethodS3("getParameters", "ProbeLevelTransform3", function(this, expand=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Getting algorithm parameters");

  # Get parameters from super class
  params <- NextMethod("getParameters", expand=expand);

  # Get local parameters
  params2 <- list(
    unitsToFit = this$.unitsToFit,
    typesToFit = this$.typesToFit,
    unitsToUpdate = this$.unitsToUpdate,
    typesToUpdate = this$.typesToUpdate,
    shift = this$shift
  );

  # Expand to unit and cell indices?
  if (expand) {
    verbose && enter(verbose, "Expanding unit indices");
    params2$unitsToFit <- getUnitsToFit(this, verbose=verbose);
    params2$unitsToUpdate <- getUnitsToUpdate(this, verbose=verbose);
    verbose && exit(verbose);
    verbose && enter(verbose, "Expanding cell indices");
    params2$cellsToFit <- getCellsToFit(this, verbose=verbose);
    params2$cellsToUpdate <- getCellsToUpdate(this, verbose=verbose);
    verbose && exit(verbose);
  }

  params <- c(params, params2);

  verbose && exit(verbose);

  params;
}, protected=TRUE)



setMethodS3("writeSignals", "ProbeLevelTransform3", function(this, pathname, cells=NULL, ..., templateFile, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'templateFile':
  templateFile <- Arguments$getInstanceOf(templateFile, "AffymetrixCelFile");

  # Argument 'cells':
  if (is.null(cells)) {
  } else {
    # Validated below...
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Writing probe signals");
  verbose && cat(verbose, "Output pathname: ", pathname);
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, list(...));

  # Create CEL file to store results, if missing
  verbose && enter(verbose, "Creating CEL file for results, if missing");

  # Write to a temporary file
  pathnameT <- pushTemporaryFile(pathname, isFile=FALSE, verbose=verbose);

  createFrom(templateFile, filename=pathnameT, path=NULL, verbose=verbose);
  verbose && exit(verbose);

  verbose && enter(verbose, "Storing normalized signals");
  .updateCel(pathnameT, indices=cells, ...);
  verbose && exit(verbose);

  # Rename temporary file
  popTemporaryFile(pathnameT, verbose=verbose);

  verbose && exit(verbose);

  invisible(pathname);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2008-07-25
# o Added protected writeSignals(), which write signals to a CEL file,
#   which is created from a template CEL file, if missing.  This method
#   cen then later be designed to be atomic.
# o Now getSubsetTo() also handles ChipEffectSet, which might use
#   so called "restructors".
# 2008-07-20
# o Extracted from BaseCountNormalization.R.
############################################################################
