setMethodS3("getFitUnitGroupFunction", "SmoothMultiarrayModel", abstract=TRUE, protected=TRUE);


###########################################################################/**
# @set "class=SmoothMultiarrayModel"
# @RdocMethod fit
#
# @title "Fits the model for one chromosome across samples"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{data}{A @data.frame with columns \code{M} (log-ratio) and
#      \code{x} (locus position).
#   }
#   \item{chromosome}{An @integer specifying the index of the chromosome to
#      be fitted.}
#   \item{...}{Additional arguments passed down to the internal fit function.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fit", "SmoothMultiarrayModel", function(this, chromosome, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'chromosome':
  chromosome <- Arguments$getIndex(chromosome);
  knownChromosomes <- getChromosomes(this);
  if (!chromosome %in% knownChromosomes) {
    throw("Argument 'chromosome' contains an unknown chromosome: ", chromosome);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  fitOneChromosome(this, chromosome=chromosome, ..., verbose=less(verbose, 5));

  fit;
}, private=TRUE) # fit()





setMethodS3("getPositionChipTypeUnit", "CopyNumberSegmentationModel", function(this, chromosome, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'chromosome':
  chromosome <- Arguments$getIndex(chromosome);
  knownChromosomes <- getChromosomes(this);
  if (!chromosome %in% knownChromosomes) {
    throw("Argument 'chromosome' contains an unknown chromosome: ", chromosome);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (position, chipType, unit) map for this chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting (position, chipType, unit) map");

  # Get the UnitNameFile:s
  unfList <- getListOfUnitNamesFiles(this, verbose=less(verbose, 10));

  # Get the genome information files
  ugpList <- lapply(unfList, FUN=getAromaUgpFile, verbose=less(verbose, 10));
  verbose && print(verbose, ugpList);

  # Get the units on the chromosome of interest
  unitsList <- lapply(ugpList, FUN=function(ugp) {
    getUnitsOnChromosome(ugp, chromosome=chromosome, ...);
  });
  verbose && str(verbose, unitsList);
  # Not needed anymore
  ugpList <- NULL;

  # Gets (position, chipType) for these units
  posList <- vector("list", length(unitsList));
  names(posList) <- names(unitsList);
  chipTypeList <- vector("list", length(unitsList));
  names(chipTypeList) <- names(unitsList);
  for (kk in seq_along(posList)) {
    ugp <- ugpList[[kk]];
    units <- unitsList[[kk]];
    pos <- getPositions(ugp, units=units);

    # Keep only units with a position
    keep <- which(is.finite(pos));
    nbrOfUnitsBefore <- length(pos);
    nbrOfUnits <- length(keep);
    nbrOfUnitsExcl <- nbrOfUnitsBefore - nbrOfUnits;
    if (nbrOfUnitsExcl > 0) {
      pos <- pos[keep];
      units <- units[keep];
      verbose && cat(verbose, "Excluded ", nbrOfUnitsExcl, " (out of", nbrOfUnitsBefore, ") units because there is no position information available for those.");
    }
    unitsList[[kk]] <- units;
    posList[[kk]] <- pos;
    chipTypeList[[kk]] <- rep(kk, length(units));
    # Not needed anymore
    ugp <- units <- keep <- NULL;
  }
  # Not needed anymore
  ugpList <- NULL;

  verbose && str(verbose, unitsList);
  verbose && str(verbose, posList);
  verbose && str(verbose, chipTypeList);

  # Unlist and order (units, position, chipType) by position
  pos <- unlist(posList, use.names=FALSE);
  # Not needed anymore
  posList <- NULL;
  o <- order(pos);
  pos <- pos[o];

  chipType <- unlist(chipTypeList, use.names=FALSE);
  # Not needed anymore
  chipTypeList <- NULL;
  chipType <- chipType[o];

  # Convert chipType into a factor
  chipTypes <- sapply(unfList, FUN=getChipType);
  attr(chipType, "levels") <- chipTypes;
  class(chipType) <- "factor";
  # Not needed anymore
  unfList <- NULL;

  units <- unlist(unitsList, use.names=FALSE);
  # Not needed anymore
  unitsList <- NULL;
  units <- units[o];
  # Not needed anymore
  o <- NULL;

  pcu <- data.frame(position=pos, chipType=chipType, unit=units);
  # Not needed anymore
  units <- pos <- chipType <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && cat(verbose, "(position, chipType, unit) map:");
  verbose && str(verbose, pcu);

  verbose && exit(verbose);

  pcu;
}, protected=TRUE)



setMethodS3("createOutputTuple", "SmoothMultiarrayModel", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  outTuple <- this$.outTuple;
  if (!force && !is.null(outTuple))
    return(outTuple);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create chip-effect sets
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating output chip-effect tuple");
  inTuple <- getSetTuple(this);
  nbrOfChipTypes <- nbrOfChipTypes(inTuple);
  verbose && cat(verbose, "Number of chip types: ", nbrOfChipTypes);
  inList <- getSets(inTuple);
  outList <- vector("list", nbrOfChipTypes);
  names(outList) <- names(inList);
  parentPath <- getParentPath(this);
  for (kk in seq_along(inList)) {
    chipType <- names(inList)[kk];
    verbose && enter(verbose, "Chip type #", kk, "'(", chipType, ")' of ", nbrOfChipTypes);
    inSet <- inList[[chipType]];
    if (length(inSet) == 0)
      throw("Cannot create output data set. The input data set is empty.");

    chipType <- getChipType(inSet, fullname=FALSE);
    path <- Arguments$getWritablePath(file.path(parentPath, chipType));
    verbose && enter(verbose, "Creating output data set using input data set as a template (by copying)");
    verbose && cat(verbose, "Path: ", path);
    nbrOfArrays <- length(inSet);
    for (jj in seq_len(nbrOfArrays)) {
      inFile <- inSet[[jj]];
      outFile <- createFrom(inFile, filename=getFilename(inFile), path=path, methods=c("create", "copy"), clear=TRUE, verbose=less(verbose, 10));
    }

    # Speed this up by inheriting parameters from input data set
    args <- list(path);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # BEGIN: AFFX methods
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Ad hoc for ChipEffectSet classes. /HB 2007-09-25
    if (inherits(inSet, "SnpChipEffectSet"))
      args$mergeStrands <- inSet$mergeStrands;
    if (inherits(inSet, "CnChipEffectSet"))
      args$combineAlleles <- inSet$combineAlleles;
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # END: AFFX methods
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    verbose && str(verbose, args);
    args$verbose <- less(verbose, 30);
    staticFcn <- inSet$fromFiles;
    outSet <- do.call(staticFcn, args=args);
    verbose && print(verbose, outSet);
    verbose && exit(verbose);

    outList[[chipType]] <- outSet;
    verbose && exit(verbose);
  } # for (kk in ...)
  verbose && exit(verbose);

  outTuple <- newInstance(inTuple, outList);

  # Store in cache
  this$.outTuple <- outTuple;

  outTuple;
}, protected=TRUE)


setMethodS3("fitOneChromosome", "SmoothMultiarrayModel", function(this, chromosome, ..., vebose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fit one chromosome");

  verbose && cat(verbose, "Chromosome: ", chromosome);

  smoothFitFcn <- getFitUnitGroupFunction(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # BEGIN: AFFX methods
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cesList <- getSets(this);
  verbose && cat(verbose, "List of input data sets:");
  verbose && print(verbose, cesList);
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # END: AFFX methods
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Getting output data set (create if missing)
  outTuple <- getOutputTuple(this, verbose=less(verbose, 10));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extracting data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extract data across arrays and chip types");
  inData <- getPcuTheta(this, chromosome=chromosome, verbose=less(verbose, 10));
  colnames(inData$theta) <- getNames(cesList[[1]]); # AD HOC /2007-09-26
  verbose && cat(verbose, "Positions:");
  verbose && str(verbose, inData$pcu[,"position"]);
  verbose && cat(verbose, "thetas:");
  verbose && str(verbose, inData$theta);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Transforming data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Shift?
  shift <- this$.shift;
  verbose && cat(verbose, "shift:");
  verbose && str(verbose, shift);
  inData$theta <- inData$theta + shift;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit along chromosome
  verbose && enter(verbose, "Fit smoothWRMA()");
  bandwidth <- getBandwidth(this);
  verbose && printf(verbose, "Bandwidth: %.2fkb\n", bandwidth/1e3);
  sd <- bandwidth;

  # Prior weights?
  if (this$.weights == "1/s2") {
    # Calculate prior weights as the inverse variance of log ratios
    Y <- inData$theta;
    # AD HOC /HB 2007-09-26
    if (chromosome == 23) {
      # Keep only diploid samples
      n23 <- as.integer(sapply(cesList[[1]], FUN=getAttribute, "n23"));
      names(n23) <- getNames(cesList[[1]]);
      isDiploid <- (n23[colnames(Y)] == 2);
      Y <- Y[,isDiploid,drop=FALSE];
    }
    YR <- rowMedians(Y, na.rm=TRUE);
    M <- log2(Y/YR);
    s <- rowMads(M, na.rm=TRUE);
    w <- 1/(s^2);
    # Not needed anymore
    Y <- YR <- M <- s <- NULL;
    gc <- gc();
  } else {
    w <- NULL;
  }

  verbose && cat(verbose, "Prior weights:");
  verbose && str(verbose, w);
  fit <- smoothFitFcn(Y=inData$theta, x=inData$pcu[,"position"], w=w,
                          sd=sd, progress=TRUE, verbose=less(verbose, 10));
  inData$theta <- fit$theta;
  inData$phi <- fit$phi;
  verbose && cat(verbose, "Results:");
  verbose && str(verbose, inData);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Storing estimates");

  outList <- getSets(outTuple);
  verbose && cat(verbose, "List of output data sets:");
  verbose && print(verbose, outList);

  for (kk in seq_len(nbrOfChipTypes(this))) {
    verbose && enter(verbose, "Chip type #", kk, " of ", nbrOfChipTypes(this));

    outSet <- outList[[kk]];
    verbose && cat(verbose, "Output data set:");
    verbose && str(verbose, outSet);

    # Extract the data for this chip type
    idxs <- which(as.integer(inData$pcu[,"chipType"]) == kk);
    theta <- inData$theta[idxs,,drop=FALSE];
    phi <- inData$phi[idxs];
    units <- inData$pcu[idxs,"unit"];
    # Not needed anymore
    idxs <- NULL;

    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);


    map <- outData <- NULL;
    for (aa in seq_len(ncol(theta))) {
      verbose && enter(verbose, "Array #", aa, " of ", ncol(theta));

      outFile <- outSet[[aa]];
      verbose && cat(verbose, "Output data file:");
      verbose && str(verbose, outFile);

      if (is.null(map)) {
        # TODO: Create a (unit,cell) map
        verbose && enter(verbose, "Retrieving (unit,cell) map for all arrays");
        map <- getUnitGroupCellMap(outFile, units=units, verbose=less(verbose,2));
        verbose && str(verbose, map);
        # Not needed anymore
        # Not needed anymore
        units <- NULL;
        verbose && exit(verbose);
      }

      if (is.null(outData)) {
        # Allocate 'outData' object for writing
        outData <- data.frame(cell=map[,"cell"], theta=rep(0, nrow(map)));
      }

      outData[,"theta"] <- theta[,aa,drop=TRUE];
      verbose && cat(verbose, "(cell, theta):");
      verbose && str(verbose, outData);

      updateDataFlat(outFile, data=outData, verbose=less(verbose));

      # Not needed anymore
      outFile <- NULL;

      verbose && exit(verbose);
    } # for (aa in ...)

    # Clean up
    # Not needed anymore
    map <- outData <- theta <- phi <- NULL;
    gc <- gc();
    verbose && print(verbose, gc);

    # Not needed anymore
    outSet <- NULL;

    verbose && exit(verbose);
  } # for (kk in ...)
  verbose && exit(verbose);

  # Clean up
  # Not needed anymore
  outList <- inData <- NULL;
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  invisible(outTuple);
}, protected=TRUE)


############################################################################
# HISTORY:
# 2009-01-26
# o Updated getPositionChipTypeUnit() of SmoothMultiarrayModel to utilize
#   the UnitNamesFile Interface instead of assuming an AffymetrixCdfFile.
#   This requires aroma.core v1.0.1.
# 2007-09-26
# o Renamed to SmoothMultiarrayModel (from SrmaModel).
# o Added support for prior weights as the inverse of the log-ratio
#   variances.
# 2007-09-25
# o Renamed to SrmaModel (from CsrmaModel).
# o For now fit() calls fitOneChromosome(), but all this should go in a
#   different model class.
# o Added fitOneChromosome().
# o Added createOuputTuple().
# 2007-09-20
# o Created from GladModel.fitOne.R.
############################################################################
