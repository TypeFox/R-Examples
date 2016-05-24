###########################################################################/**
# @RdocClass FragmentEquivalentClassNormalization
#
# @title "The FragmentEquivalentClassNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization method that corrects for systematic
#  effects between loci of different equivalent classes of pairs of
#  sequences that are recognized by the restriction enzymes that cut the
#  DNA studies.
# }
#
# @synopsis
#
# \arguments{
#   \item{dataSet}{A @see "CnChipEffectSet".}
#   \item{...}{Additional arguments passed to the constructor of
#     @see "ChipEffectTransform".}
#   \item{targetAvgs}{An optional list of @functions.
#     For each enzyme there is one target averages to which all arrays
#     should be normalized to.}
#   \item{subsetToFit}{The units from which the normalization curve should
#     be estimated.  If @NULL, all are considered.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Requirements}{
#   This class requires an UFC (Unit Fragment Class) annotation file.
# }
#
# \section{Acknowledgments}{
#   The idea of normalization signals stratified on enzyme
#   recognition sequences is credited to Jim Veitch and
#   Ben Bolstad at Affymetrix Inc. (2008) who have designed
#   a similar method for copy number estimation in the
#   Affymetrix' Genotype Console v2.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("FragmentEquivalentClassNormalization", function(dataSet=NULL, ..., targetAvgs=NULL, subsetToFit="-XY") {
  extraTags <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "CnChipEffectSet");

    if (dataSet$combineAlleles != TRUE) {
      throw("Currently only total copy-number chip effects can be normalized, i.e. 'combineAlleles' must be TRUE");
    }

#    if (dataSet$mergeStrands != TRUE) {
#      throw("Currently only non-strands specific copy-number chip effects can be normalized, i.e. 'mergeStrands' must be TRUE");
#    }
  }

  # Argument 'targetAvgs':
  if (!is.null(targetAvgs)) {
    if (!is.list(targetAvgs)) {
      throw("Argument 'targetAvgs' is not a list: ", class(targetAvgs)[1]);
    }

    # Validate each element
    for (kk in seq_along(targetAvgs)) {
    }
  }

  # Argument 'subsetToFit':
  if (is.null(subsetToFit)) {
  } else if (is.character(subsetToFit)) {
    if (subsetToFit %in% c("-X", "-Y", "-XY")) {
    } else {
      throw("Unknown value of argument 'subsetToFit': ", subsetToFit);
    }
    extraTags <- c(extraTags, subsetToFit=subsetToFit);
  } else {
    cdf <- getCdf(dataSet);
    subsetToFit <- Arguments$getIndices(subsetToFit, max=nbrOfUnits(cdf));
    subsetToFit <- unique(subsetToFit);
    subsetToFit <- sort(subsetToFit);
  }


  extend(ChipEffectTransform(dataSet, ...), "FragmentEquivalentClassNormalization",
    .subsetToFit = subsetToFit,
    .extraTags = extraTags,
    .targetAvgs = targetAvgs,
    "cached:.units" = NULL,
    "cached:.unitSets" = NULL,
    "cached:.ufc" = NULL
  )
})


setMethodS3("getAsteriskTags", "FragmentEquivalentClassNormalization", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Extra tags?
  tags <- c(tags, this$.extraTags);

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)



setMethodS3("getParameters", "FragmentEquivalentClassNormalization", function(this, expand=TRUE, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters", expand=expand);

  # Get parameters of this class
  params <- c(params, list(
    subsetToFit = this$.subsetToFit,
    .targetAvgs = this$.targetAvgs
  ));

  # Expand?
  if (expand) {
    subsetToFit <- getSubsetToFit(this);
  }

  params;
}, protected=TRUE)


setMethodS3("getCdf", "FragmentEquivalentClassNormalization", function(this, ...) {
  inputDataSet <- getInputDataSet(this);
  getCdf(inputDataSet);
})


setMethodS3("getAromaUfcFile", "FragmentEquivalentClassNormalization", function(this, force=FALSE, ...) {
  ufc <- this$.ufc;

  if (force || is.null(ufc)) {
    cdf <- getCdf(this);
    chipType <- getChipType(cdf);
    nbrOfUnits <- nbrOfUnits(cdf);
    ufc <- AromaUfcFile$byChipType(chipType, nbrOfUnits=nbrOfUnits);
    this$.ufc <- ufc;
  }

  ufc;
})


setMethodS3("getOutputDataSet00", "FragmentEquivalentClassNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting input data set");
  ces <- getInputDataSet(this);
  verbose && cat(verbose, "Class: ", class(ces)[1]);
  verbose && exit(verbose);

  verbose && enter(verbose, "Getting output data set for ", class(this)[1]);

  args <- list("getOutputDataSet");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Inherit certain arguments from the input data set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (inherits(ces, "CnChipEffectSet"))
    args$combineAlleles <- ces$combineAlleles;
  if (inherits(ces, "SnpChipEffectSet"))
    args$mergeStrands <- ces$mergeStrands;

  verbose && cat(verbose, "Calling NextMethod() with arguments:");
  verbose && str(verbose, args);

  args$verbose <- less(verbose, 10);
  res <- do.call(NextMethod, args);

  # Carry over parameters too.  AD HOC for now. /HB 2007-01-07
  if (inherits(res, "SnpChipEffectSet")) {
    verbose && enter(verbose, "Carrying down parameters for ", class(res)[1]);

    res$mergeStrands <- ces$mergeStrands;
    if (inherits(res, "CnChipEffectSet")) {
      res$combineAlleles <- ces$combineAlleles;
    }
    verbose && exit(verbose);
  }

  # Let the set update itself
  update2(res, verbose=less(verbose, 1));

  verbose && exit(verbose);

  res;
}, protected=TRUE)

setMethodS3("getSubsetToFit", "FragmentEquivalentClassNormalization", function(this, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Cached?
  units <- this$.units;
  if (!is.null(units) && !force)
    return(units);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identifying all potential units, i.e. SNPs and CN probes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying units that are SNP and CN probes");

  # Get fragment class information
  ufc <- getAromaUfcFile(this);

  # Identify all SNP and CN units (==potential units to be fitted)
  verbose && enter(verbose, "Identifying SNPs and CN probes");
  cdf <- getCdf(this);
## OLD:  units <- indexOf(cdf, "^(SNP|CN)");
  types <- getUnitTypes(cdf);
  units <- which(types == 2 | types == 5);
  # Not needed anymore
  types <- NULL;
  verbose && str(verbose, units);
  verbose && exit(verbose);

  # Keep only those for which we have fragmenth equivalent class information
  # for at least one enzyme
  verbose && enter(verbose, "Reading fragment equivalent classes");
  ufe <- getOrderedFragmentPairs(ufc, verbose=less(verbose, 5));
  keep <- rep(FALSE, nrow(ufe));
  for (ee in seq_len(ncol(ufe))) {
    keep <- keep | is.finite(ufe[,ee]);
  }
  units <- units[keep];
  verbose && printf(verbose, "Number of SNP/CN units without fragment equivalent class annotation: %d out of %d (%.1f%%)\n", sum(!keep), length(keep), 100*sum(!keep)/length(keep));

  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset with a prespecified set of units?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subsetToFit <- this$.subsetToFit;
  if (is.character(subsetToFit)) {
    if (subsetToFit %in% c("-X", "-Y", "-XY")) {
      verbose && enter(verbose, "Identify subset of units from genome information");
      verbose && cat(verbose, "subsetToFit: ", subsetToFit);

      # Look up in cache
      subset <- this$.subsetToFitExpanded;
      if (is.null(subset)) {
        # Get the genome information (throws an exception if missing)
        gi <- getGenomeInformation(cdf);
        verbose && print(verbose, gi);

        # Identify units to be excluded
        if (subsetToFit == "-X") {
          subset <- getUnitsOnChromosomes(gi, 23, .checkArgs=FALSE);
        } else if (subsetToFit == "-Y") {
          subset <- getUnitsOnChromosomes(gi, 24, .checkArgs=FALSE);
        } else if (subsetToFit == "-XY") {
          subset <- getUnitsOnChromosomes(gi, 23:24, .checkArgs=FALSE);
        }

        verbose && cat(verbose, "Units to exclude: ");
        verbose && str(verbose, subset);

        # The units to keep
        subset <- setdiff(1:nbrOfUnits(cdf), subset);

        verbose && cat(verbose, "Units to include: ");
        verbose && str(verbose, subset);

        # Store
        this$.subsetToFitExpanded <- subset;
      }

      subsetToFit <- subset;
      # Not needed anymore
      subset <- NULL;

      verbose && exit(verbose);
    }
  }

  if (!is.null(subsetToFit)) {
    # A fraction subset?
    if (length(subsetToFit) == 1 && 0 < subsetToFit && subsetToFit < 1) {
      keep <- seq(from=1, to=length(units), length=subsetToFit*length(units));
    } else {
      keep <- which(units %in% subsetToFit);
    }

    units <- units[keep];
    # Not needed anymore
    keep <- NULL;
  }

  # Sort units
  units <- sort(units);

  # Assert correctness
  units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));

  # Cache
  this$.units <- units;

  # Clear dependent fields
  this$.targetAvgs <- NULL;
  this$.unitSets <- NULL;

  verbose && exit(verbose);

  units;
}, private=TRUE)



setMethodS3("getExclusiveUnitSubsets", "FragmentEquivalentClassNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  unitSets <- this$.unitSets;
  if (is.null(unitSets) || force) {
    verbose && enter(verbose, "Identifying sets of units to be averaged over");

    # Get fragment class information
    ufc <- getAromaUfcFile(this);
    verbose && print(verbose, ufc);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    # Get units to fit
    units <- getSubsetToFit(this, verbose=less(verbose, 5));
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);

    # Get fragmenth equivalent class information for these
    ufe <- getOrderedFragmentPairs(ufc, units=units, verbose=less(verbose, 5));
    # Not needed anymore
    ufc <- NULL; # Not needed anymore
    verbose && cat(verbose, "Fragment equivalent classes data:");
    verbose && str(verbose, ufe);
    verbose && cat(verbose, "Unique equivalent classes:");
    verbose && str(verbose, sort(unique(as.vector(ufe))));

    #
    verbose && enter(verbose, "Identifying sets");
    hasUfe <- is.finite(ufe);
    nbrOfEnzymes <- ncol(ufe);
    allEnzymes <- seq_len(nbrOfEnzymes);

    unitSets <- list();
    for (ee in allEnzymes) {
      verbose && enter(verbose, "Enzyme #", ee, " of ", nbrOfEnzymes);

      # a) Fit only to units with known length...
      ok <- hasUfe[,ee];

      # b) ...exclude multi-enzyme units
      for (ff in setdiff(allEnzymes, ee))
        ok <- ok & !hasUfe[,ff];

      # Sanity check
      if (sum(ok) == 0) {
        throw("Cannot fit target function to enzyme, because there are no (finite) data points that are unique to this enzyme: ", ee);
      }

      # Fit effects to each fragment equivalent class
      data <- ufe[ok,ee];
      # Not needed anymore
      ok <- NULL;

      classes <- sort(unique(data));
      nbrOfEqClasses <- length(classes);
      verbose && cat(verbose, "Number of equivalent classes (in data): ",
                                                         nbrOfEqClasses);

      sets <- vector("list", nbrOfEqClasses);
      names(sets) <- sprintf("0x%02x", classes);
      for (cc in seq_len(nbrOfEqClasses)) {
        verbose && enter(verbose, "Equivalent class #", cc,
                                  " ('", names(sets)[cc], "')");
        subset <- which(data == classes[cc]);
        sets[[cc]] <- list(
          name   = names(sets)[cc],
          class  = classes[cc],
          subset = subset
        );
        # Not needed anymore
        subset <- NULL;
        verbose && exit(verbose);
      }

      unitSets[[ee]] <- sets;
      # Not needed anymore
      sets <- NULL;

      verbose && exit(verbose);
    }

    verbose && cat(verbose, "Identified sets: ");
    verbose && str(verbose, unitSets);

    # Remove as many promises as possible
    # Not needed anymore
    nbrOfEnzymes <- allEnzymes <- ufe <- hasUfe <- NULL;
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);

    this$.unitSets <- unitSets;
  }

  unitSets;
}, private=TRUE) # getExclusiveUnitSubsets()




setMethodS3("calculateAverages", "FragmentEquivalentClassNormalization", function(this, y, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  if (!is.numeric(y)) {
    throw("Argument 'y' is not numeric: ", mode(y));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calculating averages");

  verbose && cat(verbose, "Signals:");
  verbose && str(verbose, y);

  # Get subsets of units to be fitted
  unitSubsets <- getExclusiveUnitSubsets(this, verbose=less(verbose, 5));

  # Calculate the average for each subset
  verbose && enter(verbose, "Calculating the average for each subset");
  nbrOfEnzymes <- length(unitSubsets);
  allEnzymes <- seq_len(nbrOfEnzymes);

  avgs <- vector("list", nbrOfEnzymes);
  names(avgs) <- names(unitSubsets);
  for (ee in allEnzymes) {
    verbose && enter(verbose, "Enzyme #", ee, " of ", nbrOfEnzymes);

    subsets <- unitSubsets[[ee]];
    nbrOfSubsets <- length(subsets);

    avgsEE <- vector("list", nbrOfSubsets);
    names(avgsEE) <- names(subsets);

    for (kk in seq_len(nbrOfSubsets)) {
#      verbose && enter(verbose, "Subset #", kk, " of ", nbrOfSubsets);
      subset <- subsets[[kk]];

      # Sanity check
      if (length(subset$subset) == 0) {
         # Should never happen
         throw("Cannot calculate average, because there are no units in this set: (enzyme,class)=(", ee, ",", subset$name, ")");
      }

      # Estimate mean and standard deviation.
      data <- y[subset$subset];
      data <- data[!is.na(data)];

      mu <- median(data, na.rm=FALSE);
      sigma <- mad(data, na.rm=FALSE);
      n <- length(data);

      # Store
      avgsEE[[kk]] <- list(mu=mu, sigma=sigma, n=n);

      # Not needed anymore
      data <- mu <- sigma <- n <- subset <- NULL;

#      verbose && exit(verbose);
    } # for (kk ...)

    avgs[[ee]] <- avgsEE;
    # Not needed anymore
    avgsEE <- NULL;

    verbose && exit(verbose);
  } # for (ee ...)

  # Remove as many promises as possible
  # Not needed anymore
  nbrOfEnzymes <- allEnzymes <- y <- NULL;
  verbose && exit(verbose);

  verbose && exit(verbose);

  avgs;
}, private=TRUE)


setMethodS3("getTargetAverages", "FragmentEquivalentClassNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  log2 <- TRUE;

  avgs <- this$.targetAvgs;
  if (is.null(avgs) || force) {
    verbose && enter(verbose, "Calculating target averages");

    # Get target set
    ces <- getInputDataSet(this);
    verbose && enter(verbose, "Get average signal across arrays");
    ceR <- getAverageFile(ces, force=force, verbose=less(verbose));
    # Not needed anymore
    ces <- NULL; # Not needed anymore
    verbose && exit(verbose);

    cef <- ceR;
    # Not needed anymore
    ceR <- NULL;

    # Read signals
    units <- getSubsetToFit(this, verbose=less(verbose, 5));
    y <- getDataFlat(cef, units=units, fields="theta",
                               verbose=less(verbose, 5))[,"theta"];
    # Not needed anymore
    cef <- NULL; # Not needed anymore

    if (log2) {
      suppressWarnings({
        y <- log2(y);
      });
    }

    verbose && cat(verbose, "Signals:");
    verbose && str(verbose, y);

    avgs <- calculateAverages(this, y=y, verbose=less(verbose,5));
    # Not needed anymore
    y <- NULL;

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);

    this$.targetAvgs <- avgs;
  }

  avgs;
}, private=TRUE)



setMethodS3("normalizeOneArrayVector", "FragmentEquivalentClassNormalization", function(this, y, ufe, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  if (!is.numeric(y)) {
    throw("Argument 'y' is not numeric: ", mode(y));
  }

  # Argument 'ufe':
  if (!is.matrix(ufe)) {
    throw("Argument 'ufe' is not a matrix: ", class(ufe)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  log2 <- TRUE;

  verbose && enter(verbose, "Normalizing one array vector");

  # Get units to fit
  units <- getSubsetToFit(this, verbose=less(verbose, 5));

  # Get target averages
  verbose && enter(verbose, "Calculating target averages");
  targetAvgs <- getTargetAverages(this, verbose=less(verbose, 5));
  targetAvgs <- lapply(targetAvgs, FUN=function(enzyme) {
    sapply(enzyme, FUN=.subset2, "mu");
  })
  verbose && exit(verbose);

  if (log2) {
    verbose && enter(verbose, "Taking the log2");
    suppressWarnings({
      y <- log2(y);
    });
    verbose && exit(verbose);
  }

  verbose && cat(verbose, "All signals:");
  verbose && str(verbose, y);

  # Calculating sample averages
  verbose && enter(verbose, "Calculating sample averages");
  avgs <- calculateAverages(this, y=y[units], verbose=less(verbose,5));
  verbose && str(verbose, avgs);
  avgs <- lapply(avgs, FUN=function(enzyme) {
    sapply(enzyme, FUN=.subset2, "mu");
  })
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate correction factors for each enzyme equivalent class
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating correction factors");

  nbrOfEnzymes <- length(targetAvgs);
  deltas <- vector("list", nbrOfEnzymes);
  for (ee in seq_len(nbrOfEnzymes))
    deltas[[ee]] <- targetAvgs[[ee]] - avgs[[ee]];
  verbose && print(verbose, deltas);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalizing all units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  naValue <- as.double(NA);
  deltas2 <- matrix(naValue, nrow=length(y), ncol=nbrOfEnzymes);
  for (ee in seq_len(nbrOfEnzymes)) {
    verbose && enter(verbose, "Enzyme #", ee, " of ", nbrOfEnzymes);

    # Get correction factors
    delta <- deltas[[ee]];

    # Identify the class IDs, e.g. 0x14
    classIds <- sapply(names(delta), FUN=function(x) eval(parse(text=x)));

    # Normalize on log2-scale?
    if (log2) {
      suppressWarnings({
        delta <- 2^delta;
      })
    }

    # Get fragment equivalent class for each unit
    class <- ufe[,ee,drop=TRUE];

    # Map to correction factor
    idxs <- match(class, classIds);
    # Not needed anymore
    class <- classIds <- NULL;

    # Get correction factor
    deltas2[,ee] <- delta[idxs];
    # Not needed anymore
    delta <- idxs <- NULL;

    verbose && exit(verbose);
  } # for (ee ...)

  # Calculate the number of correction factors per unit
  counts <- integer(nrow(deltas2));
  dy <- double(nrow(deltas2));
  for (ee in seq_len(nbrOfEnzymes)) {
    ok <- is.finite(deltas2[,ee]);
    counts <- counts + as.integer(ok);
    dy[ok] <- dy[ok] + deltas2[ok,ee];
    # Not needed anymore
    ok <- NULL;
  } # for (ee ...)
  # Not needed anymore
  deltas2 <- NULL;

  dy[counts == 0] <- NA;
  dy <- dy / counts;

  # Take anti-log2?
  if (log2) {
    verbose && enter(verbose, "Taking the log2 of the correction factors");
    suppressWarnings({
      dy <- log2(dy);
    });
    verbose && exit(verbose);
  }

  verbose && cat(verbose, "Summary of correction factors:");
  verbose && summary(verbose, dy);

  verbose && cat(verbose, "Counts:");
  verbose && print(verbose, table(counts));

  # Normalize
  ok <- which(is.finite(dy));
  y[ok] <- y[ok] + dy[ok];
  # Not needed anymore
  ok <- NULL;

  # Take anti-log2?
  if (log2) {
    verbose && enter(verbose, "Taking the anti-log2");
    suppressWarnings({
      y <- 2^y;
    });
    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  y;
}, private=TRUE)




###########################################################################/**
# @RdocMethod process
#
# @title "Normalizes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already normalized is re-normalized,
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "FragmentEquivalentClassNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Normalizing set for systematic effects between fragment equivalent classes");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already normalized");
    verbose && enter(verbose, "Getting output data set");
    outputSet <- getOutputDataSet(this, verbose=less(verbose, 10));
    verbose && exit(verbose);
    verbose && exit(verbose);
    return(invisible(outputSet));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ces <- getInputDataSet(this);

  # Get SNP units
  cdf <- getCdf(ces);
# OLD:  subsetToUpdate <- indexOf(cdf, "^(SNP|CN)");
  types <- getUnitTypes(cdf);
  subsetToUpdate <- which(types == 2 | types == 5);
  # Not needed anymore
  types <- NULL;

  verbose && enter(verbose, "Retrieving fragment class annotations");
  chipType <- getChipType(cdf);
  nbrOfUnits <- nbrOfUnits(cdf);
  ufc <- AromaUfcFile$byChipType(chipType, nbrOfUnits=nbrOfUnits);
  verbose && cat(verbose, "Number of enzymes: ", nbrOfEnzymes(ufc));
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying the subset used to fit normalization function(s)");
  # Get subset to fit
  subsetToFit <- getSubsetToFit(this, verbose=less(verbose));
  verbose && str(verbose, subsetToFit);
  verbose && exit(verbose);

  # Get (and create) the output path
  path <- getPath(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ufe <- map <- NULL;
  nbrOfArrays <- length(ces);
  for (kk in seq_len(nbrOfArrays)) {
    ce <- ces[[kk]];
    verbose && enter(verbose, sprintf("Array #%d of %d ('%s')",
                                            kk, nbrOfArrays, getName(ce)));

    filename <- getFilename(ce);
    pathname <- filePath(path, filename);
    if (isFile(pathname)) {
      verbose && cat(verbose, "Already normalized. Skipping.");
      ceN <- fromFile(ce, pathname);

      # Carry over parameters too.  AD HOC for now. /HB 2007-01-07
      if (inherits(ce, "SnpChipEffectFile")) {
        ceN$mergeStrands <- ce$mergeStrands;
        if (inherits(ce, "CnChipEffectFile")) {
          ceN$combineAlleles <- ce$combineAlleles;
        }
      }

      # CDF inheritance
      setCdf(ceN, cdf);

      verbose && exit(verbose);
      next;
    }

    if (is.null(map)) {
      # Only loaded if really needed.
      verbose && enter(verbose, "Retrieving unit-to-cell map for all arrays");
      map <- getUnitGroupCellMap(ce, units=subsetToUpdate, verbose=less(verbose));
      verbose && str(verbose, map);
      verbose && exit(verbose);
    }

    if (is.null(ufe)) {
      verbose && enter(verbose, "Reading fragment equivalent classes");
      ufc <- getAromaUfcFile(this);
      ufe <- getOrderedFragmentPairs(ufc, units=subsetToUpdate,
                                                verbose=less(verbose));
      # Not needed anymore
      ufc <- NULL;
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Reading signals");
    data <- getDataFlat(ce, units=map, fields="theta", verbose=less(verbose));
    y <- data[,"theta"];
    verbose && str(verbose, y);
    verbose && exit(verbose);

    verbose && enter(verbose, "Normalizing signals");
    yN <- normalizeOneArrayVector(this, y=y, ufe=ufe, ..., verbose=less(verbose));
    # Not needed anymore
    y <- NULL;
    verbose && str(verbose, yN);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    # Create CEL file to store results, if missing
    verbose && enter(verbose, "Creating CEL file for results, if missing");
    ceN <- createFrom(ce, filename=pathname, path=NULL, verbose=less(verbose));
    verbose && exit(verbose);

    # Carry over parameters too.  AD HOC for now. /HB 2007-01-07
    if (inherits(ce, "SnpChipEffectFile")) {
      ceN$mergeStrands <- ce$mergeStrands;
      if (inherits(ce, "CnChipEffectFile")) {
        ceN$combineAlleles <- ce$combineAlleles;
      }
    }

    # CDF inheritance
    setCdf(ceN, cdf);

    verbose && enter(verbose, "Storing normalized signals");
    data[,"theta"] <- yN;
    # Not needed anymore
    yN <- NULL;
    updateDataFlat(ceN, data=data, verbose=less(verbose));
    # Not needed anymore
    data <- NULL;

    ## Create checksum file
    ceNZ <- getChecksumFile(ceN)

    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # for (kk in ...)

  # Create the output set
  outputSet <- getOutputDataSet(this, verbose=less(verbose,5));

  verbose && exit(verbose);

  outputSet;
})

############################################################################
# HISTORY:
# 2008-09-19
# o BUG FIX: process() of FragmentEquivalentClassNormalization did not
#   return a data set for which the sample attributes has been updated
#   according to optional sample annotation files (SAFs).
# o MEMORY OPTIMIZATION: process() no longer records each normalized array.
# o CLEANUP: process() no longer sets (unused) .outputSet field.
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: normalizeOneArrayVector().
# 2008-01-20
# o Created from FragmentLengthNormalization.R.
############################################################################
