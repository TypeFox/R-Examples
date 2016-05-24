###########################################################################/**
# @RdocClass FragmentLengthNormalization
#
# @title "The FragmentLengthNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization method that corrects for PCR
#  fragment length effects on copy-number chip-effect estimates.
# }
#
# @synopsis
#
# \arguments{
#   \item{dataSet}{A @see "SnpChipEffectSet".}
#   \item{...}{Additional arguments passed to the constructor of
#     @see "ChipEffectTransform".}
#   \item{target}{(Optional) A @character string or a list of @functions
#     specifying what to normalize toward.
#     For each enzyme there is one target function to which all arrays
#     should be normalized to.}
#   \item{subsetToFit}{The units from which the normalization curve should
#     be estimated.  If @NULL, all are considered.}
#   \item{lengthRange}{If given, a @numeric @vector of length 2 specifying
#     the range of fragment lengths considered.  All fragments with lengths
#     outside this range are treated as if they were missing.}
#   \item{onMissing}{Specifies how to normalize units for which the
#     fragment lengths are unknown.}
#   \item{shift}{An optional amount the data points should be shifted
#      (translated).}
#   \item{targetFunctions}{Deprecated.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Requirements}{
#   This class requires a SNP information annotation file for the
#   chip type to be normalized.
# }
#
# \details{
#   For SNPs, the normalization function is estimated based on the total
#   chip effects, i.e. the sum of the allele signals.  The normalizing
#   is done by rescale the chip effects on the intensity scale such that
#   the mean of the total chip effects are the same across samples for
#   any given fragment length.  For allele-specific estimates, both alleles
#   are always rescaled by the same amount.  Thus, when normalizing
#   allele-specific chip effects, the total chip effects is change, but not
#   the relative allele signal, e.g. the allele B frequency.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("FragmentLengthNormalization", function(dataSet=NULL, ..., target=targetFunctions, subsetToFit="-XY", lengthRange=NULL, onMissing=c("median", "ignore"), shift=0, targetFunctions=NULL) {
  extraTags <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "SnpChipEffectSet");

#    if (dataSet$combineAlleles != TRUE) {
#      throw("Currently only total copy-number chip effects can be normalized, i.e. 'combineAlleles' must be TRUE");
#    }

#    dataSet <- Arguments$getInstanceOf(dataSet, "CnChipEffectSet");

#    if (dataSet$mergeStrands != TRUE) {
#      throw("Currently only non-strands specific copy-number chip effects can be normalized, i.e. 'mergeStrands' must be TRUE");
#    }
  }

  # Argument 'target':
  if (!is.null(target)) {
    if (is.character(target)) {
      if (target == "zero") {
      } else {
        throw("Unknown value of argument 'target': ", target);
      }
    } else if (is.list(target)) {
      # Validate each element
      for (kk in seq_along(target)) {
        if (!is.function(target[[kk]])) {
          throw("One element in 'target' is not a function: ",
                                          class(target[[kk]])[1]);
        }
      }
    } else {
      throw("Unknown value of argument 'target': ", class(target)[1]);
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

  # Argument 'shift':
  shift <- Arguments$getDouble(shift, disallow=c("NA", "NaN", "Inf"));

  # Argument 'lengthRange':
  if (!is.null(lengthRange)) {
    lengthRange <- Arguments$getDoubles(lengthRange);
    stopifnot(lengthRange[1] <= lengthRange[2]);
  }

  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);


  extend(ChipEffectTransform(dataSet, ...), "FragmentLengthNormalization",
    "cached:.targetFunctions" = NULL,
    "cached:.target" = target,
    .subsetToFit = subsetToFit,
    .lengthRange = lengthRange,
    .onMissing = onMissing,
    .extraTags = extraTags,
    shift = shift
  )
})


setMethodS3("getAsteriskTags", "FragmentLengthNormalization", function(this, collapse=NULL, ...) {
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



setMethodS3("getParameters", "FragmentLengthNormalization", function(this, expand=TRUE, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters", expand=expand);

  # Get parameters of this class
  params <- c(params, list(
    subsetToFit = this$.subsetToFit,
    lengthRange = this$.lengthRange,
    onMissing = this$.onMissing,
    .target = this$.target,
    shift = this$shift
  ));


  # Expand?
  if (expand) {
    subsetToFit <- getSubsetToFit(this);
  }

  params;
}, protected=TRUE)


setMethodS3("getCdf", "FragmentLengthNormalization", function(this, ...) {
  inputDataSet <- getInputDataSet(this);
  getCdf(inputDataSet);
})


setMethodS3("getOutputDataSet00", "FragmentLengthNormalization", function(this, ..., verbose=FALSE) {
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


setMethodS3("getFilteredFragmentLengths", "FragmentLengthNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading and filtering fragment lengths");

  verbose && enter(verbose, "Reading fragment lengths");
  # Get SNP information
  cdf <- getCdf(this);
  si <- getSnpInformation(cdf);
  verbose && print(verbose, si);

  fl <- getFragmentLengths(si, ...);
  verbose && cat(verbose, "Summary of non-filtered fragment lengths:");
  verbose && str(verbose, fl);
  verbose && summary(verbose, fl);
  verbose && exit(verbose);

  verbose && enter(verbose, "Filtering fragment lengths");
  # Get the range fragment lengths to be considered
  params <- getParameters(this, expand=FALSE);
  range <- params$lengthRange;

  if (!is.null(range)) {
    naValue <- as.double(NA);
    for (ee in seq_len(ncol(fl))) {
      flEE <- fl[,ee];
      ok <- (!is.na(flEE));

      tooSmall <- (ok & (flEE < range[1]));
      n <- sum(tooSmall);
      if (n > 0) {
        verbose && printf(verbose, "Detected %d fragments on enzyme %d that are too short (< %.0g)\n", n, ee, range[1]);
      }

      tooLarge <- (ok & (flEE > range[2]));
      n <- sum(tooLarge);
      if (n > 0) {
        verbose && printf(verbose, "Detected %d fragments on enzyme %d that are too long (> %.0g)\n", n, ee, range[2]);
      }

      tooExtreme <- (tooSmall | tooLarge);
      n <- sum(tooExtreme);

      verbose && printf(verbose, "Detected %d fragments on enzyme %d with lengths outside of filtered range [%.0g,%.0g]\n", n, ee, range[1], range[2]);

      fl[tooExtreme,ee] <- naValue;
    } # for (ee ...)

    verbose && cat(verbose, "Summary of filtered fragment lengths:");
    verbose && summary(verbose, fl);

    verbose && exit(verbose);
  } # if (!is.null(range))

  verbose && exit(verbose);

  fl;
}, protected=TRUE)


setMethodS3("getSubsetToFit", "FragmentLengthNormalization", function(this, force=FALSE, ..., verbose=FALSE) {
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

  # Get SNP information
  cdf <- getCdf(this);
  si <- getSnpInformation(cdf);
  verbose && print(verbose, si);

  # Identify all SNP and CN units (==potential units to be fitted)
  verbose && enter(verbose, "Identifying SNPs and CN probes");
## OLD:
## units <- indexOf(cdf, "^(SNP|CN)");
  types <- getUnitTypes(cdf, verbose=less(verbose,1));
  verbose && print(verbose, table(types));
  units <- which(types == 2 | types == 5);
  # Not needed anymore
  types <- NULL;
  verbose && cat(verbose, "units:");
  verbose && str(verbose, units);
  verbose && exit(verbose);

  onMissing <- this$.onMissing;
  if (onMissing == "ignore") {
    # Keep only those for which we have PCR fragmenth-length information
    # for at least one enzyme
    verbose && enter(verbose, "Identifying finite fragment lengths");
    fl <- getFilteredFragmentLengths(this, units=units, verbose=less(verbose,3));
    keep <- rep(FALSE, nrow(fl));
    for (ee in seq_len(ncol(fl))) {
      keep <- keep | is.finite(fl[,ee]);
    }
    units <- units[keep];
    verbose && printf(verbose, "Number of SNP/CN units without fragment-length details: %d out of %d (%.1f%%)\n", sum(!keep), length(keep), 100*sum(!keep)/length(keep));
    verbose && exit(verbose);
    # Not needed anymore
    fl <- NULL;
  }


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

    verbose && enter(verbose, "Reading fragment lengths");
    fl <- getFilteredFragmentLengths(this, units=units, verbose=less(verbose,3));
    verbose && exit(verbose);

    # Make sure to keep data points at the tails too
    extremeUnits <- c();
    for (ee in seq_len(ncol(fl))) {
      extremeUnits <- c(extremeUnits, which.min(fl[,ee]), which.max(fl[,ee]));
    }
    # Not needed anymore
    fl <- NULL;

    keep <- c(keep, extremeUnits);
    keep <- unique(keep);

    # Now filter
    units <- units[keep];
    # Not needed anymore
    keep <- NULL;
  }

  # Sort units
  units <- sort(units);

  # Assert correctness
  units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));

  # Cache
  this$.units <- units;

  verbose && exit(verbose);

  units;
}, private=TRUE)




setMethodS3("getTargetFunctions", "FragmentLengthNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  target <- this$.target;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Handling argument 'force'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # AD HOC: If forced, make sure to use the original target argument
  if (force) {
    if (!is.null(target)) {
      targetType <- attr(target, "targetType");
      if (!is.null(targetType))
        target <- targetType;
      # Not needed anymore
      targetType <- NULL;
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get SNP information
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);
  si <- getSnpInformation(cdf);
  verbose && print(verbose, si);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Predefined target functions?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.character(target)) {
    verbose && enter(verbose, "Setting up predefined target functions");
    targetType <- target;
    verbose && cat(verbose, "Target type: ", targetType);

    # Infer the number of enzymes
    fl <- getFilteredFragmentLengths(this, units=1:5);
    nbrOfEnzymes <- ncol(fl);
    # Not needed anymore
    fl <- NULL;

    if (identical(targetType, "zero")) {
      target <- rep(list(function(...) log2(2200)), nbrOfEnzymes);
    } else {
      throw("Unknown target function: ", targetType);
    }

    # Store the original target type in an attribute
    attr(target, "targetType") <- targetType;

    this$.target <- target;
    force <- FALSE;
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Target functions based on the average array?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (force || is.null(target)) {
    # Get target set
    ces <- getInputDataSet(this);
    verbose && enter(verbose, "Getting average signal across arrays");
    ceR <- getAverageFile(ces, force=force, verbose=less(verbose));
    # Not needed anymore
    ces <- NULL; # Not needed anymore
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    # Get units to fit
    units <- getSubsetToFit(this);
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);

    # Get target signals for SNPs
    yR <- extractTheta(ceR, units=units, verbose=less(verbose, 5));
    verbose && cat(verbose, "(Allele-specific) thetas:");
    verbose && str(verbose, yR);

    # If more than one theta per unit, sum them up to get the total signal
    if (ncol(yR) > 1) {
      # Row sums with na.rm=TRUE => NAs are treated as zeros.
      yR[is.na(yR)] <- 0;
      for (cc in 2:ncol(yR)) {
        yR[,1] <- yR[,1] + yR[,cc];
      }
    }
    yR <- yR[,1,drop=TRUE];
    verbose && cat(verbose, "Total thetas:");
    verbose && str(verbose, yR);

##     data <- getDataFlat(ceR, units=units, fields="theta", verbose=less(verbose));
    # Not needed anymore
    ceR <- NULL; # Not needed anymore
##     units <- data[,"unit"];
##     yR <- data[,"theta"];
##     # Not needed anymore
##     data <- NULL; # Not needed anymore

    verbose && cat(verbose, "Summary of total signals (on the intensity scale):");
    verbose && summary(verbose, yR);

    # Shift?
    shift <- this$shift;
    if (shift != 0) {
      verbose && printf(verbose, "Applying shift: %+f\n", shift);
      yR <- yR + shift;
      verbose && cat(verbose, "Summary of shifted signals (on the intensity scale):");
      verbose && summary(verbose, yR);
    }

    yR <- log2(yR);
    verbose && cat(verbose, "Signals:");
    verbose && str(verbose, yR);
    verbose && cat(verbose, "Summary of signals (on the log2 scale):");
    verbose && summary(verbose, yR);

    # Get PCR fragment lengths for these
    fl <- getFilteredFragmentLengths(this, units=units, verbose=less(verbose, 3));
    # Not needed anymore
    units <- NULL; # Not needed anymore

    verbose && cat(verbose, "Fragment lengths:");
    verbose && str(verbose, fl);
    verbose && cat(verbose, "Summary of fragment lengths:");
    verbose && summary(verbose, fl);

    # Fit lowess function
    verbose && enter(verbose, "Fitting target prediction function to each enzyme exclusively");
    okYR <- is.finite(yR);
    verbose && cat(verbose, "Distribution of log2 signals that are finite:");
    verbose && summary(verbose, okYR);

    hasFL <- is.finite(fl);
    verbose && cat(verbose, "Distribution of units with known fragment lengths:");
    verbose && summary(verbose, hasFL);

    nbrOfEnzymes <- ncol(fl);
    allEnzymes <- seq_len(nbrOfEnzymes);

    fits <- list();
    for (ee in allEnzymes) {
      verbose && enter(verbose, "Enzyme #", ee, " of ", nbrOfEnzymes);

      # Fit only to units with known length and non-missing data points.
      ok <- (hasFL[,ee] & okYR);

      verbose && cat(verbose, "Distribution of units with known fragment lengths and finite signals:");
      verbose && summary(verbose, ok);

      # Exclude multi-enzyme units
      for (ff in setdiff(allEnzymes, ee)) {
        ok <- ok & !hasFL[,ff];
      }

      verbose && cat(verbose, "Distribution of units with known fragment lengths and finite signals that are exclusively on this enzyme:");
      verbose && summary(verbose, ok);


      # Sanity check
      if (sum(ok) == 0) {
        throw(sprintf("Cannot fit target function to enzyme, because there are no (finite) data points that are unique to this enzyme (%d). Consider using argument 'lengthRange' when setting up the FragmentLengthNormalization.", ee));
      }

      # Fit fragment-length effect to single-enzyme units
      fit <- lowess(fl[ok,ee], yR[ok]);
      class(fit) <- "lowess";

      # Not needed anymore
      ok <- NULL;

      fits[[ee]] <- fit;

      # Not needed anymore
      fit <- NULL;

      verbose && exit(verbose);
    } # for (ee in allEnzymes)

    # Remove as many promises as possible
    # Not needed anymore
    target <- nbrOfEnzymes <- allEnzymes <- fl <- yR <- okYR <- hasFL <- NULL;

    # Create a target prediction function for each enzyme
    fcns <- vector("list", length(fits));
    for (ee in seq_along(fits)) {
      fcns[[ee]] <- function(x, ...) {
        predict(fits[[ee]], x, ...);  # Dispatched predict.lowess().
      }
    }
    verbose && str(verbose, fcns);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    target <- fcns;
    # Not needed anymore
    fcns <- NULL;
    this$.target <- target;
    force <- FALSE;

    verbose && exit(verbose);
  } # if (force || ...)


  target;
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
#   \item{...}{Additional arguments passed to
#     @see "aroma.light::normalizeFragmentLength" (only for advanced users).}
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
setMethodS3("process", "FragmentLengthNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Normalizing set for PCR fragment-length effects");

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

  verbose && enter(verbose, "Identifying SNP and CN units");
  # Get SNP & CN units
  cdf <- getCdf(ces);
## OLD:
##  subsetToUpdate <- indexOf(cdf, "^(SNP|CN)");
  types <- getUnitTypes(cdf, verbose=less(verbose,1));
  verbose && print(verbose, table(types));
  subsetToUpdate <- which(types == 2 | types == 5);
  # Not needed anymore
  types <- NULL;
  verbose && cat(verbose, "subsetToUpdate:");
  verbose && str(verbose, subsetToUpdate);
  verbose && exit(verbose);


  verbose && enter(verbose, "Retrieving SNP information annotations");
  si <- getSnpInformation(cdf);
  verbose && print(verbose, si);
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying the subset used to fit normalization function(s)");
  # Get subset to fit
  subsetToFit <- getSubsetToFit(this, verbose=less(verbose));
  verbose && str(verbose, subsetToFit);
  verbose && exit(verbose);

  shift <- this$shift;
  verbose && cat(verbose, "Shift: ", shift);

  onMissing <- this$.onMissing;
  verbose && cat(verbose, "onMissing: ", onMissing);

  # Get (and create) the output path
  path <- getPath(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fl <- NULL;
  targetFcns <- NULL;
#  map <- NULL;
  cellMatrixMap <- NULL;
  nbrOfArrays <- length(ces);

  res <- listenv()

  for (kk in seq_len(nbrOfArrays)) {
    ce <- ces[[kk]];
    verbose && enter(verbose, sprintf("Array #%d of %d ('%s')",
                                            kk, nbrOfArrays, getName(ce)));

    filename <- getFilename(ce);
    pathname <- filePath(path, filename);
    if (isFile(pathname)) {
      verbose && cat(verbose, "Already normalized. Skipping.")

      ## Assert validity of file
      ceN <- fromFile(ce, pathname)

      # Carry over parameters too.  AD HOC for now. /HB 2007-01-07
      if (inherits(ce, "SnpChipEffectFile")) {
        ceN$mergeStrands <- ce$mergeStrands
        if (inherits(ce, "CnChipEffectFile")) {
          ceN$combineAlleles <- ce$combineAlleles
        }
      }

      # CDF inheritance
      setCdf(ceN, cdf)

      res[[kk]] <- pathname

      verbose && exit(verbose)
      next
    }

    # Get unit-to-cell? (for optimized reading)
#    if (is.null(map)) {
#      # Only loaded if really needed.
#      verbose && enter(verbose, "Retrieving unit-to-cell map for all arrays");
#      map <- getUnitGroupCellMap(ce, units=subsetToUpdate, verbose=less(verbose));
#      verbose && str(verbose, map);
#      verbose && exit(verbose);
#    }

    if (is.null(fl)) {
      # For the subset to be fitted, get PCR fragment lengths (for all enzymes)
      fl <- getFilteredFragmentLengths(this, units=subsetToUpdate, verbose=less(verbose, 3));
      verbose && summary(verbose, fl);

      # Get the index in the data vector of subset to be fitted.
      # Note: match() only returns first match, which is why we do
      # it this way.
      subset <- match(subsetToUpdate, subsetToFit);
      subset <- subset[!is.na(subset)];
      subset <- match(subsetToFit[subset], subsetToUpdate);
      verbose && str(verbose, subset);
    }

    if (is.null(targetFcns)) {
      # Only loaded if really needed.
      # Retrieve/calculate the target function
      targetFcns <- getTargetFunctions(this, verbose=less(verbose));
    }

    if (is.null(cellMatrixMap)) {
      verbose && enter(verbose, "Getting cell matrix map");
      cellMatrixMap <- getUnitGroupCellMatrixMap(ce, units=subsetToUpdate, verbose=less(verbose, 10));
      verbose && str(verbose, cellMatrixMap);
      verbose && exit(verbose);
    }


    res[[kk]] %<=% {
      # Get target log2 signals for all SNPs to be updated
      verbose && enter(verbose, "Getting theta estimates");
      theta <- extractTheta(ce, units=cellMatrixMap, drop=FALSE, verbose=less(verbose, 5));
      verbose && str(verbose, theta);
      verbose && summary(verbose, theta);
      verbose && exit(verbose);

      verbose && enter(verbose, "Calculating total signals");
      # Get the total locus signals?
      if (ncol(theta) > 1) {
        # Row sums with na.rm=TRUE => NAs are treated as zeros.
        y <- theta;
        y[is.na(y)] <- 0;
        for (cc in 2:ncol(y)) {
          y[,1] <- y[,1] + y[,cc];
        }
        y <- y[,1,drop=TRUE];
      } else {
        y <- theta[,1,drop=TRUE];
      }
      verbose && cat(verbose, "Total thetas:");
      verbose && str(verbose, y);
      verbose && exit(verbose);

  #    data <- getDataFlat(ce, units=map, fields="theta", verbose=less(verbose));
  #    verbose && str(verbose, data);
  #    y0 <- data[,"theta",drop=TRUE];
  #    stopifnot(identical(y,y0));
  #    verbose && str(verbose, y);
  #    verbose && exit(verbose);


      # Extract the values to fit the normalization function
      verbose && enter(verbose, "Normalizing log2 signals");

      # Shift?
      if (shift != 0)
        y <- y + shift;

      # Fit on the log2 scale
      y <- log2(y);

      verbose && cat(verbose, "Log2 signals:");
      verbose && str(verbose, y);
      yN <- .normalizeFragmentLength(y, fragmentLengths=fl,
                      targetFcns=targetFcns, subsetToFit=subset,
                      onMissing=onMissing, ...);
      verbose && cat(verbose, "Normalized log2 signals:");
      verbose && summary(verbose, yN);

      # Normalization scale factor for each unit (on the log2 scale)
      rho <- y-yN;
      # Not needed anymore
      y <- yN <- NULL;
      # On the intensity scale
      rho <- 2^rho;
      verbose && cat(verbose, "Normalization scale factors:");
      verbose && summary(verbose, rho);

      # Sanity check
      stopifnot(length(rho) == nrow(theta));

      # Normalize the theta:s (on the intensity scale)
      ok <- which(is.finite(rho));
      verbose && str(verbose, ok);
      theta[ok,] <- theta[ok,]/rho[ok];
      # Not needed anymore
      ok <- rho <- NULL;

      verbose && cat(verbose, "Normalized thetas:");
      verbose && str(verbose, theta);
      verbose && summary(verbose, theta);

      verbose && exit(verbose);

      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");

      # Write to a temporary file (allow rename of existing one if forced)
      isFile <- isFile(pathname);
      pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose);

      ceN <- createFrom(ce, filename=pathnameT, path=NULL, verbose=less(verbose));
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
  #    data[,"theta"] <- yN;
  #    # Not needed anymore
  #    yN <- NULL;
  #    updateDataFlat(ceN, data=data, verbose=less(verbose));
  #    # Not needed anymore
  #    data <- NULL;
      ok <- which(is.finite(cellMatrixMap));
      cells <- cellMatrixMap[ok];
      data <- theta[ok];
      # Not needed anymore
      ok <- theta <- NULL;

      verbose2 <- -as.integer(verbose) - 5;
      pathnameN <- getPathname(ceN);
      .updateCel(pathnameN, indices=cells, intensities=data, verbose=verbose2);
      # Not needed anymore
      cells <- data <- ceN <- NULL;
      verbose && exit(verbose);

      # Rename temporary file
      popTemporaryFile(pathnameT, verbose=verbose);

      ## Create checksum file
      dfZ <- getChecksumFile(pathname)

      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);

      pathname
    } ## %<=%

    verbose && exit(verbose);
  } # for (kk in ...)

  fl <- NULL  ## Not needed anymore

  ## Resolve futures
  res <- as.list(res)
  res <- NULL

  # Garbage collect
#  clearCache(this);
  gc <- gc();
  verbose && print(verbose, gc);

  # Create the output set
  outputSet <- getOutputDataSet(this, verbose=less(verbose,5));

  verbose && exit(verbose);

  outputSet;
})

############################################################################
# HISTORY:
# 2011-03-28
# o CLARIFICATION: Now the error message thrown when there are not
#   enough data points for a unique enzyme gives a hint on what to do.
# 2011-02-08
# o GENERALIZATION: Now it is possible to specify the range of fragment
#   lengths to be considered when normalizing for PCR fragment-length
#   effects.  See argument 'lengthRange' to FragmentLengthNormalization.
# 2010-02-15
## o MEMORY OPTIMIZATION: Now process() of FragmentLengthNormalization
##   clears the in-memory cache when finished.
# 2008-12-03
# o BUG FIX: Missing 'si' object.
# 2008-12-01
# o BUG FIX: For allele-specific estimates, FragmentLengthNormalization
#   would correctly estimate normalization scale factors, but due to a
#   typo, it effectively only update the signals for allele A.
#   Looking at the SVN history, this has always been the case.
# 2008-11-28
# o Now constructor argument 'targetFunctions' can also be "zero".
# 2008-09-19
# o BUG FIX: process() of FragmentLengthNormalization did not return a
#   data set for which the sample attributes has been updated according
#   to optional sample annotation files (SAFs).
# o MEMORY OPTIMIZATION: process() no longer records each normalized array.
# o CLEANUP: process() no longer sets (unused) .outputSet field.
# 2008-09-12
# o Added argument 'onMissing' to FragmentLengthNormalization, which is
#   passed down to normalizeFragmentLength() [req aroma.light v1.9.2] to
#   make it possible to normalize also units for which fragment lenghts
#   are unknown.
# 2008-07-07
# o Typo: The constructor error message for validating 'dataSet' outputted
#   vector of messages because the whole class vector was pasted.
# 2008-06-09
# o Updated process() to normalize allele-specific estimates as well.
#   We note that the different alleles are rescaled with the same factor,
#   that is, it is only the total chip effect that is changed whereas
#   for instance freqB = thetaB/(thetaA+thetaB) remains constant.
# 2008-06-06
# o getTargetFunctions() now sums allele-specific thetas and fits the
#   normalization function on the total thetas, if allele specific.
# o Now getTargetFunctions() utilizes extractTheta() and no longer
#   getDataFlat().
# 2008-03-29
# o Added more verbose output for the getTargetFunctions() in order to
#   simplify troubleshooting.
# 2008-02-18
# o Added 'shift' to FragmentLengthNormalization, cf ProbeLevelModel.
# 2007-12-01
# o Added getAsteriskTag() to FragmentLengthNormalization.
# o Similar to AllelicCrosstalkCalibration, the constructor argument
#   'subsetToFit' of FragmentLengthNormalization accept "-XY" (and "-X" and
#   "-Y") to specify the set of units to fit the model over to be all units
#   that are not on ChrX or ChrY.
# 2007-11-19
# o Updated getSubsetToFit() to handle chip types with multiple enzymes.
# o Updated methods to give an error if chip types with more than one
#   enzyme is tried to be normalized.
# 2007-09-16
# o Added clearCache() to FragmentLengthNormalization such that cached
#   target distributions can be cleared.
# 2007-09-12
# o BUG FIX: getSubsetToFit() of FragmentLengthNormalization would only
#   return SNP units, but not CN units which are available on the newer
#   chip types.  Similarly, process() would only update SNPs, but not
#   CN units.
# o Now getOutputDataSet() of FragmentLengthNormalization set and pass down
#   'mergeStrands' and 'combineAlleles' to ditto of the super class, if
#   applicable.  This way we avoid having to infer those arguments from
#   the contents of the files.
# 2007-02-20
# o Now FragmentLengthNormalization should handle cases with more than one
#   chip effect per unit, e.g. when mergeStrands=FALSE.
# 2007-01-16
# o BUG FIX: Forgot to clear the cache after cloning data set in process().
#   This would cause getAverage() to return a cached averaged from the
#   non-normalized data set.
# 2007-01-07
# o Now chip-effect parameters are carried over to the output set too.
# o BUG FIX: process(): Forgot to skip to next array in for loop if an
#   array was detected to be already normalized. Generated a "file already
#   exists" error.
# o Added garbage collection after each array been normalized.
# 2007-01-04
# o BUG FIX: process() gave an error if the data set was already done.
# 2006-12-08
# o Now this class inherits from the ChipEffectPreprocessing class.
# o Now this pre-processor output results to plmData/.
# 2006-11-28
# o Created from QuantileNormalizer.R.
############################################################################
