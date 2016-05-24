###########################################################################/**
# @RdocClass AdditiveCovariatesNormalization
#
# @title "The AdditiveCovariatesNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization method that corrects for
#  GC-content effects on copy-number chip-effect estimates.
# }
#
# @synopsis
#
# \arguments{
#   \item{dataSet}{A @see "SnpChipEffectSet".}
#   \item{...}{Additional arguments passed to the constructor of
#     @see "ChipEffectTransform".}
#   \item{target}{(Optional) A @character string or a @function
#     specifying what to normalize toward.}
#   \item{subsetToFit}{The units from which the normalization curve should
#     be estimated.  If @NULL, all are considered.}
#   \item{onMissing}{Specifies how to normalize units for which the
#     GC contents are unknown.}
#   \item{shift}{An optional amount the data points should be shifted
#      (translated).}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   For SNPs, the normalization function is estimated based on the total
#   chip effects, i.e. the sum of the allele signals.  The normalizing
#   is done by rescale the chip effects on the intensity scale such that
#   the mean of the total chip effects are the same across samples for
#   any given GC content.  For allele-specific estimates, both alleles
#   are always rescaled by the same amount.  Thus, when normalizing
#   allele-specific chip effects, the total chip effects is change, but not
#   the relative allele signal, e.g. the allele B frequency.
# }
#
# @author
#*/###########################################################################
setConstructorS3("AdditiveCovariatesNormalization", function(dataSet=NULL, ..., target=NULL, subsetToFit="-XY", shift=0, onMissing=c("median", "ignore")) {
  extraTags <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSets <- Arguments$getInstanceOf(dataSet, "SnpChipEffectSet");
  }

  # Argument 'target':
  if (!is.null(target)) {
    if (is.function(target)) {
      target <- list(target);
    }
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
    unf <- getUnitNamesFile(dataSet);
    subsetToFit <- Arguments$getIndices(subsetToFit, max=nbrOfUnits(unf));
    subsetToFit <- unique(subsetToFit);
    subsetToFit <- sort(subsetToFit);
  }

  # Argument 'shift':
  shift <- Arguments$getDouble(shift, disallow=c("NA", "NaN", "Inf"));

  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);


  extend(ChipEffectTransform(dataSet, ...), "AdditiveCovariatesNormalization",
    .subsetToFit = subsetToFit,
    "cached:.target" = target,
    .onMissing = onMissing,
    .extraTags = extraTags,
    shift = shift
  )
})



setMethodS3("getAsteriskTags", "AdditiveCovariatesNormalization", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=collapse);

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


setMethodS3("getParameters", "AdditiveCovariatesNormalization", function(this, expand=TRUE, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters", expand=expand);

  # Get parameters of this class
  params <- c(params, list(
    subsetToFit = this$.subsetToFit,
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


setMethodS3("getCdf", "AdditiveCovariatesNormalization", function(this, ...) {
  inputDataSet <- getInputDataSet(this);
  getCdf(inputDataSet);
})


setMethodS3("getOutputDataSet00", "AdditiveCovariatesNormalization", function(this, ..., verbose=FALSE) {
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

  args <- list(generic="getOutputDataSet", object=this, ...);

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

setMethodS3("getCovariates", "AdditiveCovariatesNormalization", abstract=TRUE, protected=TRUE);


setMethodS3("getSubsetToFit", "AdditiveCovariatesNormalization", function(this, force=FALSE, ..., verbose=FALSE) {
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
  nbrOfUnits <- nbrOfUnits(cdf);

  # Identify all SNP and CN units (==potential units to be fitted)
  verbose && enter(verbose, "Identifying SNPs and CN probes");
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
    # Keep only those for which we have annotation data covariates on
    # for at least one covariate
    verbose && enter(verbose, "Reading annotation data covariates");
    X <- getCovariates(this, units=units, verbose=less(verbose,5));
    keep <- rep(FALSE, nrow(X));
    for (ee in seq_len(ncol(X))) {
      keep <- keep | is.finite(X[,ee]);
    }
    units <- units[keep];
    verbose && printf(verbose, "Number of SNP/CN units without finite covariates: %d out of %d (%.1f%%)\n", sum(!keep), length(keep), 100*sum(!keep)/length(keep));
    verbose && exit(verbose);
    # Not needed anymore
    keep <- X <- NULL;
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
        subset <- setdiff(seq_len(nbrOfUnits), subset);

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

    verbose && enter(verbose, "Reading annotation data covariates");
    X <- getCovariates(this, units=units, verbose=less(verbose,5));
    verbose && str(verbose, X);
    verbose && exit(verbose);

    # Make sure to keep data points at the tails too
    extremeUnits <- c();
    for (ee in seq_len(ncol(X))) {
      extremeUnits <- c(extremeUnits, which.min(X[,ee]), which.max(X[,ee]));
    }
    # Not needed anymore
    X <- NULL;

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
  units <- Arguments$getIndices(units, max=nbrOfUnits);

  # Cache
  this$.units <- units;

  verbose && exit(verbose);

  units;
}, private=TRUE)




setMethodS3("getTargetFunctions", "AdditiveCovariatesNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
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
  # Predefined target functions?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.character(target)) {
    verbose && enter(verbose, "Setting up predefined target functions");
    targetType <- target;
    verbose && cat(verbose, "Target type: ", targetType);

    # Infer the number of covariates
    X <- getCovariates(this, units=1:5);
    nbrOfCovariates <- ncol(X);
    # Not needed anymore
    X <- NULL;

    if (identical(targetType, "zero")) {
      target <- rep(list(function(...) log2(2200)), nbrOfCovariates);
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
    # Not needed anymore
    ceR <- NULL; # Not needed anymore
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
    X <- getCovariates(this, units=units, verbose=less(verbose,5));
    # Not needed anymore
    units <- NULL; # Not needed anymore
    verbose && cat(verbose, "Annotation data covariates:");
    verbose && str(verbose, X);
    verbose && cat(verbose, "Summary of covariates:");
    verbose && summary(verbose, X);

    # Fit lowess function
    verbose && enter(verbose, "Fitting target prediction function to each covariate exclusively");
    okYR <- is.finite(yR);
    verbose && cat(verbose, "Distribution of log2 signals that are finite:");
    verbose && summary(verbose, okYR);

    hasX <- is.finite(X);
    verbose && cat(verbose, "Distribution of units with known covariates:");
    verbose && summary(verbose, hasX);

    nbrOfCovariates <- ncol(X);
    allCovariates <- seq_len(nbrOfCovariates);

    fits <- list();
    for (ee in allCovariates) {
      verbose && enter(verbose, "Covariate #", ee, " of ", nbrOfCovariates);

      # Fit only to units with known length and non-missing data points.
      ok <- (hasX[,ee] & okYR);

      verbose && cat(verbose, "Distribution of units with known covariates and finite signals:");
      verbose && summary(verbose, ok);

      # Exclude multi-covariate units
      for (ff in setdiff(allCovariates, ee)) {
        ok <- ok & !hasX[,ff];
      }

      verbose && cat(verbose, "Distribution of units with known covariates and finite signals that are exclusively for this covarite:");
      verbose && summary(verbose, ok);


      # Sanity check
      if (sum(ok) == 0) {
        throw("Cannot fit target function to covariate, because there are no (finite) data points that are unique to this covariate: ", ee);
      }

      # Fit effect to single-covariate units
      fit <- lowess(X[ok,ee], yR[ok]);
      class(fit) <- "lowess";

      # Not needed anymore
      ok <- NULL;

      fits[[ee]] <- fit;

      # Not needed anymore
      fit <- NULL;

      verbose && exit(verbose);
    } # for (ee ...)

    # Remove as many promises as possible
    # Not needed anymore
    target <- nbrOfCovariates <- allCovariates <- X <- hasX <- yR <- okYR <- NULL;

    # Create a target prediction function for each covariate
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
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "AdditiveCovariatesNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
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
  types <- getUnitTypes(cdf, verbose=less(verbose,1));
  verbose && print(verbose, table(types));
  subsetToUpdate <- which(types == 2 | types == 5);
  # Not needed anymore
  types <- NULL;
  verbose && cat(verbose, "subsetToUpdate:");
  verbose && str(verbose, subsetToUpdate);
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
  X <- NULL;
  targetFcns <- NULL;
#  map <- NULL;
  cellMatrixMap <- NULL;
  nbrOfArrays <- length(ces);
  for (kk in seq_len(nbrOfArrays)) {
    ce <- ces[[kk]];
    verbose && enter(verbose, sprintf("Array #%d of %d ('%s')",
                                            kk, nbrOfArrays, getName(ce)));

    filename <- getFilename(ce);
    pathname <- filePath(path, filename);
    if (!force && isFile(pathname)) {
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

    # Get unit-to-cell? (for optimized reading)
#    if (is.null(map)) {
#      # Only loaded if really needed.
#      verbose && enter(verbose, "Retrieving unit-to-cell map for all arrays");
#      map <- getUnitGroupCellMap(ce, units=subsetToUpdate, verbose=less(verbose));
#      verbose && str(verbose, map);
#      verbose && exit(verbose);
#    }

    if (is.null(X)) {
      # For the subset to be fitted, get PCR fragment lengths (for all covariates)
      X <- getCovariates(this, units=subsetToUpdate, verbose=less(verbose,5));
      verbose && summary(verbose, X);

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

    # Extract the values to fit the normalization function
    verbose && enter(verbose, "Normalizing log2 signals");

    # Shift?
    if (shift != 0)
      y <- y + shift;

    # Fit on the log2 scale
    y <- log2(y);

    verbose && cat(verbose, "Log2 signals:");
    verbose && str(verbose, y);
    yN <- .normalizeFragmentLength(y, fragmentLengths=X,
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

    # Write to a temporary file (allow rename of existing one if forced)
    isFile <- (force && isFile(pathname));
    pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose);

    # Create CEL file to store results, if missing
    verbose && enter(verbose, "Creating CEL file for results, if missing");
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
    ok <- which(is.finite(cellMatrixMap));
    cells <- cellMatrixMap[ok];
    data <- theta[ok];
    # Not needed anymore
    ok <- theta <- NULL;

    verbose2 <- -as.integer(verbose) - 5;
    pathnameN <- getPathname(ceN);
    .updateCel(pathnameN, indices=cells, intensities=data, verbose=verbose2);
    # Not needed anymore
    cells <- data <- NULL;
    verbose && exit(verbose);

    # Rename temporary file
    popTemporaryFile(pathnameT, verbose=verbose);

    ## Create checksum
    ceNZ <- getChecksumFile(pathname)

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
# 2013-06-01
# o Argument process(..., force=TRUE) for AdditiveCovariatesNormalization
#   did not force reprocessing.
# 2009-03-22
# o Generalized by creating abstract getCovariates() method.
# o Created from FragmentLengthNormalization.R.
############################################################################
