###########################################################################/**
# @RdocClass QuantileNormalization
#
# @title "The QuantileNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization function that transforms the
#  probe-level signals towards the same empirical distribution.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "ProbeLevelTransform".}
#   \item{subsetToUpdate}{The probes to be updated.
#     If @NULL, all probes are updated.}
#   \item{typesToUpdate}{Types of probes to be updated.}
#   \item{targetDistribution}{A @numeric @vector.  The empirical
#     distribution to which all arrays should be normalized to.}
#   \item{subsetToAvg}{The probes to calculate average empirical
#     distribution over.  If a single @numeric in (0,1), then this
#     fraction of all probes will be used.
#     If @NULL, all probes are considered.}
#   \item{typesToAvg}{Types of probes to be used when calculating the
#     average empirical distribution.
#     If \code{"pm"} and \code{"mm"} only perfect-match and mismatch
#     probes are used, respectively. If \code{"pmmm"} both types are used.
#   }
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \examples{\dontrun{
#   @include "../incl/QuantileNormalization.Rex"
# }}
#
# @author "HB"
#*/###########################################################################
setConstructorS3("QuantileNormalization", function(..., subsetToUpdate=NULL, typesToUpdate=NULL, targetDistribution=NULL, subsetToAvg=subsetToUpdate, typesToAvg=typesToUpdate) {
  extend(ProbeLevelTransform(...), "QuantileNormalization",
    .subsetToUpdate = subsetToUpdate,
    .typesToUpdate = typesToUpdate,
    "cached:.targetDistribution" = targetDistribution,
    .subsetToAvg = subsetToAvg,
    .typesToAvg = typesToAvg
  )
})


setMethodS3("getSubsetToUpdate", "QuantileNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  subsetToUpdate <- this$.subsetToUpdate;

  # Done?
  if (identical(attr(subsetToUpdate, "adjusted"), TRUE))
    return(subsetToUpdate);

  # Ad hoc solution for ChipEffectSet:s for now. /HB 2007-04-11
  ds <- getInputDataSet(this);
  if (inherits(ds, "ChipEffectSet")) {
    verbose && enter(verbose, "Identifying possible cells in ", class(ds)[1]);
    df <- getOneFile(ds);
    possibleCells <- getCellIndices(df, verbose=less(verbose));
    possibleCells <- unlist(possibleCells, use.names=FALSE);
    possibleCells <- sort(possibleCells);
    verbose && str(verbose, possibleCells);

    verbose && cat(verbose, "'subsetToUpdate' (before): ");
    verbose && str(verbose, subsetToUpdate);

    if (is.null(subsetToUpdate)) {
      subsetToUpdate <- possibleCells;
    } else {
      subsetToUpdate <- intersect(subsetToUpdate, possibleCells);
    }
    # Not needed anymore
    possibleCells <- NULL;
    verbose && cat(verbose, "'subsetToUpdate' (after): ");
    verbose && str(verbose, subsetToUpdate);

    attr(subsetToUpdate, "adjusted") <- TRUE;
    this$.subsetToUpdate <- subsetToUpdate;

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);
  }

  subsetToUpdate;
}, private=TRUE)



setMethodS3("getExclCells", "QuantileNormalization", function(this, ...) {
  NULL;
}, private=TRUE)


setMethodS3("getSubsetToAvg", "QuantileNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  subsetToAvg <- this$.subsetToAvg;

  # Done?
  if (identical(attr(subsetToAvg, "adjusted"), TRUE))
    return(subsetToAvg);

  # Ad hoc solution for ChipEffectSet:s for now. /HB 2007-04-11
  ds <- getInputDataSet(this);
  if (inherits(ds, "ChipEffectSet")) {
    verbose && enter(verbose, "Identifying possible cells in ", class(ds)[1]);
    df <- getOneFile(ds);
    possibleCells <- getCellIndices(df, verbose=less(verbose));
    possibleCells <- unlist(possibleCells, use.names=FALSE);
    possibleCells <- sort(possibleCells);
    verbose && str(verbose, possibleCells);

    verbose && cat(verbose, "'subsetToAvg' (before): ");
    verbose && str(verbose, subsetToAvg);

    if (is.null(subsetToAvg)) {
      subsetToAvg <- possibleCells;
    } else {
      subsetToAvg <- intersect(subsetToAvg, possibleCells);
    }
    # Not needed anymore
    possibleCells <- NULL;
    verbose && cat(verbose, "'subsetToAvg' (after): ");
    verbose && str(verbose, subsetToAvg);

    attr(subsetToAvg, "adjusted") <- TRUE;
    this$.subsetToAvg <- subsetToAvg;

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);
  }

  subsetToAvg;
}, private=TRUE)




setMethodS3("getParameters", "QuantileNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  # Get parameters of this class
  params2 <- list(
    subsetToUpdate = getSubsetToUpdate(this),
    typesToUpdate = this$.typesToUpdate,
    subsetToAvg = getSubsetToAvg(this),
    typesToAvg = this$.typesToAvg,
    .targetDistribution = this$.targetDistribution
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, protected=TRUE)



setMethodS3("getTargetDistribution", "QuantileNormalization", function(this, sort=TRUE, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting target distribution");

  yTarget <- this$.targetDistribution;
  if (inherits(yTarget, "AffymetrixCelFile")) {
    df <- yTarget;
    verbose && enter(verbose, "Reading distribution from baseline array");
    verbose && cat(verbose, "Array: ", getFullName(df));
    yTarget <- getData(df, field="intensities")$intensities;
    this$.targetDistribution <- yTarget;
    verbose && exit(verbose);
  } else if (force || is.null(yTarget)) {
    pathname <- findTargetDistributionFile(this, verbose=less(verbose));
    verbose && print(verbose, pathname);

    if (isFile(pathname)) {
      verbose && enter(verbose, "Reading saved distribution: ", pathname);
      yTarget <- readApd(pathname)$quantiles;
      verbose && exit(verbose);
    } else {
      verbose && enter(verbose, "Calculating");
      yTarget <- calculateTargetDistribution(this, verbose=less(verbose));
      verbose && exit(verbose);
    }
    attr(yTarget, "identifier") <- getTargetDistributionIdentifier(this);
    this$.targetDistribution <- yTarget;
  } else {
    verbose && cat(verbose, "Was specified or cached in-memory.");
  }

  if (sort)
    yTarget <- sort(yTarget);
  verbose && str(verbose, yTarget);

  verbose && exit(verbose);

  yTarget;
}, private=TRUE)


setMethodS3("getTargetDistributionIdentifier", "QuantileNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting identifier for target distribution");

  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);
  params <- getParameters(this);

  # Get the parameters used for averaging
  indices <- params$subsetToAvg;

  # Speed up for digest() in case there are many indices
  nbrOfCells <- nbrOfCells(cdf);
  if (length(indices) > nbrOfCells/2) {
    indices <- -setdiff(1:nbrOfCells, indices);
  }

  key <- list(
    identifier=getIdentifier(ds),
    indices=indices,
    types=params$typesToAvg
  );
  id <- getChecksum(key);
  verbose && exit(verbose);

  id;
}, private=TRUE)

setMethodS3("getTargetDistributionPathname", "QuantileNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Getting the pathname for target distribution file to be created");

  ds <- getInputDataSet(this);
  path <- getPath(ds);

  if (getOption(aromaSettings, "devel/dropRootPathTags", TRUE)) {
    path <- dropRootPathTags(path, depth=2, verbose=less(verbose, 5));
  }
  verbose && cat(verbose, "Path without root-path tags: ", path);

  id <- getTargetDistributionIdentifier(this, verbose=less(verbose));
  filename <- sprintf(".averageQuantile-%s.apq", id);
  verbose && cat(verbose, "Filename: ", filename);

  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
  verbose && cat(verbose, "Pathname:");
  verbose && print(verbose, pathname);

  verbose && exit(verbose);

  pathname;
}, private=TRUE) # getTargetDistributionPathname()



setMethodS3("findTargetDistributionFile", "QuantileNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating the target distribution file");

  ds <- getInputDataSet(this);
  path <- getPath(ds);

  if (getOption(aromaSettings, "devel/dropRootPathTags", TRUE)) {
    path <- dropRootPathTags(path, depth=2, verbose=less(verbose, 5));
  }

  depth <- 2;

  # Search all possible root paths
  rootPath <- getParent(path, depth=depth);
  rootRootPath <- dirname(rootPath);
  rootPath <- basename(rootPath);
  pattern <- sprintf("^%s(|,.*)$", rootPath);
  rootPaths <- list.files(path=rootRootPath, pattern=pattern, full.names=FALSE);
  if (rootRootPath != ".") {
    rootPaths <- file.path(rootRootPath, rootPaths);
  }
  verbose && cat(verbose, "Root paths to be searched:");
  verbose && print(verbose, rootPaths);

  # Identify subdirectories
  subdirs <- sapply(seq_len(depth), FUN=function(d) {
    basename(getParent(path, depth=d-1L));
  });
  subdirs <- rev(subdirs);
  subdirs <- do.call(file.path, args=as.list(subdirs));
  verbose && cat(verbose, "Subdirectories: ", subdirs);

  id <- getTargetDistributionIdentifier(this, verbose=less(verbose));
  filename <- sprintf(".averageQuantile-%s.apq", id);
  verbose && cat(verbose, "Filename: ", filename);

  # All potential paths
  paths <- file.path(rootPaths, subdirs);
  verbose && cat(verbose, "Paths to be considered:");
  verbose && print(verbose, paths);

  # Keep only existing paths
  paths <- paths[sapply(paths, FUN=isDirectory)];
  verbose && cat(verbose, "Existing paths:");
  verbose && print(verbose, paths);

  pathnames <- sapply(paths, FUN=function(path) {
    Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
  });

  # Keep only existing pathnames
  pathnames <- pathnames[sapply(pathnames, FUN=isFile)];
  verbose && cat(verbose, "Existing pathnames:");
  verbose && print(verbose, pathnames);

  if (length(pathnames) > 0) {
    pathname <- pathnames[1];
    verbose && cat(verbose, "Keeping first pathname:");
    verbose && print(verbose, pathname);
  } else {
    verbose && cat(verbose, "Could not locate a matching pathname.");
    pathname <- NULL;
  }

  verbose && exit(verbose);

  pathname;
}, protected=TRUE) # findTargetDistributionFile()



setMethodS3("calculateTargetDistribution", "QuantileNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##   # Argument 'targetDataSet':
##   if (is.null(targetDataSet)) {
##     targetDataSet <- getInputDataSet(this);
##   } else if (inherits(targetDataSet, "AffymetrixCelSet")) {
##     cdf <- getCdf(targetDataSet);
##     dataSet <- getInputDataSet(this);
##     if (getChipType(cdf) != getChipType(getCdf(dataSet))) {
##       throw("Argument 'targetDataSet' does not have the same chip type as the input data set: ", getChipType(cdf), " != ", getChipType(getCdf(dataSet)));
##     }
##   } else {
##     throw("Argument 'targetDataSet' is not an AffymetrixCelSet: ",
##                                                   class(targetDataSet)[1]);
##   }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calculating target distribution");
  verbose && cat(verbose, "Method: average empirical distribution");

  # Get pathname where to store the target distribution
  pathname <- getTargetDistributionPathname(this, verbose=less(verbose));
  pathname <- Arguments$getWritablePathname(pathname);

  targetDataSet <- getInputDataSet(this);
  cdf <- getCdf(targetDataSet);

  params <- getParameters(this);

  cellsToSearch <- params$subsetToAvg;
  probes <- identifyCells(cdf, indices=cellsToSearch,
                         types=params$typesToAvg, verbose=less(verbose));
  # Not needed anymore
  params <- cellsToSearch <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && cat(verbose, "Using ", length(probes), " probes");

  # Exclude certain cells when *fitting* the normalization function?
  # If so, exclude them also when calculating the average distribution
  excl <- getExclCells(this);
  verbose && cat(verbose, "Excluding some cells when calculating the average distribution");
  verbose && str(verbose, excl);

  verbose && cat(verbose, "Calculating target distribution from the ",
                length(targetDataSet), " arrays in the input data set");
  # Calculate the average quantile
  yTarget <- averageQuantile(targetDataSet, probes=probes,
                              excludeCells=excl, verbose=less(verbose));
  # Not needed anymore
  probes <- NULL;

  # Write the result to file
  verbose && cat(verbose, "Saving distribution: ", pathname);
  writeApd(pathname, data=yTarget, name="quantiles");

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  invisible(yTarget);
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
setMethodS3("process", "QuantileNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Quantile normalizing data set");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already normalized");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve/calculate the target distribution
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving target distribution");
  getTargetDistribution(this, verbose=less(verbose));
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get algorithm parameters (including the target distribution above)
  params <- getParameters(this);

  # Get the output path
  outputPath <- getPath(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing data towards target distribution");
  names(params) <- gsub(".targetDistribution", "xTarget", names(params));
  args <- c(list(ds, path=outputPath), params);

  # Garbage collect
  # Not needed anymore
  params <- NULL; gc();

  verbose && cat(verbose, "Calling normalizeQuantile() with arguments:");
  verbose && str(verbose, args);
  args$verbose <- verbose;
  outputDataSet <- do.call(normalizeQuantile, args=args);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  # Update the output data set
  this$outputDataSet <- outputDataSet;

  verbose && exit(verbose);

  outputDataSet;
})

############################################################################
# HISTORY:
# 2011-02-24
# o Added findTargetDistributionFile() to QuantileNormalization for
#   locating an existing target-distribution file.  The previously used
#   getTargetDistributionPathname(), which returns a hardwired pathname,
#   is now only used for creating a target-distribution file.
# 2008-07-03
# o Now process() calls normalizeQuantileRank(), which is the new updated
#   name for normalizeQuantile().
# 2007-04-19
# o BUG FIX: Added missing getExclCells() to QuantileNormalization.
#   Thanks Elizabeth Purdom for the report.
# 2007-04-11
# o Added clearCache() for this class.
# 2007-02-04
# o Now QuantileNormalization() takes an AffymetrixCelFile as a target
#   distribution too, cf argument 'targetDistribution'.
# 2006-12-08
# o Now this class inherits from the ProbePreprocessor class.
# o Now this pre-processor output results to probeData/.
# o Renamed from QuantileNormalizer.
# 2006-11-18
# o Removed version and subversion tags, and related functions.
#   Now getTags() returns the tags of the input data set plus any tags
#   of this instance.
# 2006-10-30
# o Created.
############################################################################
