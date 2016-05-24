###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod normalizeQuantile
#
# @title "Normalizes samples to have the same empirical distribution"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The path where to save the normalized data files.
#     If @NULL, a default name is used.}
#   \item{name}{The name of the normalized data set, which will also be
#     part of the default path.}
#   \item{subsetToUpdate}{The probes to be updated.
#     If @NULL, all probes are updated.}
#   \item{typesToUpdate}{Types of probes to be updated.}
#   \item{xTarget}{A @numeric @vector.  The empirical distribution
#     to which all arrays should be normalized to.}
#   \item{subsetToAvg}{The probes to calculate average empirical
#     distribution over.  If a single @numeric in (0,1), then this
#     fraction of all probes will be used.
#     If @NULL, all probes are considered.}
#   \item{typesToAvg}{Types of probes to be used when calculating the
#     average empirical distribution.
#     If \code{"pm"} and \code{"mm"} only perfect-match and mismatch
#     probes are used, respectively. If \code{"pmmm"} both types are used.
#   }
#   \item{...}{Additional arguments passed to \code{normalizeQuantile()}.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# @author "HB"
#
# \seealso{
#   @see "aroma.light::normalizeQuantileRank.numeric"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("normalizeQuantile", "AffymetrixCelSet", function(this, path=NULL, name="normQuantile", subsetToUpdate=NULL, typesToUpdate=NULL, xTarget=NULL, subsetToAvg=subsetToUpdate, typesToAvg=typesToUpdate, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'path':
  if (is.null(path)) {
    # Path structure: /normQuantile/<data set name>/chip_data/<chip type>/
    path <- file.path(name, getName(this), "chip_data", getChipType(cdf));
  }
  if (!is.null(path)) {
    # Verify this path (and create if missing)
    path <- Arguments$getWritablePath(path);
  }

  if (identical(getPath(this), path)) {
    throw("Cannot calibrate data file. Argument 'path' refers to the same path as the path of the data file to be calibrated: ", path);
  }

  # Argument 'xTarget':
  if (is.null(xTarget)) {
    throw("DEPRECATED: normalizeQuantile() must no longer be called with xTarget=NULL.");
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the subset of probes to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying the probes to be updated");
  subsetToUpdate <- identifyCells(cdf, indices=subsetToUpdate,
                                                     types=typesToUpdate);
  verbose && exit(verbose);

  verbose && cat(verbose, "Normalizing ", length(subsetToUpdate), " probes");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(this)
  verbose && enter(verbose, "Normalizing ", nbrOfArrays, " arrays");

  dataFiles <- listenv()

  for (kk in seq_len(nbrOfArrays)) {
    df <- this[[kk]]
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, getName(df), nbrOfArrays))

    dataFiles[[kk]] %<=% {
      verbose && print(verbose, df)
      normalizeQuantile(df, path=path,
                        subsetToUpdate=subsetToUpdate, typesToUpdate=NULL,
                        xTarget=xTarget, ..., verbose=less(verbose))
    }

    # Garbage collect
    gc()

    verbose && exit(verbose)
  }

  ## Resolve futures
  dataFiles <- as.list(dataFiles)

  verbose && exit(verbose)

  ## Setup output data set
  res <- newInstance(this, dataFiles)
  setCdf(res, cdf)

  res
}, protected=TRUE) # normalizeQuantile()




###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod averageQuantile
#
# @title "Calculates the average empirical distribution across all samples"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{probes}{An optional @numeric @vector specifying what subset of
#      probes to be used to calculate the empirical distribution.
#      If @NULL, all probes are used.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector of length N.
# }
#
# \section{Missing values}{
#   If @NAs are detected in a sample, these are excluded and the
#   \code{approx()} function (@see "stats::approx") is used to "expand"
#   the @vector of the remaining values so that the sorted @vector
#   (still) have length N.
# }
#
# \details{
#   This methods implements Step A2-A3 in the algorithm for quantile
#   normalization proposed by Bengtsson et al. (2008).
# }
#
# @author "HB"
#
# \references{
#   [1] H. Bengtsson, R. Irizarry, B. Carvalho, & T.P. Speed.
#       Estimation and assessment of raw copy numbers at the single
#       locus level, Bioinformatics, 2008.
# }
#
# \seealso{
#   @see "aroma.light::averageQuantile.list"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("averageQuantile", "AffymetrixCelSet", function(this, probes=NULL, excludeCells=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'probes':
  probes <- identifyCells(cdf, indices=probes);  # TODO!
  # "TODO" since when? ;) /HB 2007-04-11

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arrays of interest
  arrays <- getNames(this);
  nbrOfChannels <- length(arrays);

  if (is.null(probes)) {
    nbrOfObservations <- nbrOfCells(this);
  } else {
    nbrOfObservations <- length(probes);
  }

  # Construct the sample quantiles
  quantiles <- (0:(nbrOfObservations-1))/(nbrOfObservations-1);

  # Create a vector to hold the target distribution
  xTarget <- vector("double", length=nbrOfObservations);

  verbose && enter(verbose, "Calculating the average empircal distribution across ", nbrOfChannels, " arrays");

  verbose && printf(verbose, "Number of probes: %d (%.1f%%)\n",
                   nbrOfObservations, 100*nbrOfObservations/nbrOfCells(cdf));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the sample quantile for all channels (columns)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (cc in 1:nbrOfChannels) {
    verbose && enter(verbose, "Array #", cc);

    verbose && printf(verbose, "reading...\n");
    df <- this[[cc]];
    Xcc <- getData(df, indices=probes, fields="intensities", ..., verbose=less(verbose, 2));
    Xcc <- as.vector(Xcc$intensities);
#    verbose && str(verbose, Xcc);

    # Exclude cells?
    if (!is.null(excludeCells))
      Xcc[excludeCells] <- NA;

    # Garbage collect
    gc();

    # Order and sort the values
    verbose && printf(verbose, "sorting...\n");
    Scc <- sort(Xcc);
#    verbose && str(verbose, Scc);

    # Garbage collect
    gc();

    # The number of non-NA observations
    nobs <- length(Scc);

    # Has NAs?
    nbrOfNAs <- (nbrOfObservations - nobs);
    if(nbrOfNAs > 0) {
      verbose && printf(verbose, "Detected %d NAs (%.2f%%),\n",
                           nbrOfNAs, 100*nbrOfNAs/nbrOfObservations);
      tt <- !is.na(Xcc);  # TODO?!? /HB 2006-07-22
      # Not needed anymore
      Xcc <- tt <- NULL;

      # Get the sample quantiles for those values
      bins <- (0:(nobs-1))/(nobs-1);

      # Interpolate to get the values at positions specified by
      # 'quantile' using data points given by 'bins' and 'Scc'.
      Scc <- approx(x=bins, y=Scc, xout=quantiles, ties="ordered")$y;
      # Not needed anymore
      bins <- NULL;
    } else {
      # Not needed anymore
      Xcc <- NULL;
    }

    # Incremental mean
    verbose && printf(verbose, "summing...\n");
    xTarget <- xTarget + Scc;
    # Not needed anymore
    Scc <- NULL;

    # Garbage collect
    gc();

    verbose && exit(verbose);
  }

  xTarget <- xTarget/nbrOfChannels;

  verbose && exit(verbose);

  xTarget;
}, protected=TRUE) # averageQuantile()





setMethodS3("transformAffine", "AffymetrixCelSet", function(this, outPath=file.path("transformed", getChipType(getCdf(this))), offset=0, scale=1, subsetToUpdate=NULL, typesToUpdate=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'outPath':
  outPath <- Arguments$getReadablePathname(outPath, mustExist=FALSE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'offset':
  offset <- Arguments$getDouble(offset);

  # Argument 'scale':
  scale <- Arguments$getDouble(scale, range=c(0,Inf));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the subset of probes to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying the probes to be updated");
  cdf <- getCdf(this);
  subsetToUpdate <- identifyCells(cdf, indices=subsetToUpdate,
                                                     types=typesToUpdate);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Transforming ", length(subsetToUpdate),
                               " probes on ", length(this), " arrays");
  dataFiles <- list();
  for (kk in seq_along(this)) {
    df <- this[[kk]];
    verbose && enter(verbose, "Array #", kk, " (", getName(df), ")");
    dataFiles[[kk]] <- transformAffine(df, outPath=outPath,
                  offset=offset, scale=scale, subsetToUpdate=subsetToUpdate,
                                                ..., verbose=less(verbose));
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  # CDF inheritance
  res <- newInstance(this, dataFiles);
  setCdf(res, cdf);
  res;
}, private=TRUE) # transformAffine()


############################################################################
# HISTORY:
# 2012-10-21 [HB]
# o CLEANUP: Dropped unneeded mkdirs(), because they were all preceeded
#   by an Arguments$getWritablePath().
# 2007-04-11
# o Added more verbose output to averageQuantile() in the case when NAs
#   are detected.  Added some Rdoc comments on how NAs are handles.
# 2006-09-15
# o Modified some argument names for normalizeQuantile().
# 2006-09-14
# o Recreated from old AffymetrixDataSet.NORM.R.
# 2006-07-27
# o Added transformAffine().
# o BUG FIX: The 'outPath' argument of normalizeQuantile() in the
#   AffymetrixDataSet class was not recognized.
# 2006-05-15
# o Extracted from AffymetrixDataSet.R.
# 2006-03-18
# o Added argument 'subset' to calcAvgProbeSignals() & normalizeQuantile().
############################################################################
