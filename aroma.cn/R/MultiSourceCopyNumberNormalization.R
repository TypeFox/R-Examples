###########################################################################/**
# @RdocClass MultiSourceCopyNumberNormalization
#
# @title "The MultiSourceCopyNumberNormalization class"
#
# \description{
#  @classhierarchy
#
#  The multi-source copy-number normalization (MSCN) method [1] is a
#  normalization method that normalizes copy-number estimates measured
#  by multiple sites and/or platforms for common samples.  It normalizes the
#  estimates toward a common scale such that for any copy-number level
#  the mean level of the normalized data are the same.
# }
#
# @synopsis
#
# \arguments{
#  \item{dsList}{A @list of K @see "aroma.core::AromaUnitTotalCnBinarySet":s.}
#  \item{fitUgp}{An @see "aroma.core::AromaUgpFile" that specifies the
#    common set of loci used to normalize the data sets at.}
#  \item{subsetToFit}{The subset of loci (as mapped by the \code{fitUgp}
#    object) to be used to fit the normalization functions.
#    If @NULL, loci on chromosomes 1-22 are used, but not on ChrX and ChrY.
#  }
#  \item{targetDimension}{A @numeric index specifying the data set in
#    \code{dsList} to which each platform in standardize towards.
#    If @NULL, the arbitrary scale along the fitted principal curve
#    is used.  This always starts at zero and increases.}
#  \item{align}{A @character specifying type of alignment applied, if any.
#    If \code{"none"}, no alignment is done.
#    If \code{"byChromosome"}, the signals are shifted chromosome
#    by chromosome such the corresponding smoothed signals have the same
#    median signal across sources.
#    For more details, see below.
#  }
#  \item{tags}{(Optional) Sets the tags for the output data sets.}
#  \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   The multi-source normalization method is by nature a single-sample method,
#   that is, it normalizes arrays for one sample at the time and independently
#   of all other samples/arrays.
#
#   However, the current implementation is such that it first generates
#   smoothed data for \emph{all} samples/arrays.  Then, it normalizes the
#   sample one by one.
# }
#
# \section{Different preprocessing methods normalize ChrX \& ChrY differently}{
#    Some preprocessing methods estimate copy numbers on sex chromosomes
#    differently from the autosomal chromosomes.  The way this is done may
#    vary from method to method and we cannot assume anything about what
#    approach is.  This is the main reason why the estimation of the
#    normalization  function is by default based on signals from autosomal
#    chromosomes only; this protects the estimate of the function from
#    being biased by specially estimated sex-chromosome signals.
#    Note that the normalization function is still applied to all chromosomes.
#
#    This means that if the transformation applied by a particular
#    preprocessing method is not the same for the sex chromosomes as the
#    autosomal chromosomes, the normalization applied on the sex
#    chromosomes is not optimal one.  This is why multi-source
#    normalization sometimes fails to bring sex-chromosome signals
#    to the same scale across sources.  Unfortunately, there is no
#    automatic way to handle this.
#    The only way would be to fit a specific normalization function to each
#    of the sex chromosomes, but that would require that there exist
#    copy-number abberations on those chromosomes, which could be a too
#    strong assumption.
#
#    A more conservative approach is to normalize the signals such that
#    afterward the median of the smoothed copy-number levels are the same
#    across sources for any particular chromosome.
#    This is done by setting argument \code{align="byChromosome"}.
# }
#
# \references{
#   [1] H. Bengtsson, A. Ray, P. Spellman & T.P. Speed,
#       \emph{A single-sample method for normalizing and combining
#         full-resolution copy numbers from multiple platforms,
#         labs and analysis methods},
#       Bioinformatics 2009. \cr
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("MultiSourceCopyNumberNormalization", function(dsList=NULL, fitUgp=NULL, subsetToFit=NULL, targetDimension=1, align=c("byChromosome", "none"), tags="*", ...) {
  if (!is.null(dsList)) {
    # aroma.light::fitPrincipalCurve()
    .requirePkg("aroma.light", quietly=TRUE);

    # Arguments 'dsList':
    if (is.list(dsList)) {
      K <- length(dsList);

      className <- "AromaUnitTotalCnBinarySet";
      for (kk in seq_len(K)) {
        ds <- dsList[[kk]];
        ds <- Arguments$getInstanceOf(ds, className, .name="dsList");
      }
      if (length(dsList) < 2L) {
        throw("Argument 'dsList' must contain more than one ",
                                                         className, ": ", K);
      }
    } else {
      throw("Argument 'dsList' is not a list: ", class(dsList)[1L]);
    }

    # Arguments 'fitUgp':
    fitUgp <- Arguments$getInstanceOf(fitUgp, "AromaUgpFile");

    # Argument 'subsetToFit':
    if (is.null(subsetToFit)) {
    } else if (is.character(subsetToFit)) {
      throw("Yet not implemented: Argument 'subsetToFit' is of type character.");
    } else {
      subsetToFit <- Arguments$getIndices(subsetToFit, max=nbrOfUnits(fitUgp));
    }

    # Argument 'align'
    align <- match.arg(align);

    # Argument 'targetDimension'
    targetDimension <- Arguments$getIndex(targetDimension, max=K);
  }

  # Arguments '...':
  args <- list(...);
  if (length(args) > 0L) {
    argsStr <- paste(names(args), collapse=", ");
    throw("Unknown arguments: ", argsStr);
  }

  extend(Object(), c("MultiSourceCopyNumberNormalization", uses("ParametersInterface")),
    .tags = tags,
    .dsList = dsList,
    .fitUgp = fitUgp,
    .subsetToFit = subsetToFit,
    .align = align,
    .targetDimension = targetDimension,
    "cached:.dsSmoothList" = NULL
  )
})


setMethodS3("as.character", "MultiSourceCopyNumberNormalization", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1L]);

  # Tags:
  tags <- getTags(this, collapse=", ");
  s <- c(s, sprintf("Tags: %s", tags));

  # Data sets:
  dsList <- getInputDataSets(this);
  s <- c(s, sprintf("Data sets (%d):", length(dsList)));
  for (kk in seq_along(dsList)) {
    ds <- dsList[[kk]];
    s <- c(s, as.character(ds));
  }

  # All common array names:
  names <- getAllNames(this);
  n <- length(names);
  s <- c(s, sprintf("Number of common array names: %d", n));
  s <- c(s, sprintf("Names: %s [%d]", hpaste(names), n));

  # Parameters:
  s <- c(s, sprintf("Parameters: %s", getParametersAsString(this)));
  GenericSummary(s);
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getInputDataSets
#
# @title "Gets the list of data sets to be normalized"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @list.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getInputDataSets", "MultiSourceCopyNumberNormalization", function(this, ...) {
  this$.dsList;
})



setMethodS3("nbrOfDataSets", "MultiSourceCopyNumberNormalization", function(this, ...) {
  length(getInputDataSets(this));
});


setMethodS3("getAsteriskTags", "MultiSourceCopyNumberNormalization", function(this, ...) {
  tags <- "mscn";

  # Align-by-chromosome tag?
  align <- this$.align;
  if (align != "none") {
    tags <- c(tags, align);
  }

  tags <- paste(tags, collapse=",");
  tags;
})

setMethodS3("getTags", "MultiSourceCopyNumberNormalization", function(this, collapse=NULL, ...) {
  tags <- this$.tags;

  # Split tags
  tags <- unlist(strsplit(tags, split=","));

  # Asterisk tags
  tags[tags == "*"] <- getAsteriskTags(this);

  # Split tags
  tags <- unlist(strsplit(tags, split=","));

  # Collapse?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
})

setMethodS3("getOutputPaths", "MultiSourceCopyNumberNormalization", function(this, ...) {
  dsList <- getInputDataSets(this);
  tags <- getTags(this);

  paths <- lapply(dsList, FUN=function(ds) {
    path <- getPath(ds);
    path <- getParent(path, 2L);
    rootPath <- basename(path);
    path <- getParent(path);
    rootPath <- "cnData";
    path <- Arguments$getWritablePath(rootPath);

    fullname <- getFullName(ds);
    fullname <- paste(c(fullname, tags), collapse=",");
    chipType <- getChipType(ds);
    file.path(path, fullname, chipType);
  });
  paths <- unlist(paths, use.names=FALSE);
  paths;
}, protected=TRUE)


setMethodS3("getOutputDataSets", "MultiSourceCopyNumberNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving list of output data sets");

  dsList <- getInputDataSets(this);
  paths <- getOutputPaths(this);
  dsOutList <- list();
  for (kk in seq_along(dsList)) {
    ds <- dsList[[kk]];
    verbose && enter(verbose, sprintf("Data set %d ('%s') of %d",
                                    kk, getFullName(ds), length(dsList)));

    path <- paths[[kk]];
    if (isDirectory(path)) {
      verbose && enter(verbose, "Scanning directory for matching data files");
      verbose && cat(verbose, "Path: ", path);

      dsOut <- byPath(ds, path=path, ..., verbose=less(verbose, 10));

      verbose && enter(verbose, "Keeping output data files matching input data files");
      # Identify output data files that match the input data files
      fullnames <- getFullNames(ds);
      df <- getFile(ds, 1L);
      translator <- getFullNameTranslator(df);
      setFullNamesTranslator(dsOut, translator);
      fullnamesOut <- getFullNames(dsOut);
      idxs <- match(fullnames, fullnamesOut);
      verbose && str(verbose, idxs);
      if (anyMissing(idxs)) {
        throw("Should not happen.");
      }
      verbose && cat(verbose, "Number of files dropped: ", length(dsOut) - length(idxs));
      verbose && cat(verbose, "Number of files kept: ", length(idxs));
      dsOut <- extract(dsOut, idxs);
      verbose && exit(verbose);

      verbose && exit(verbose);
    } else {
      dsOut <- NA;
    }
    dsOutList[[kk]] <- dsOut;

    verbose && exit(verbose);
  } # for (kk ...)

  verbose && exit(verbose);

  dsOutList;
})



###########################################################################/**
# @RdocMethod getFitAromaUgpFile
#
# @title "Gets the UGP file specifying the common set of loci to normalize at"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "aroma.core::AromaUgpFile".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFitAromaUgpFile", "MultiSourceCopyNumberNormalization", function(this, ...) {
  this$.fitUgp;
}, protected=TRUE)




###########################################################################/**
# @RdocMethod getAllNames
#
# @title "Gets the names of all unique samples across all sources"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Passed to \code{getNames(...)} of each data set.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAllNames", "MultiSourceCopyNumberNormalization", function(this, ...) {
  # Identify all array names across all sources
  dsList <- getInputDataSets(this);
  allNames <- lapply(dsList, getNames, ...);
  allNames <- unlist(allNames, use.names=FALSE);
  allNames <- unique(allNames);
  allNames <- sort(allNames);
  allNames;
})



###########################################################################/**
# @RdocMethod extractTupleOfDataFiles
#
# @title "Gets a list of data files for a particular name across several data sets"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{A @character string specifying the sample name of interest.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @list of K @see "aroma.core::AromaUnitTotalCnBinarySet":s.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractTupleOfDataFiles", "MultiSourceCopyNumberNormalization", function(this, dsList, name, ..., na.rm=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'dsList':
  if (is.list(dsList)) {
    className <- "AromaUnitTotalCnBinarySet";
    for (kk in seq_along(dsList)) {
      ds <- dsList[[kk]];
      ds <- Arguments$getInstanceOf(ds, className, .name="dsList");
    }
    if (length(dsList) < 2L) {
      throw("Argument 'dsList' must contain more than one ", className,
                                                     ": ", length(dsList));
    }
  } else {
    throw("Argument 'dsList' is not a list: ", class(dsList)[1L]);
  }

  # Argument 'name':
  name <- Arguments$getCharacter(name);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Getting list tuple of data files for one sample");
  verbose && cat(verbose, "Sample name: ", name);

  dfList <- lapply(dsList, function(ds) {
    idx <- indexOf(ds, name);
    df <- NA;
    if (!is.na(idx)) {
      if (length(idx) > 1L) {
        throw("Multiple occurances identified for this sample: ",
                           getName(ds), " => ", paste(idx, collapse=", "));
      }
      df <- getFile(ds, idx);
    }
    df;
  });

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Filter out missing data files in order to identify the set of files
  # to fit the model on
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (na.rm) {
    keep <- sapply(dfList, FUN=function(df) !identical(df, NA));
    dfList <- dfList[keep];
  }

  verbose && cat(verbose, "Number of arrays: ", length(dfList));

  verbose && exit(verbose);

  dfList;
}, protected=TRUE)





###########################################################################/**
# @RdocMethod getSmoothedDataSets
#
# @title "Gets the data sets smoothed toward the UGP file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @list of K @see "aroma.core::AromaUnitTotalCnBinarySet":s.
# }
#
# \details{
#   This method smooth each data set, each array, and each chromosome
#   toward the target (smoothing) UGP file independently of everything else.
#
#   The resulting data sets are stored in a separate location where they
#   will be located automatically in future sessions.
# }
#
# @author
#
# \seealso{
#   @seemethod "getFitAromaUgpFile".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getSmoothedDataSets", "MultiSourceCopyNumberNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  dsSmoothList <- this$.dsSmoothList;
  if (is.null(dsSmoothList)) {
    verbose && enter(verbose, "Smoothing all data sets to the same set of loci");
    dsList <- getInputDataSets(this);
    verbose && cat(verbose, "Number of data sets: ", length(dsList));

    targetUgp <- getFitAromaUgpFile(this);
    verbose && print(verbose, targetUgp);

    kernel <- "gaussian";
    sd <- 50e3;
    verbose && printf(verbose, "Kernel: %s\n", kernel);
    verbose && printf(verbose, "Bandwidth (sd): %.2f\n", sd);

    dsSmoothList <- list();
    for (kk in seq_along(dsList)) {
      ds <- dsList[[kk]];
      verbose && enter(verbose, sprintf("Data set %d ('%s') of %d",
                                         kk, getFullName(ds), length(dsList)));
      sm <- TotalCnKernelSmoothing(ds, targetUgp=targetUgp,
                                       kernel=kernel, bandwidth=sd);
      verbose && print(verbose, sm);
      dsSmoothList[[kk]] <- process(sm, verbose=less(verbose, 1));
      verbose && exit(verbose);
    }
    names(dsSmoothList) <- names(dsList);

    # Cache in memory
    this$.dsSmoothList <- dsSmoothList;

    verbose && exit(verbose);
  }


  dsSmoothList;
}, protected=TRUE)







###########################################################################/**
# @RdocMethod getSubsetToFit
#
# @title "Gets subset of (smoothing) units for fitting the model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @integer @vector of unit indices.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getSubsetToFit", "MultiSourceCopyNumberNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  units <- this$.subsetToFit;
  if (is.null(units)) {
    verbose && enter(verbose, "Identify subset of (smoothed) units for fitting the model");

    ugp <- getFitAromaUgpFile(this);
    verbose && print(verbose, ugp);

    verbose && enter(verbose, "Querying UGP for units on chromosomes of interest");
    chromosomes <- 1:22;
    verbose && cat(verbose, "Chromosomes to fit: ",
                                               seqToHumanReadable(chromosomes));
    units <- sapply(chromosomes, FUN=function(cc) {
      getUnitsOnChromosome(ugp, cc);
    });
    units <- unlist(units, use.names=FALSE);
    units <- unique(units);
    units <- sort(units);
    verbose && str(verbose, units);
    verbose && exit(verbose);

    this$.subsetToFit <- units;

    verbose && exit(verbose);
  }


  units;
}, protected=TRUE)



setMethodS3("getParameters", "MultiSourceCopyNumberNormalization", function(this, ...) {
  params <- NextMethod("getParameters");

  params$subsetToFit <- getSubsetToFit(this, ...);
  params$fitUgp <- getFitAromaUgpFile(this, ...);
  params$align <- this$.align;
  params$targetDimension <- this$.targetDimension;
  params$pcBandwidth <- this$.pcBandwidth;

  params;
}, protected=TRUE)


setMethodS3("getPrincipalCurveEstimator", "MultiSourceCopyNumberNormalization", function(this, ...) {
  # aroma.light::fitPrincipalCurve()
  .requirePkg("aroma.light", quietly=TRUE);

  params <- getParameters(this);
  df <- params$pcBandwidth;
  if (is.null(df)) {
    df <- 5;
  }
  df <- Arguments$getDouble(df);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  smoother <- function(lambda, xj, ...) {
    o <- order(lambda);
    lambda <- lambda[o];
    xj <- xj[o];
    fit <- smooth.spline(lambda, xj, ..., df=df, keep.data=FALSE);
    predict(fit, x=lambda)$y;
  }

  robustSmoother <- function(lambda, xj, ...) {
    o <- order(lambda);
    lambda <- lambda[o];
    xj <- xj[o];
    fit <- .robustSmoothSpline(lambda, xj, ..., df=df);
    predict(fit, x=lambda)$y;
  }

  # Create principal curve estimator
  fcn <- function(Y, ...) {
    .fitPrincipalCurve(Y, smoother=smoother, ...);
  }
  attr(fcn, "smoother") <- smoother;
  attr(fcn, "df") <- df;

  fcn;
}, protected=TRUE);



###########################################################################/**
# @RdocMethod fitOne
#
# @title "Fits the multi-source model for one sample"
#
# \description{
#  @get "title".
#  The model is fitted on the subset of units returned
#  by @seemethod "getSubsetToFit".
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{A @character string specifying the sample name of interest.}
#   \item{...}{Not used.}
#   \item{force}{If @FALSE, cached model fits are returned, otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @list of transforms.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fitOne", "MultiSourceCopyNumberNormalization", function(this, dfList, ..., force=FALSE, .retData=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'force':
  force <- Arguments$getLogical(force);



  verbose && enter(verbose, "Fitting one sample across multiple sources");

  if (is.character(dfList)) {
    verbose && enter(verbose, "Extracting list of input data files");
    name <- dfList;
    verbose && cat(verbose, "Sample name: ", name);
    dsList <- getInputDataSets(this);
    dfList <- extractTupleOfDataFiles(this, dsList=dsList, name=name,
                                                verbose=less(verbose, 1));
    verbose && print(verbose, dfList);
    verbose && exit(verbose);
  }

  nbrOfArrays <- length(dfList);
  verbose && cat(verbose, "Number of arrays: ", nbrOfArrays);

  # Get name of the sample from the tuple of input arrays
  # (We do it this way so that we at some stage can process() one sample
  #  at the time without first smoothing all samples. /HB 2008-08-18)
  df <- dfList[[1]];
  name <- getName(df);
  verbose && cat(verbose, "Sample name: ", name);
  # Not needed anymore
  dfList <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get model parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- getParameters(this, verbose=less(verbose, 1));
  verbose && str(verbose, params);
  subsetToFit <- params$subsetToFit;
  align <- params$align;
  targetDimension <- params$targetDimension;
  pcEstimator <- getPrincipalCurveEstimator(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify list of data files to fit model to
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Smooth data towards target UGP, which specifies the common set of loci
  dsSmooth <- getSmoothedDataSets(this, verbose=less(verbose, 1));
  dfSList <- extractTupleOfDataFiles(this, dsList=dsSmooth, name=name,
                                                 verbose=less(verbose, 1));
  # Not needed anymore
  dsSmooth <- NULL;
  verbose && str(verbose, dfSList);

  # Identify and exlude missing data sets
  keep <- sapply(dfSList, FUN=function(df) !identical(df, NA));
  keep <- which(keep);
  dfSList <- dfSList[keep];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already fitted?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fullnames <- sapply(dfSList, getFullName);
  fullnames <- unname(fullnames);

  chipTypes <- sapply(dfSList, getChipType);
  chipTypes <- unname(chipTypes);

  checkSums <- sapply(dfSList, getChecksum);
  checkSums <- unname(checkSums);

  df <- params$pcBandwidth;

  key <- list(method="fitOne", class="MultiSourceCopyNumberNormalization",
           fullnames=fullnames, chipTypes=chipTypes, checkSums=checkSums,
           subsetToFit=subsetToFit, align=align, df=df,
           .retData=.retData, version="2010-01-14");
  dirs <- c("aroma.cn", "MultiSourceCopyNumberNormalization");
  if (!force) {
    fit <- loadCache(key=key, dirs=dirs);
    if (!is.null(fit)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(fit);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract smoothed data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting data");
  verbose && cat(verbose, "Subset of units used for fitting:");
  verbose && str(verbose, subsetToFit);
  # Extracting data for sample to be normalized
  Y <- lapply(dfSList, FUN=function(df) {
    extractMatrix(df, rows=subsetToFit, column=1, drop=TRUE);
  });

  # Not needed anymore
  subsetToFit <- NULL;

  Y <- as.data.frame(Y);
  colnames(Y) <- NULL;
  Y <- as.matrix(Y);
  dimnames(Y) <- NULL;
  dim <- dim(Y);
  verbose && str(verbose, Y);
  verbose && summary(verbose, Y);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit principal curve to smoothed data (Y[,1], Y[,2], ..., Y[,K])
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Fitting across-source normalization function");
  verbose && cat(verbose, "Estimator for principal curves:");
  verbose && str(verbose, pcEstimator);
  t <- system.time({
    fit <- pcEstimator(Y);
  });
  verbose && cat(verbose, "Fitting time:");
  verbose && print(verbose, t);

  # Flip direction of the curve ('lambda')?
  rho <- cor(fit$lambda, Y[,1], use="complete.obs");
  flip <- (rho < 0);
  if (flip) {
    fit$lambda <- max(fit$lambda, na.rm=TRUE) - fit$lambda;
    verbose && cat(verbose, "Direction of fitted curve ('lambda') was flipped such that it increases with the signal.");
  }

  verbose && printf(verbose, "Processing time: %.1f seconds\n",
                                                          as.double(t[3L]));

  if (.retData) {
    fit$Y <- Y;
  }
  # Not needed anymore
  Y <- NULL;

  # Sanity check
  if (!identical(dim(fit$s), dim)) {
    throw("Internal error: The fitted data has a different dimension that the input data: ",
                         paste(dim(fit$s), collapse="x"), " != ", paste(dim, collapse="x"));
  }
  verbose && str(verbose, fit);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Standardize the channels to a target channel?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  targetChannel <- NULL;
  if (!is.null(targetChannel)) {
##     for (kk in seq_len(dim[2])) {
##       if (kk == targetChannel) {
##         targetTransform <- function(x, ...) x;
##       } else {
##         targetTransform <- makeSmoothSplinePredict(Yn[,kk], Yn[,targetChannel]);
##       }
##     } # for (kk ...)
  }

#  class(fit) <- c("MultiSourceCopyNumberNormalizationFit", class(fit));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Shift each chromosome?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.element(align, c("byChromosome"))) {
    verbose && enter(verbose, "Calculating shift for each chromosome");
    verbose && cat(verbose, "align=", align);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Grouping units by chromosome
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ugpS <- getAromaUgpFile(dfSList[[1L]]);
    chromosomes <- getChromosomes(ugpS);
    verbose && cat(verbose, "Chromosomes: ", seqToHumanReadable(chromosomes));

    verbose && enter(verbose, "Grouping units by chromosome");
    values <- ugpS[,1L,drop=TRUE];
    unitsS <- list();
    for (chr in chromosomes) {
      chrStr <- sprintf("Chr%02d", chr);
      unitsS[[chrStr]] <- which(values == chr);
    }
    # Not needed anymore
    values <- NULL;
#    verbose && str(verbose, unitsS);
    # Dropping chromosomes with too few units
    ns <- sapply(unitsS, FUN=length);
    unitsS <- unitsS[ns > 5L];
    verbose && str(verbose, unitsS);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculating means of each chromosome in each source
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Allocating matrix for smooth data");
    dfS <- dfSList[[1L]];
    naValue <- as.double(NA);
    YSN <- matrix(naValue, nrow=nbrOfUnits(dfS), ncol=nbrOfArrays);
    verbose && cat(verbose, "RAM: ", objectSize(YSN), " bytes");
    verbose && exit(verbose);

    verbose && enter(verbose, "Loading and backtransforming *smoothed* data");
    for (kk in seq_len(nbrOfArrays)) {
      dfS <- dfSList[[kk]];
      verbose && enter(verbose, sprintf("Source #%d ('%s') of %d", kk,
                                        getFullName(dfS), nbrOfArrays));

      verbose && enter(verbose, "Loading smoothed data");
      yS <- extractMatrix(dfS, column=1, drop=TRUE);
      verbose && str(verbose, yS);
      verbose && exit(verbose);

      verbose && enter(verbose, "Backtransforming smoothed data");
      ySN <- .backtransformPrincipalCurve(yS, fit=fit, dimensions=kk,
                                      targetDimension=targetDimension);
      ySN <- ySN[,1L,drop=TRUE];
      verbose && str(verbose, ySN);
      # Not needed anymore
      yS <- NULL;
      verbose && exit(verbose);

      # Storing
      YSN[,kk] <- ySN;
      # Not needed anymore
      ySN <- NULL;

      verbose && exit(verbose);
    } # for (kk ...)
    verbose && summary(verbose, YSN);
    verbose && str(verbose, YSN);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating shifts chromosome by chromosome");

    # Allocate matrices to hold all mean and shift values
    nbrOfChromosomes <- length(unitsS);
    naValue <- as.double(NA);
    mus <- matrix(naValue, nrow=nbrOfChromosomes, ncol=nbrOfArrays);
    rownames(mus) <- names(unitsS);
    dmus <- mus;

    for (chr in seq_len(nbrOfChromosomes)) {
      chrStr <- sprintf("Chr%02d", chr);
      verbose && enter(verbose, sprintf("Chromosome #%d of %d",
                                                     chr, nbrOfChromosomes));
      # Get the units
      unitsCC <- unitsS[[chrStr]];

      verbose && enter(verbose, "Extracting backtransformed *smoothed* data");
      yList <- list();
      for (kk in seq_len(nbrOfArrays)) {
        yList[[kk]] <- YSN[unitsCC,kk,drop=TRUE];
      } # for (kk ...)
      verbose && str(verbose, yList);
      verbose && exit(verbose);

      verbose && enter(verbose, "Estimating averages and shifts toward targetDimension");
      verbose && cat(verbose, "Target dimension: ", targetDimension);
      # Estimate averages and shifts toward targetDimension
      yNList <- .normalizeDifferencesToAverage(yList, baseline=targetDimension);
      alignFit <- attr(yNList, "fit");
      verbose && str(verbose, alignFit);
      verbose && exit(verbose);

      mus[chrStr,] <- alignFit$mus;
      dmus[chrStr,] <- alignFit$deltas;

      # Not needed anymore
      alignFit <- yList <- yNList <- NULL;
      verbose && exit(verbose);
    } # for (chr ...)
    verbose && exit(verbose);

    # Not needed anymore
    YSN <- NULL;

    verbose && cat(verbose, "Overall averages:");
    verbose && print(verbose, mus);
    verbose && cat(verbose, "Overall shifts:");
    verbose && print(verbose, dmus);
    verbose && cat(verbose, "Target dimension: ", targetDimension);
    verbose && exit(verbose);

    fit$alignParams <- list(
      dmus=dmus,
      mus=mus
    );

    verbose && exit(verbose);
  } # if (align ...)

  verbose && str(verbose, fit);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  saveCache(key=key, dirs=dirs, fit);

  fit;
}, protected=TRUE)  # fitOne()






setMethodS3("normalizeOne", "MultiSourceCopyNumberNormalization", function(this, dfList, fit, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dfList':

  # Argument 'fit':


  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'force':
  force <- Arguments$getLogical(force);



  verbose && enter(verbose, "Normalize one sample across multiple sources");

  # Get name of the sample from the tuple of input arrays
  # (We do it this way so that we at some stage can process() one sample
  #  at the time without first smoothing all samples. /HB 2008-08-18)
  df <- dfList[[1]];
  name <- getName(df);
  verbose && cat(verbose, "Sample name: ", name);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get model parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- getParameters(this, verbose=less(verbose, 1));
  verbose && str(verbose, params);
  subsetToUpdate <- params$subsetToUpdate;
  targetDimension <- params$targetDimension;
  align <- params$align;

  # Get (and create) the output paths
  outputPaths <- getOutputPaths(this);

  if (is.element(align, c("byChromosome"))) {
    verbose && enter(verbose, "Estimate alignment parameters");
    verbose && cat(verbose, "align=", align);

    verbose && enter(verbose, "Extracting align-by-chromosome parameters");
    alignParams <- fit$alignParams;
    verbose && str(verbose, alignParams);
    # Sanity check
    if (is.null(alignParams)) {
      throw("Internal error: No shift estimates found.");
    }

    dmus <- alignParams$dmus;
    verbose && print(verbose, dmus);
    # Sanity check
    if (is.null(dmus)) {
      throw("Internal error: No shift estimates found.");
    }
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # if (align ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalizing array by array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing source by source (array by array)");
  verbose && cat(verbose, "Units to be updated:");
  verbose && str(verbose, subsetToUpdate);

  nbrOfArrays <- length(dfList);
  dfNList <- vector("list", length=nbrOfArrays);
  for (kk in seq_len(nbrOfArrays)) {
    df <- dfList[[kk]];
    verbose && enter(verbose, sprintf("Source #%d ('%s') of %d", kk,
                                            getFullName(df), nbrOfArrays));

    outputPath <- outputPaths[[kk]];
    # Here we really should use the fullname /HB 2009-05-05
    filename <- getFilename(df);
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Already normalized.");
      dfN <- newInstance(df, pathname);
    } else {
      verbose && enter(verbose, "Normalizing");

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Reading data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Reading data");
      y <- extractMatrix(df, rows=subsetToUpdate, column=1, drop=TRUE);
      verbose && str(verbose, y);
      verbose && exit(verbose);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Normalizing data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Backtransforming data");
      yN <- .backtransformPrincipalCurve(y, fit=fit, dimensions=kk,
                                        targetDimension=targetDimension);
      yN <- yN[,1L,drop=TRUE];
      verbose && str(verbose, yN);
      verbose && exit(verbose);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Aligning signals chromosome by chromosome?
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (is.element(align, c("byChromosome"))) {
        verbose && enter(verbose, "Align genomic signals");
        verbose && cat(verbose, "align=", align);

        verbose && enter(verbose, "Aligning signals for each chromosome");
        ugp <- getAromaUgpFile(df);
        chromosomes <- getChromosomes(ugp);
        verbose && cat(verbose, "Chromosomes: ", seqToHumanReadable(chromosomes));

        verbose && enter(verbose, "Grouping units by chromosome");
        values <- ugp[subsetToUpdate,1L,drop=TRUE];
        # Sanity check
        stopifnot(length(values) == length(yN));

        listOfUnits <- list();
        for (chr in chromosomes) {
          chrStr <- sprintf("Chr%02d", chr);
          subset <- which(values == chr);
          listOfUnits[[chrStr]] <- subset;
        }
        # Not needed anymore
        values <- NULL;
        verbose && str(verbose, listOfUnits);

        # Dropping chromosomes with too few units
        ns <- sapply(listOfUnits, FUN=length);
        listOfUnits <- listOfUnits[ns > 5L];
        verbose && str(verbose, listOfUnits);
        verbose && exit(verbose);

        # Dropping chromosomes for which there is no shift estimate
        idxs <- match(names(listOfUnits), rownames(dmus));
        if (anyMissing(idxs)) {
          verbose && cat(verbose, "Shift estimates are not available for some chromosomes, which are skipped:");
          verbose && print(verbose, names(listOfUnits[!is.finite(idxs)]));
          listOfUnits <- listOfUnits[is.finite(idxs)];
        }

        # Aligning mean signals chromosome by chromosome
        for (chrStr in names(listOfUnits)) {
          subset <- listOfUnits[[chrStr]];
          dmu <- dmus[chrStr,kk];
          yN[subset] <- yN[subset] - dmu;
        } # for (chrStr ...)

        # Not needed anymore
        listOfUnits <- NULL;

        verbose && str(verbose, yN);

        verbose && exit(verbose);

        verbose && exit(verbose);
      } # if (align ...)


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Writing normalized data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Storing normalized data");
      verbose && cat(verbose, "Output pathname: ", pathname);

      verbose && enter(verbose, "Create output file");
      pathnameT <- sprintf("%s.tmp", pathname);
      verbose && cat(verbose, "Temporary pathname: ", pathnameT);
      pathnameT <- Arguments$getWritablePathname(pathnameT, mustNotExist=TRUE);

      file.copy(getPathname(df), pathnameT);
      dfN <- newInstance(df, pathnameT);
      verbose && print(verbose, dfN);
      verbose && exit(verbose);

      verbose && enter(verbose, "Writing data");
      if (is.null(subsetToUpdate)) {
        dfN[,1L] <- yN;
      } else {
        dfN[subsetToUpdate,1L] <- yN;
      }
      # Not needed anymore
      yN <- NULL;
      verbose && exit(verbose);

      verbose && enter(verbose, "Updating file footer");
      footer <- readFooter(dfN);
      srcFile <- df;
      footer$srcFile <- list(
        filename = getFilename(srcFile),
        filesize = getFileSize(srcFile),
        checksum = getChecksum(srcFile)
      );
      pkg <- aroma.cn;
      footer$createdBy <- list(
        class=class(this)[1],
        package = getName(pkg),
        version = getVersion(pkg)
      );
      footer$createdOn <- format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE);
      writeFooter(dfN, footer);
      verbose && exit(verbose);

      verbose && enter(verbose, "Renaming temporary filename");
      file.rename(pathnameT, pathname);
      if (isFile(pathnameT) || !isFile(pathname)) {
        throw("Failed to rename temporary file: ",
                        pathnameT, " -> ", pathname);
      }
      # Not needed anymore
      pathnameT <- NULL;
      verbose && exit(verbose);
      dfN <- newInstance(df, pathname);
      # Not needed anymore
      pathname <- NULL;

      verbose && exit(verbose);


      verbose && exit(verbose);
    }

    verbose && print(verbose, dfN);
    dfNList[[kk]] <- dfN;
    # Not needed anymore
    dfN <- NULL;

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && print(verbose, dfNList);

  # Not needed anymore
  subsetToUpdate <- NULL;

  verbose && exit(verbose);

  # Return normalized arrays
  invisible(dfNList);
}, protected=TRUE)  # normalizeOne()






###########################################################################/**
# @RdocMethod process
#
# @title "Normalizes all samples"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @list of K @see "aroma.core::AromaUnitTotalCnBinarySet":s.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "MultiSourceCopyNumberNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit normalization functions
  #
  # This is a multi-source (same sample across sources) whole-genome method.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Multi-source normalize all samples");
  allNames <- getAllNames(this);
  nbrOfSamples <- length(allNames);
  verbose && cat(verbose, "Number of unique samples in all sets: ",
                                                               nbrOfSamples);
  verbose && str(verbose, allNames);

  # Get the input data sets
  dsList <- getInputDataSets(this);

  # Get (and create) the output paths
  outputPaths <- getOutputPaths(this);


  verbose && enter(verbose, "Processing each array");
  for (kk in seq_len(nbrOfSamples)) {
    name <- allNames[kk];
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d",
                                                    kk, name, nbrOfSamples));


    verbose && enter(verbose, "Identifying source data files");
    dfList <- extractTupleOfDataFiles(this, dsList=dsList, name=name,
                                                   verbose=less(verbose, 1));
    verbose && print(verbose, dfList);
    verbose && exit(verbose);


    verbose && enter(verbose, "Check if all arrays are already normalized");
    isDone <- TRUE;
    for (jj in seq_along(dfList)) {
      df <- dfList[[jj]];
      outputPath <- outputPaths[[jj]];
      filename <- getFilename(df);
      pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);
      isDone <- isDone && isFile(pathname);
      if (!isDone)
        break;
    }
    verbose && cat(verbose, "Is done: ", isDone);
    verbose && exit(verbose);

    if (!force && isDone) {
      verbose && cat(verbose, "Normalized data files already exist");
    } else {
      verbose && enter(verbose, "Fitting model");
      fit <- fitOne(this, dfList=dfList, ..., force=force,
                                           verbose=less(verbose, 1))
      verbose && str(verbose, fit);
      verbose && exit(verbose);


      verbose && enter(verbose, "Normalizing");
      dfNList <- normalizeOne(this, dfList=dfList, fit=fit, ...,
                             force=force, verbose=less(verbose, 1));
      # Not needed anymore
      fit <- NULL;
      verbose && print(verbose, dfNList);

      # Sanity check
      if (length(dfNList) != length(dfList)) {
        throw("The number of normalized arrays does not match the number of source arrays: ", length(dfNList), " != ", length(dfList));
      }

      verbose && exit(verbose);
      # Not needed anymore
      dfNList <- NULL;
    }

    # Not needed anymore
    dfList <- NULL;

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  # Not needed anymore
  dsList <- NULL;

  outputDataSets <- getOutputDataSets(this, force=TRUE, verbose=less(verbose, 1));

  verbose && exit(verbose);

  invisible(outputDataSets);
})


###########################################################################
# HISTORY:
# 2013-01-07
# o BUG FIX: getPrincipalCurveEstimator() used non-existing 'verbose'.
# 2012-11-21
# o Now class utilizes the new ParametersInterface.
# 2012-11-13
# o CLEANUP/FIX: Used "cache:" field modified instead of "cached:".
#   After correction, all clearCache() methods could be dropped.
# 2012-04-16
# o MultiSourceCopyNumberNormalization() now explicitly requires the
#   'aroma.light' package, instead of assuming it is loaded.
# 2010-04-04
# o Added citation for MSCN to the Rdocs.
# 2010-01-14
# o Added protected getPrincipalCurveEstimator() for the
#   MultiSourceCopyNumberNormalization class.  This is done in order to
#   one day support custom principal-curve estimators.
# 2009-09-30
# o Now the alignment is done using normalizeDifferencesToAverage(),
#   which is robust against outliers and waviness etc.  The previous
#   method which normalized towards the same overall median is dropped.
# o Renamed argument 'alignByChromosome' to "align" in order to allow for
#   more types of aligned.
# o BUG FIX: getTags() of MultiSourceCopyNumberNormalization would return
#   all asterisk tags as merged, e.g. c("mscn,align", "tagA", "tagB").
# 2009-05-17
# o Now the constructor of MultiSourceCopyNumberNormalization asserts that
#   there are no stray arguments.
# 2009-05-06
# o Now the 'alignByChromosome' is corrected using median estimates.
# 2009-05-05
# o Now getOutputDataSets() of  MultiSourceCopyNumberNormalization only
#   returns output data files with a matching fullname in the input set.
# o Added a clearCache() to MultiSourceCopyNumberNormalization.
# 2009-05-03
# o Added argument 'alignByChromosome'.  If used, the signals are shifted
#   per chromosome such that the mean of the normalized smoothed signals
#   is the same for all sources.
# o Added getAsteriskTags() and getTags().
# o Now normalized data is first written to a temporary file, which is
#   then renamed.
# 2009-02-08
# o Fixed verbose output of gc(); used cat() instead of print().
# o Updated to make use of TotalCnKernelSmoothing().
# 2009-01-26
# o Adapted to new aroma.core::AromaUnitTotalCnBinary{Set|File} classes.
# 2008-10-08
# o Added argument 'targetDimension' to the constructor.
# o Now fitOne() makes sure the fitted curve has a "positive" direction.
# o Added argument 'subsetToFit' with some support, but still incomplete.
# 2008-10-07
# o Updated fitOne() and normalizeOne() to make use of the updated/new
#   fit- and backtransformPrincipalCurve() functions.
# 2008-08-18
# o Added normalizeOne() and process().
# o Added utility function extractTupleOfDataFiles().
# o Added Rdoc comments.
# 2008-07-04
# o Added as.character().
# o BUG FIX: getAllNames() did return duplicated names.
# 2008-06-24
# o Created first stub from existing "manual" scripts.
# 2008-05-27
# o Created "manual" script.
###########################################################################
