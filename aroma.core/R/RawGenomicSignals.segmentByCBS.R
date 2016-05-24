###########################################################################/**
# @set "class=RawGenomicSignals"
# @RdocMethod segmentByCBS
#
# @title "Segment copy numbers using the CBS method"
#
# \description{
#  @get "title" of the \pkg{DNAcopy} package.
#  For more details on the Circular Binary Segmentation (CBS) method
#  see [1,2].
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments passed to the segmentation function.}
#   \item{seed}{An (optional) @integer specifying the random seed to be
#     set before calling the segmentation method.  The random seed is
#     set to its original state when exiting.  If @NULL, it is not set.}
#   \item{cache}{If @TRUE, results are cached to file, otherwise not.}
#   \item{force}{If @TRUE, cached results are ignored.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the fit object.
# }
#
# \details{
#   Internally @see "DNAcopy::segment" is used to segment the signals.
#   This segmentation method support weighted segmentation.
#
#   The "DNAcopy::segment" implementation of CBS uses approximation
#   through random sampling for some estimates.  Because of this,
#   repeated calls using the same signals may result in slightly
#   different results.
# }
#
# @examples "../incl/RawGenomicSignals.SEG.Rex"
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \references{
#  [1] A.B. Olshen, E.S. Venkatraman (aka Venkatraman E. Seshan),
#      R. Lucito and M. Wigler, \emph{Circular binary segmentation for
#      the analysis of array-based DNA copy number data},
#      Biostatistics, 2004.\cr
#  [2] E.S. Venkatraman and A.B. Olshen, \emph{A faster circular binary
#      segmentation algorithm for the analysis of array CGH data}.
#      Bioinformatics, 2007.\cr
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("segmentByCBS", "RawGenomicSignals", function(this, ..., seed=NULL, cache=FALSE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Segmenting");
  verbose && cat(verbose, "Chromosomes: ", hpaste(getChromosomes(this)));

  # This is a single-chromosome method. Assert that is the case.
  assertOneChromosome(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving the fit function");
  pkgName <- "DNAcopy";
  # Assert that package is installed
  isPackageInstalled(pkgName) || throw("Package is not installed: ", pkgName);
  pkg <- packageDescription(pkgName);
  pkgVer <- pkg$Version;
  pkgDetails <- sprintf("%s v%s", pkgName, pkgVer);

  methodName <- "segment";
  verbose && cat(verbose, "Method: ", methodName);
  verbose && cat(verbose, "Package: ", pkgDetails);

  # We need to load package
  require(pkgName, character.only=TRUE) || throw("Package not loaded: ", pkgName);

  # Get the fit function for the segmentation method
  fitFcn <- getExportedValue(pkgName, methodName);
  fitFcn <- getFromNamespace(methodName, pkgName);
  verbose && str(verbose, "Function: ", fitFcn);
  formals <- formals(fitFcn);
  verbose && cat(verbose, "Formals:");
  verbose && str(verbose, formals);
  verbose && exit(verbose);

  signatures <- list();
  signatures$fitFcn <- list(
    pkgName=pkgName,
    methodName=methodName,
    formals=formals,
    pkgDetails=pkgDetails
  );

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting data of interest");
  data <- extractDataForSegmentation(this, ..., verbose=less(verbose, 5));
  verbose && str(verbose, data);
  verbose && exit(verbose);

  sampleName <- attr(data, "sampleName");
  chromosome <- data$chromosome[1];
  nbrOfLoci <- nrow(data);
  hasWeights <- !is.null(data$w);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (hasWeights) {
    # Verify that weights are supported (from DNAcopy v1.18.0);
    if (!is.element("weights", names(formals))) {
      hasWeights <- FALSE;
      msg <- paste("Weights detected but ignored, because the available segmentation function ('", methodName, "()') does not support weights. Check with a more recent version of the package: ", pkgDetails, sep="");
      verbose && cat(verbose, msg);
      warning(msg);
    }
  } # if (hasWeights)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up arguments to pass to segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up method arguments");

  verbose && enter(verbose, "Setting up ", pkgName, " data structure");
  cnData <- DNAcopy::CNA(
    genomdat  = data$y,
    chrom     = data$chromosome,
    data.type = "logratio",
    maploc    = data$x,
    sampleid  = sampleName
  );
  verbose && str(verbose, cnData);
  names(cnData)[3] <- sampleName;
  verbose && str(verbose, cnData);
  verbose && exit(verbose);

  params <- list();
  if (hasWeights) {
    params$weights <- data$w;
    verbose && cat(verbose, "Segmentation parameters:");
    verbose && str(verbose, params);
  }

  userArgs <- list(...);
  if (length(userArgs) > 0) {
    verbose && cat(verbose, "User and segmentation arguments:");
    verbose && str(verbose, userArgs);
    # Assign/overwrite by user arguments
    for (ff in names(userArgs)) {
      params[[ff]] <- userArgs[[ff]];
    }
  }

  # Cleaning out unknown parameters?
  if (!any(names(formals) == "...")) {
    keep <- (names(params) %in% names(formals));
    params <- params[keep];
  }

  signatures$data <- cnData;
  signatures$params <- params;
  signatures$seed <- seed;

  args <- c(list(cnData), params, verbose=as.logical(verbose));
  verbose && cat(verbose, "Final arguments:");
  verbose && str(verbose, args);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Now, check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Looking for cached results");
  key <- list(method="segmentByCBS", class=class(this)[1],
                                                signatures=signatures);
  dirs <- c("aroma.cn", class(this)[1]);
  if (!force) {
    res <- loadCache(key, dirs=dirs);
    if (!is.null(res)) {
      verbose && cat(verbose, "Found cached results.");
      verbose && exit(verbose);
      return(res);
    }
  }
  verbose && exit(verbose);
  # Not needed anymore
  signatures <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set the random seed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(seed)) {
    verbose && enter(verbose, "Setting (temporary) random seed");
    oldRandomSeed <- NULL;
    if (exists(".Random.seed", mode="integer")) {
      oldRandomSeed <- get(".Random.seed", mode="integer");
    }
    on.exit({
      if (!is.null(oldRandomSeed)) {
        .Random.seed <<- oldRandomSeed;
      }
    }, add=TRUE);
    verbose && cat(verbose, "The random seed will be reset to its original state afterward.");
    verbose && cat(verbose, "Seed: ", seed);
    set.seed(seed);
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calling segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, sprintf("Calling %s() of %s", methodName, pkgName));
  # In case the method writes to stdout, we capture it
  # Note: DNAcopy::segment() *does* this.
  stdout <- capture.output({
    # Does not work, because some internal function of the fit function
    # may only be accessible from within the namespace
    # How to do this for DNAcopy::segment()? /HB
##    fit <- do.call(fitFcn, args);
    # This works, but requires that one loads the package and that the
    # function is not masked in the search() path.
    t <- system.time({
      fit <- do.call(methodName, args);
    });
    # Drop the 'call' (because it will be huge due to the do.call() call)
    fit$call <- NULL;
  });
  attr(fit, "processingTime") <- t;
  attr(fit, "pkgDetails") <- pkgDetails;
  attr(fit, "randomSeed") <- seed;

  verbose && cat(verbose, "Captured output that was sent to stdout:");
  stdout <- paste(stdout, collapse="\n");
  verbose && cat(verbose, stdout);

  verbose && cat(verbose, "Fitting time (in seconds):");
  verbose && print(verbose, t);

  verbose && cat(verbose, "Fitting time per 1000 loci (in seconds):");
  verbose && print(verbose, 1000*t/nbrOfLoci);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimating aroma parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Estimating aroma parameters");
  # Estimate the standard deviation
  sigma <- estimateStandardDeviation(this);

  # Estimate the standard *error* for each segment
  cnr <- extractCopyNumberRegions(fit);
  cnrData <- as.data.frame(cnr);
  regions <- as.matrix(cnrData[,c("start", "stop")]);
  nbrOfRegions <- nrow(regions);
  # Not needed anymore
  cnr <- cnrData <- NULL;
  x <- data$x;
  y <- data$y;
  naValue <- as.double(NA);
  sigmas <- rep(naValue, times=nbrOfRegions);
  for (kk in seq_len(nbrOfRegions)) {
    keep <- which(regions[kk,1] < x & x <= regions[kk,2]);
    t <- y[keep];
    t <- diff(t);
    t <- median(t, na.rm=TRUE)/sqrt(2);
    sigmas[kk] <- t;
  } # for (kk ...)
  # Not needed anymore
  x <- y <- t <- keep <- NULL;
  aromaEstimates <- list(
    stddevAll = sigma,
    stddevRegions = sigmas
  );
  attr(fit, "aromaEstimates") <- aromaEstimates;
  verbose && exit(verbose);

  verbose && cat(verbose, "Results object:");
  verbose && str(verbose, fit);

  verbose && exit(verbose);


  # Save cached results?
  if (cache) {
    saveCache(fit, key=key, dirs=dirs);
  }

  verbose && exit(verbose);

  fit;
}) # segmentByCBS()


############################################################################
# HISTORY:
# 2014-02-17
# o Now unknown user arguments are only dropped if the segmentation method
#   does not have a formal argument '...'.
# 2010-04-05
# o Added argument 'seed' to segmentByCBS().
# o Added Rdoc references to CBS papers.
# 2009-05-10
# o Created.
############################################################################
