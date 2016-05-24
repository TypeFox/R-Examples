###########################################################################/**
# @set "class=RawGenomicSignals"
# @RdocMethod segmentByHaarSeg
#
# @title "Segment copy numbers using the HaarSeg method"
#
# \description{
#  @get "title" of the \pkg{HaarSeg} package.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments passed to the segmentation function.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the fit object.
# }
#
# \details{
#   Internally \code{haarSeg()} of the \pkg{HaarSeg} is used to segment
#   the signals.
#   This segmentation method support weighted segmentation.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("segmentByHaarSeg", "RawGenomicSignals", function(this, ..., cache=FALSE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
  pkgName <- "HaarSeg";
  # Assert that package is installed
  isPackageInstalled(pkgName) || throw("Package is not installed: ", pkgName);
  pkg <- packageDescription(pkgName);
  pkgVer <- pkg$Version;
  pkgDetails <- sprintf("%s v%s", pkgName, pkgVer);

  methodName <- "haarSeg";
  verbose && cat(verbose, "Method: ", methodName);
  verbose && cat(verbose, "Package: ", pkgDetails);

  # We need to load the method
  require(pkgName, character.only=TRUE) || throw("Package not loaded: ", pkgName);

  # Get the fit function for the segmentation method
  envir <- as.environment(sprintf("package:%s", pkgName));
  fitFcn <- get(methodName, mode="function", envir=envir);
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

  verbose && cat(verbose, "Sample name: ", sampleName);
  verbose && cat(verbose, "Chromosome: ", chromosome);
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (hasWeights) {
    # Verify that weights are supported (not yet)
    if (!is.element("W", names(formals))) {
      hasWeights <- FALSE;
      msg <- paste("Weights detected but ignored, because the available segmentation function ('", methodName, "()') does not support weights. Check with a more recent version of the package: ", pkgDetails, sep="");
      verbose && cat(verbose, "WARNING: ", msg);
      warning(msg);
    }
  } # if (hasWeights)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up arguments to pass to segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up method arguments");

  verbose && enter(verbose, "Setting up ", pkgName, " data structure");
  cnData <- data$y;
  verbose && str(verbose, cnData);
  verbose && exit(verbose);

  params <- list();
  if (hasWeights) {
    params$W <- data$w;
    verbose && cat(verbose, "Additional segmentation arguments:");
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

  args <- c(list(I=cnData), params);
  verbose && cat(verbose, "Final arguments:");
  verbose && str(verbose, args);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Now, check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Looking for cached results");
  key <- list(method="segmentByHaarSeg", class=class(this)[1],
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calling segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, sprintf("Calling %s() of %s", methodName, pkgName));
  # In case the method writes to stdout, we capture it
  stdout <- capture.output({
    # This works, but requires that one loads the package and that the
    # function is not masked in the search() path.
    t <- system.time({
      fit <- do.call(methodName, args);
    });
    attr(fit, "processingTime") <- t;
    attr(fit, "pkgDetails") <- pkgDetails;
  });

  verbose && cat(verbose, "Captured output that was sent to stdout:");
  stdout <- paste(stdout, collapse="\n");
  verbose && cat(verbose, stdout);

  verbose && cat(verbose, "Fitting time (in seconds):");
  verbose && print(verbose, t);

  verbose && cat(verbose, "Fitting time per 1000 loci (in seconds):");
  verbose && print(verbose, 1000*t/nbrOfLoci);

  verbose && cat(verbose, "Results object:");
  verbose && str(verbose, fit);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup a HaarSeg fit object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up return HaarSeg object");
  fit <- list(
    output = fit,
    data   = list(M=data$y, x=data$x, chromosome=chromosome)
  );
  class(fit) <- "HaarSeg";
  attr(fit, "processingTime") <- t;
  attr(fit, "pkgDetails") <- pkgDetails;
  verbose && exit(verbose);


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
# 2009-05-10
# o Created.
############################################################################
