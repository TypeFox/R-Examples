###########################################################################/**
# @set "class=RawGenomicSignals"
# @RdocMethod segmentByMPCBS
#
# @title "Segment copy numbers using the multi-platform CBS (mpCBS) method"
#
# \description{
#  @get "title" of the \pkg{mpcbs} package.
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
#   Internally \code{mpcbs.mbic()} of the \pkg{mpcbs} package is used
#   for segmenting the signals.
#   This segmentation method does not support weighted segmentation.
# }
#
# @examples "../incl/RawGenomicSignals.SEG,MP.Rex"
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("segmentByMPCBS", "RawGenomicSignals", function(this, ..., cache=FALSE, force=FALSE, verbose=FALSE) {
  # To please R CMD check in case 'mpcbs' is not available
  # This might be need in order for it to work on CRAN. /HB 2010-01-04
  if (!isPackageInstalled("mpcbs")) {
    merge.pos <- function(...) {};
    rm(list="merge.pos");
  }

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
  pkgName <- "mpcbs";
  # Assert that package is installed
  isPackageInstalled(pkgName) || throw("Package is not installed: ", pkgName);
  pkg <- packageDescription(pkgName);
  pkgVer <- pkg$Version;
  pkgDetails <- sprintf("%s v%s", pkgName, pkgVer);

  methodName <- "mpcbs.mbic";
  verbose && cat(verbose, "Method: ", methodName);
  verbose && cat(verbose, "Package: ", pkgDetails);

  # We need to load package
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

  # Sanity check
  if (is.null(data$id)) {
    data$id <- 1L;
    warning(sprintf("%s did not contain a 'id' locus field specifying source/platform. Assuming only one source/platform.", class(this)[1]));
  }

  sampleName <- attr(data, "sampleName");
  chromosome <- data$chromosome[1];
  nbrOfLoci <- nrow(data);
  hasWeights <- !is.null(data$w);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (hasWeights) {
    msg <- paste("Weights detected but ignored, because the available segmentation function ('", methodName, "()') does not support weights. Check with a more recent version of the package: ", pkgDetails, sep="");
    verbose && cat(verbose, msg);
    warning(msg);
    hasWeights <- FALSE;
  } # if (hasWeights)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up arguments to pass to segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up method arguments");

  verbose && enter(verbose, "Setting up ", pkgName, " data structure");
  ids <- sort(unique(data$id));
  dataById <- split(data, data$id);
  nbrOfIds <- length(dataById);

  y <- lapply(dataById, FUN=function(df) df$y);
  pos <- lapply(dataById, FUN=function(df) df$x);
  anchor <- merge.pos(pos);
  cnData <- list(y=y, pos=pos, anchor=anchor);
  verbose && str(verbose, cnData);
  verbose && exit(verbose);

  params <- list();

  # Override default arguments of mpcbs::mpcbs.mbic()
  # Do not plot
  params$plots <- FALSE;

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

  args <- c(cnData, params);
  verbose && cat(verbose, "Final arguments:");
  verbose && str(verbose, args);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Now, check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Looking for cached results");
  key <- list(method="segmentByMPCBS", class=class(this)[1],
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
  # Note: DNAcopy::segment() *does* this.
  stdout <- capture.output({
    # This works, but requires that one loads the package and that the
    # function is not masked in the search() path.
    t <- system.time({
      fit <- do.call(methodName, args);
    });
  });
  fit <- list(fit=fit);
  fit$chromosome <- chromosome;
  class(fit) <- "MPCBS";
  attr(fit, "processingTime") <- t;
  attr(fit, "pkgDetails") <- pkgDetails;

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
##  verbose && enter(verbose, "Estimating aroma parameters");
##  # Estimate the standard deviation (should be done per id)
##  sigma <- estimateStandardDeviation(this);
##
##  # Estimate the standard *error* for each segment
##  cnr <- extractCopyNumberRegions(fit);
##  cnrData <- as.data.frame(cnr);
##  regions <- as.matrix(cnrData[,c("start", "stop")]);
##  nbrOfRegions <- nrow(regions);
##  # Not needed anymore
##  cnr <- cnrData <- NULL;
##  x <- data$x;
##  y <- data$y;
##  naValue <- as.double(NA);
##  sigmas <- rep(naValue, times=nbrOfRegions);
##  for (kk in seq_len(nbrOfRegions)) {
##    keep <- which(regions[kk,1] < x & x <= regions[kk,2]);
##    t <- y[keep];
##    t <- diff(t);
##    t <- median(t, na.rm=TRUE)/sqrt(2);
##    sigmas[kk] <- t;
##  } # for (kk ...)
##  # Not needed anymore
##  x <- y <- t <- keep <- NULL;
##  aromaEstimates <- list(
##    stddevAll = sigma,
##    stddevRegions = sigmas
##  );
##  attr(fit, "aromaEstimates") <- aromaEstimates;
##  verbose && exit(verbose);

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
# 2010-03-02
# o Now segmentByMPCBS() can be used to segment data from a single source.
#   Added this to the (single-platform) segmentation example.
# 2010-02-28
# o Added more sanity checks.
# 2010-02-18
# o Added an MPCBS example().
# 2010-01-02
# o Created from segmentByCBS.R.
############################################################################
