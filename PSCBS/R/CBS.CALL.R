###########################################################################/**
# @set "class=CBS"
# @RdocMethod callGainsAndLosses
#
# @title "Calls gains and losses"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{adjust}{A positive scale factor adjusting the sensitivity of the
#    caller, where a value less (greater) than 1.0 makes the caller
#    less (more) sensitive.}
#  \item{method}{A @character string specifying the calling algorithm to use.}
#  \item{...}{Additional/optional arguments used to override the default
#    parameters used by the caller.}
# }
#
# \value{
#  Returns a @see "PSCBS::CBS" object where @logical columns
#  'lossCall' and 'gainCall' have been appended to the segmentation table.
# }
#
# \section{The UCSF caller}{
#   If \code{method == "ucsf-mad"}, then segments are called using [1], i.e.
#   a segment is called gained or lost if its segment level is
#   at least two standard deviations away from the median segment level
#   on Chr1-22, where standard deviation is estimated using MAD.
#   Then same is done for \code{method == "ucsf-dmad"} with the difference
#   that the standard deviation is estimated using a robust first order
#   variance estimator.
# }
#
# \examples{
#   @include "../incl/segmentByCBS.Rex"
#   @include "../incl/segmentByCBS,calls.Rex"
# }
#
# @author "HB"
#
# \references{
#   [1] Fridlyand et al. \emph{Breast tumor copy number aberration
#       phenotypes and genomic instability}, BMC Cancer, 2006. \cr
# }
#
# \seealso{
#   @seemethod "callAmplifications".
#   @seemethod "callOutliers".
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("callGainsAndLosses", "CBS", function(fit, adjust=1.0, method=c("ucsf-mad", "ucsf-dmad"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'adjust':
  adjust <- Arguments$getDouble(adjust, range=c(0, Inf));

  # Argument 'method':
  method <- match.arg(method);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling segments that are gained or lost");

  userArgs <- list(...);

  params <- list();

  # Allocate calls
  naValue <- as.logical(NA);
  nbrOfSegments <- nbrOfSegments(fit, splitters=TRUE);
  segs <- getSegments(fit, splitters=TRUE);
  nbrOfRows <- nrow(segs);
  gainCalls <- lossCalls <- rep(naValue, times=nbrOfRows);

  verbose && cat(verbose, "Number of segments to be called: ", nbrOfSegments);
  verbose && cat(verbose, "Call method: ", method);

  if (is.element(method, c("ucsf-mad", "ucsf-dmad"))) {
    # Default arguments
    args <- list(
      chromosomes = intersect(getChromosomes(fit), c(0L, 1:22)),
      scale = 2.0
    );

    # Override by (optional) user-specified arguments
    for (key in names(userArgs)) {
      args[[key]] <- userArgs[[key]];
    }

    # Extract arguments
    chromosomes <- args$chromosomes;
    scale <- args$scale;

    # Argument check
    chromosomes <- Arguments$getVector(chromosomes, lengths=c(1,Inf));
    scale <- Arguments$getDouble(scale, range=c(0,Inf));

    # Estimate the whole-genome standard deviation of the TCNs
    if (method == "ucsf-mad") {
      sigma <- estimateStandardDeviation(fit, chromosomes=chromosomes,
                                         method="res", estimator="mad");
      sigmaKey <- "sigmaMAD";
    } else if (method == "ucsf-dmad") {
      sigma <- estimateStandardDeviation(fit, chromosomes=chromosomes,
                                         method="diff", estimator="mad");
      sigmaKey <- "sigmaDelta";
    } else {
      throw("INTERNAL ERROR: Unknown method: ", method);
    }

    # Sanity check
    sigma <- Arguments$getDouble(sigma, range=c(0,Inf));

    # Calculate the threshold
    tau <- scale * sigma;

    # Make more or less sensitive
    tau <- adjust * tau;

    verbose && cat(verbose, "Call parameters:");
    verbose && str(verbose, list(sigma=sigma, scale=scale, adjust=adjust));

    # Calculate segment levels using the median estimator
    fitT <- updateMeans(fit, avg="median")
    segsT <- getSegments(fitT, splitters=TRUE);
    mu <- segsT$mean;
    fitT <- segsT <- NULL; # Not needed anymore

    # The median segmented level
    muR <- median(mu, na.rm=TRUE);

    # The threshold for losses
    tauLoss <- muR - tau;

    # The threshold for gains
    tauGain <- muR + tau;

    # Call
    lossCalls <- (mu <= tauLoss);   # Losses
    gainCalls <- (mu >= tauGain);   # Gains

    # Call parameters used
    params$method <- method;
    params$adjust <- adjust;
    params[[sigmaKey]] <- sigma;
    params$scale <- scale;
    params$muR <- muR;
    params$tau <- tau;
    params$tauLoss <- tauLoss;
    params$tauGain <- tauGain;
  }

  verbose && cat(verbose, "Number of called segments: ", length(lossCalls));


  # Sanity check
  stopifnot(length(lossCalls) == nbrOfRows);
  stopifnot(length(gainCalls) == nbrOfRows);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update 'DNAcopy' object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (a) segmentation table
  segs <- getSegments(fit, splitters=TRUE);
  segs$lossCall <- lossCalls;
  segs$gainCall <- gainCalls;
  fit$output <- segs;

  # (b) parameters
  allParams <- fit$params;
  if (is.null(allParams)) {
    allParams <- list();
  }
  allParams$callGainsAndLosses <- params;
  fit$params <- allParams;

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return the updated 'CBS' object.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit;
}, private=TRUE) # callGainsAndLosses()





###########################################################################/**
# @RdocMethod callAmplifications
#
# @title "Calls (focal) amplifications"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{adjust}{A positive scale factor adjusting the sensitivity of the
#    caller, where a value less (greater) than 1.0 makes the caller
#    less (more) sensitive.}
#  \item{maxLength}{A @double scalar specifying the maximum length of a segment
#    in order for it to be considered a focal amplification.}
#  \item{method}{A @character string specifying the calling algorithm to use.}
#  \item{...}{Additional/optional arguments used to override the default
#    parameters used by the caller.}
#  \item{verbose}{@see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "PSCBS::CBS" object where @logical column
#  'amplificationCall' has been appended to the segmentation table.
# }
#
# \section{The UCSF caller}{
#   If \code{method == "ucsf-exp"}, then segments are called using [1], i.e.
#   a segment is called an amplification if ...
# }
#
# @author
#
# \references{
#   [1] Fridlyand et al. \emph{Breast tumor copy number aberration
#       phenotypes and genomic instability}, BMC Cancer, 2006. \cr
# }
#
# \seealso{
#   @seemethod "callGainsAndLosses".
#   @seemethod "callOutliers".
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("callAmplifications", "CBS", function(fit, adjust=1.0, maxLength=20e6, method=c("ucsf-exp"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'adjust':
  adjust <- Arguments$getDouble(adjust, range=c(0, Inf));

  # Argument 'maxLength':
  maxLength <- Arguments$getDouble(maxLength, range=c(0, Inf));

  # Argument 'method':
  method <- match.arg(method);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling segments that are amplified");

  userArgs <- list(...);

  params <- list();

  # Allocate calls
  naValue <- as.logical(NA);
  nbrOfSegments <- nbrOfSegments(fit, splitters=TRUE);
  calls <- rep(naValue, times=nbrOfSegments);

  verbose && cat(verbose, "Number of segments to be called: ", nbrOfSegments);
  verbose && cat(verbose, "Call method: ", method);

  if (method == "ucsf-exp") {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Call arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Default arguments
    args <- list(
      minLevel = 0.0,
      lambda   = 1.0,
      degree   = 3
    );

    # Override by (optional) user-specified arguments
    for (key in names(userArgs)) {
      args[[key]] <- userArgs[[key]];
    }

    # Extract arguments
    minLevel <- args$minLevel;
    lambda <- args$lambda;
    degree <- args$degree;

    # Validate arguments
    minLevel <- Arguments$getDouble(minLevel, range=c(-Inf, Inf));
    lambda <- Arguments$getDouble(lambda, range=c(0, Inf));
    degree <- Arguments$getDouble(degree, range=c(1, Inf));

    verbose && cat(verbose, "Call parameters:");
    verbose && str(verbose, list(minLevel=minLevel, lambda=lambda,
                                                  degree=degree));

    segs <- getSegments(fit, splitters=TRUE);

    verbose && cat(verbose, "Segments:");
    verbose && str(verbose, segs);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Rule #1: Only consider segments that are short enough
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # The lengths (in bp) of the segments
    start <- segs$start;
    end <- segs$end;
    length <- end - start; ## + 1L;
    keep1 <- (length <= maxLength);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Rule #2: Only consider segments that have a mean level
    #          that is large enough.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # The mean levels of the segments
    mu <- segs$mean;
    keep2 <- (mu >= minLevel);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Rule #3: Only consider segments that have a mean level
    #          that is much larger than either of the
    #          flanking segments.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # The mean levels of the flanking segments
    muL <- c(NA, mu[-nbrOfSegments]);
    muR <- c(mu[-1], NA);

    # The difference in mean levels to the flanking segments
    deltaL <- mu - muL;
    deltaR <- mu - muR;

    # The maximum difference to either of the flanking segments
    delta <- pmax(deltaL, deltaR, na.rm=TRUE);

    # The threshold for calling segments amplified
    tau <- exp(-lambda * mu^degree);

    # Make more or less sensitive
    tau <- adjust * tau;

    keep3 <- (delta >= tau);

    # Amplification calls
    calls <- (keep1 & keep2 & keep3);

    # Call parameters used
    params$method <- method;
    params$adjust <- adjust;
    params$maxLength <- maxLength;
    params$minLevel <- minLevel;
    params$lambda <- lambda;
    params$degree <- degree;
    params$tau <- tau;
  }

  verbose && cat(verbose, "Number of called segments: ", length(calls));

  # Sanity check
  stopifnot(length(calls) == nbrOfSegments);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update 'DNAcopy' object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (a) segmentation table
  segs <- getSegments(fit, splitters=TRUE);
  segs$amplificationCall <- calls;
  fit$output <- segs;

  # (b) parameters
  allParams <- fit$params;
  if (is.null(allParams)) {
    allParams <- list();
  }
  allParams$callAmplifications <- params;
  fit$params <- allParams;


  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return the updated 'CBS' object.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit;
}, private=TRUE) # callAmplifications()



###########################################################################/**
# @RdocMethod callOutliers
#
# @title "Calls outliers"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{adjust}{A positive scale factor adjusting the sensitivity of the
#    caller, where a value less (greater) than 1.0 makes the caller
#    less (more) sensitive.}
#  \item{method}{A @character string specifying the calling algorithm to use.}
#  \item{...}{Additional/optional arguments used to override the default
#    parameters used by the caller.}
# }
#
# \value{
#  Returns a @see "PSCBS::CBS" object where @logical columns
#  'negOutlierCall' and 'posOutlierCall' have been appended
#  to the segmentation table.
# }
#
# \section{The UCSF caller}{
#   If \code{method == "ucsf-mad"}, then loci are called using [1];
#  "Finally, to identify single technical or biological outliers such
#   as high level amplifications, the presence of the outliers within
#   a segment was allowed by assigning the original observed log2ratio
#   to the clones for which the observed values were more than four
#   tumor-specific MAD away from the smoothed values." [1; Suppl. Mat.]
# }
#
# @author "HB"
#
# \references{
#   [1] Fridlyand et al. \emph{Breast tumor copy number aberration
#       phenotypes and genomic instability}, BMC Cancer, 2006. \cr
# }
#
# \seealso{
#   @seemethod "callGainsAndLosses".
#   @seemethod "callAmplifications".
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("callOutliers", "CBS", function(fit, adjust=1.0, method=c("ucsf-mad"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'adjust':
  adjust <- Arguments$getDouble(adjust, range=c(0, Inf));

  # Argument 'method':
  method <- match.arg(method);


  userArgs <- list(...);

  params <- list();

  # Allocate calls
  nbrOfLoci <- nbrOfLoci(fit);
  naValue <- as.logical(NA);
  negOutlierCall <- posOutlierCall <- rep(naValue, times=nbrOfLoci);

  if (method == "ucsf-mad") {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Call arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Default arguments
    args <- list(
      scale = 4.0
    );

    # Override by (optional) user-specified arguments
    for (key in names(userArgs)) {
      args[[key]] <- userArgs[[key]];
    }

    # Extract arguments
    scale <- args$scale;

    # Validate arguments
    scale <- Arguments$getDouble(scale, range=c(0, Inf));


    # Genomic annotations
    data <- getLocusData(fit);
    chromosome <- data$chromosome;
    x <- data$x;

    # CN signals
    y <- data[,3];

    # Segmented CN signals
    yS <- extractSegmentMeansByLocus(fit);

    # CN residuals (relative to segment means)
    dy <- y - yS;

    segs <- getSegments(fit, splitters=TRUE);

    # Allocate per-segment SD estimates
    nbrOfSegments <- nbrOfSegments(fit);
    naValue <- NA_real_;
    sds <- rep(naValue, times=nbrOfSegments);

    naValue <- NA_real_;
    for (ss in seq(length=nbrOfSegments)) {
      seg <- segs[ss,];

      # Identify loci in current segment
      idxs <- which(seg$chromosome == chromosome &
                    seg$start <= x & x <= seg$end);

      # Sanity check
      idxs <- Arguments$getIndices(idxs, max=nbrOfLoci);

      # Extract CN residuals
      dySS <- dy[idxs];

      # Calculate MAD for segment
      sdSS <- mad(dySS, na.rm=TRUE);

      # Threshold for outliers
      tau <- scale * sdSS;

      # Make more or less sensitive
      tau <- adjust * tau;

      # Call outliers
      naValue <- as.logical(NA);
      callsSS <- rep(naValue, times=length(dySS));
      callsSS[-tau <= dySS & dySS <= +tau] <- 0L;
      callsSS[dySS > +tau] <- +1L;
      callsSS[dySS < -tau] <- -1L;

      # Record
      negOutlierCall[idxs] <- (callsSS < 0L);
      posOutlierCall[idxs] <- (callsSS > 0L);

      sds[ss] <- sdSS;
    } # for (ss ...)

    params$method <- method;
    params$adjust <- adjust;
    params$scale <- scale;
    params$sds <- sds;
  }


  # Sanity check
  stopifnot(length(negOutlierCall) == nbrOfLoci);
  stopifnot(length(posOutlierCall) == nbrOfLoci);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update 'DNAcopy' object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (a) segmentation table
  data <- getLocusData(fit);
  data$negOutlierCall <- negOutlierCall;
  data$posOutlierCall <- posOutlierCall;
  fit$data <- data;

  # (b) parameters
  allParams <- fit$params;
  if (is.null(allParams)) {
    allParams <- list();
  }
  allParams$callOutliers <- params;
  fit$params <- allParams;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return the updated 'CBS' object.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit;
}, private=TRUE) # callOutliers()



setMethodS3("extractCallsByLocus", "CBS", function(fit, ...) {
  # Extract locus data
  data <- getLocusData(fit, ...);

  nbrOfLoci <- nrow(data);

  # Extract segment data
  segs <- getSegments(fit, splitters=TRUE);

  # Identify segment calls
  callCols <- grep("Call$", colnames(segs));
  nbrOfCalls <- length(callCols);


  chromosome <- data$chromosome;
  x <- data$x;
  y <- data[,3];

  # Allocate locus calls
  naValue <- as.logical(NA);
  callsL <- matrix(naValue, nrow=nbrOfLoci, ncol=nbrOfCalls);
  colnames(callsL) <- colnames(segs)[callCols];
  callsL <- as.data.frame(callsL);

  # For each segment...
  for (ss in seq(length=nrow(segs))) {
    seg <- segs[ss,];
    idxs <- which(chromosome == seg$chromosome &
                  seg$start <= x & x <= seg$end);
    idxs <- Arguments$getIndices(idxs, max=nbrOfLoci);
    # Sanity check
##    stopifnot(length(idxs) == seg$nbrOfLoci);

    callsSS <- seg[callCols];
    for (cc in seq(length=nbrOfCalls)) {
      callsL[idxs,cc] <- callsSS[,cc];
    }
  } # for (ss ...)

  # The calls for loci that have missing annotations or observations,
  # should also be missing, i.e. NA.
  nok <- (is.na(chromosome) | is.na(x) | is.na(y));
  callsL[nok,] <- as.logical(NA);

  # Sanity check
  stopifnot(nrow(callsL) == nbrOfLoci);
  stopifnot(ncol(callsL) == nbrOfCalls);

  callsL;
}, private=TRUE) # extractCallsByLocus()



###########################################################################/**
# @RdocMethod getCallStatistics
#
# @title "Calculates various call statistics per chromosome"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{regions}{An optional @data.frame with columns "chromosome",
#     "start", and "end" specifying the regions of interest to calculate
#     statistics for.  If @NULL, all of the genome is used.}
#  \item{shrinkRegions}{If @TRUE, regions are shrunk to the support of
#     the data.}
#  \item{...}{Not used.}
#  \item{verbose}{@see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a CxK @data.frame, where C is the number of regions that
#  meet the criteria setup by argument \code{regions}
#  and (K-4)/2 is the number of call types.
#  The first column is the chromosome index, the second and the third
#  are the first and last position, and the fourth the length
#  (=last-first+1) of the chromosome.
#  The following columns contains call summaries per chromosome.
#  For each chromosome and call type, the total length of such calls
#  on that chromosome is reported together how large of a fraction
#  of the chromosome such calls occupy.
# }
#
# \details{
#   The estimators implemented here are based solely on the
#   segmentation results, which is very fast.
#   In the original proposal by Fridlyand et al. [1], the authors
#   estimates the parameters by converting segment-level calls back
#   to locus-level calls and there do the calculations.
#   The difference between the two approaches should be minor,
#   particularly for large density arrays.
# }
#
# @author "HB"
#
# \references{
#   [1] Fridlyand et al. \emph{Breast tumor copy number aberration
#       phenotypes and genomic instability}, BMC Cancer, 2006. \cr
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("getCallStatistics", "CBS", function(fit, regions=NULL, shrinkRegions=TRUE, ..., verbose=FALSE) {
  # To please R CMD check, cf. subset()
  chromosome <- NULL; rm(list="chromosome");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'regions':
  if (is.null(regions)) {
    # Get chromosome lengths
    regions <- getChromosomeRanges(fit)[,c("chromosome", "start", "end")];
  }
  regions <- as.data.frame(regions);
  stopifnot(all(is.element(c("chromosome", "start", "end"), colnames(regions))));
  stopifnot(!any(duplicated(regions$chromosome)));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calculating call statistics");
  segs <- getSegments(fit, splitters=FALSE);
  callTypes <- grep("Call$", colnames(segs), value=TRUE);
  verbose && cat(verbose, "Call types: ", hpaste(callTypes));
  if (length(callTypes) == 0) {
    throw("Cannot calculate call statistics. No calls have been made.");
  }

  verbose && cat(verbose, "Regions of interest:");
  verbose && str(verbose, regions);


  verbose && enter(verbose, "Filtering out segments within the requested regions");

  # Filter out segments within the requested regions
  segsT <- NULL;
  verbose && cat(verbose, "Number of segments (before): ", nrow(segs));

  for (rr in seq(length=nrow(regions))) {
    regionRR <- regions[rr,];
    chrRR <- regionRR[,"chromosome"];
    startRR <- regionRR[,"start"];
    endRR <- regionRR[,"end"];
    if (is.na(chrRR) || is.na(startRR) || is.na(endRR)) {
      next;
    }

    verbose && enter(verbose, sprintf("Region #%d of %d", rr, nrow(regions)));

    # Select regions that (at least) overlapping with the region
    segsRR <- subset(segs, chromosome == chrRR & start <= endRR & end >= startRR);

    verbose && cat(verbose, "Number of segments within region: ", nrow(segsRR));

    # Special case
    if (nrow(segsRR) == 0) {
      segsRR <- segs[1,][NA,];
      segsRR$chromosome <- chrRR;
      segsRR$start <- startRR;
      segsRR$end <- endRR;
      segsRR$nbrOfLoci <- 0L;
    }

    if (shrinkRegions) {
      range <- range(c(segsRR$start, segsRR$end), na.rm=TRUE);
      startRR <- max(startRR, range[1], na.rm=TRUE);
      endRR <- min(endRR, range[2], na.rm=TRUE);
      regions[rr,"end"] <- endRR;
      regions[rr,"start"] <- startRR;
    }

    # Adjust ranges
    segsRR$start[segsRR$start < startRR] <- startRR;
    segsRR$end[segsRR$end > endRR] <- endRR;

    segsRR$fullLength <- endRR - startRR; ## + 1L;

    segsT <- rbind(segsT, segsRR);

    verbose && exit(verbose);
  } # for (rr ...)

  segs <- segsT;

  # Order by chromosome
  o <- order(segs$chromosome);
  segs <- segs[o,];

  verbose && cat(verbose, "Number of segments (after): ", nrow(segs));
  verbose && str(verbose, segs);

  verbose && exit(verbose);


  verbose && enter(verbose, "Calculating total length per call and chromosome");
  # Sum length of calls per type and chromosome
  segs$length <- segs[,"end"] - segs[,"start"];  ## + 1L;
  res <- lapply(callTypes, FUN=function(type) {
    coeffs <- as.integer(segs[,type]);
    lens <- coeffs * segs$length;
    lens <- by(lens, INDICES=segs$chromosome, FUN=sum, na.rm=TRUE);
    as.vector(lens);
  });
  names(res) <- gsub("Call$", "Length", callTypes);
  res1 <- as.data.frame(res);
  verbose && str(verbose, res);
  verbose && exit(verbose);

  # Extract selected regions
  idxs <- match(unique(segs$chromosome), regions$chromosome);
  regionsT <- regions[idxs,];

  # Sanity check
  stopifnot(nrow(regionsT) == nrow(res1));


  verbose && enter(verbose, "Calculating fractions per region");
  # Calculate lengths
  regionsT$length <- regionsT[,"end"] - regionsT[,"start"]; ## + 1L;
  stopifnot(all(regionsT$length >= 0));

  res2 <- res1 / regionsT[,"length"];
  names(res2) <- gsub("Call$", "Fraction", callTypes);
  verbose && exit(verbose);

  res3 <- cbind(res1, res2);

  res <- regionsT;
  if (nrow(res3) > 0) {
    res <- cbind(res, res3);
  }
  rownames(res) <- NULL;

  res <- cbind(label=I(sprintf("chr%d", res[,"chromosome"])), res);

  # Sanity checks
  resT <- res[,grep("Fraction", colnames(res))];
  for (key in colnames(resT)) {
    rho <- resT[,key];
    stopifnot(all(rho >= 0, na.rm=TRUE));
    stopifnot(all(rho <= 1, na.rm=TRUE));
  }

  stopifnot(nrow(res) == nrow(regions));

  verbose && str(verbose, res);
  verbose && exit(verbose);

  res;
}, protected=TRUE) # getCallStatistics()



###########################################################################/**
# @RdocMethod getFractionOfGenomeLost
# @aliasmethod getFractionOfGenomeGained
# @aliasmethod getFractionOfGenomeAltered
# @aliasmethod getFGL
# @aliasmethod getFGG
# @aliasmethod getFGA
#
# @title "Calculates the fraction of the genome lost, gained, or aberrant either way"
#
# \description{
#  @get "title" (in sense of total copy numbers),
#  using definitions closely related to those presented in [1].
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @double in [0,1].
# }
#
# @author "HB"
#
# \references{
#   [1] Fridlyand et al. \emph{Breast tumor copy number aberration
#       phenotypes and genomic instability}, BMC Cancer, 2006. \cr
# }
#
# \seealso{
#   Internally, @seemethod "getCallStatistics" is used.
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("getFractionOfGenomeLost", "CBS", function(fit, ...) {
  stats <- getCallStatistics(fit, ...);
  mean(stats$lossFraction, na.rm=TRUE);
}, protected=TRUE)

setMethodS3("getFractionOfGenomeGained", "CBS", function(fit, ...) {
  stats <- getCallStatistics(fit, ...);
  mean(stats$gainFraction, na.rm=TRUE);
}, protected=TRUE)

setMethodS3("getFractionOfGenomeAltered", "CBS", function(fit, ...) {
  getFractionOfGenomeLost(fit, ...) + getFractionOfGenomeGained(fit, ...);
}, protected=TRUE)

# Shortcuts
setMethodS3("getFGL", "CBS", function(fit, ...) {
  getFractionOfGenomeLost(fit, ...);
}, protected=TRUE)

setMethodS3("getFGG", "CBS", function(fit, ...) {
  getFractionOfGenomeGained(fit, ...);
}, protected=TRUE)

setMethodS3("getFGA", "CBS", function(fit, ...) {
  getFractionOfGenomeAltered(fit, ...);
}, protected=TRUE)




setMethodS3("isWholeChromosomeGained", "CBS", function(fit, minFraction=0.99, ...) {
  # Argument 'minFraction':
  minFraction <- Arguments$getDouble(minFraction, range=c(0,1));

  stats <- getCallStatistics(fit, ...);
  calls <- stats$gainFraction;
  if (is.null(calls)) {
    return(rep(NA, times=nbrOfChromosomes(fit)));
  }

  res <- (calls >= minFraction);
  names(res) <- stats$chromosome;
  attr(res, "minFraction") <- minFraction;

  res;
}, protected=TRUE) # isWholeChromosomeGained()


setMethodS3("isWholeChromosomeLost", "CBS", function(fit, minFraction=0.99, ...) {
  # Argument 'minFraction':
  minFraction <- Arguments$getDouble(minFraction, range=c(0,1));

  stats <- getCallStatistics(fit, ...);
  calls <- stats$lossFraction;
  if (is.null(calls)) {
    return(rep(NA, times=nbrOfChromosomes(fit)));
  }

  res <- (calls >= minFraction);
  names(res) <- stats$chromosome;
  attr(res, "minFraction") <- minFraction;

  res;
}, protected=TRUE) # isWholeChromosomeLost()


setMethodS3("nbrOfLosses", "CBS", function(fit, ...) {
  stats <- getSegments(fit, ...);
  calls <- stats$lossCall;
  if (is.null(calls)) {
    return(NA_integer_);
  }
  sum(calls, na.rm=TRUE);
}, protected=TRUE)


setMethodS3("nbrOfGains", "CBS", function(fit, ...) {
  stats <- getSegments(fit, ...);
  calls <- stats$gainCall;
  if (is.null(calls)) {
    return(NA_integer_);
  }
  sum(calls, na.rm=TRUE);
}, protected=TRUE)


setMethodS3("nbrOfAmplifications", "CBS", function(fit, ...) {
  stats <- getSegments(fit, ...);
  calls <- stats$amplificationCall;
  if (is.null(calls)) {
    return(NA_integer_);
  }
  sum(calls, na.rm=TRUE);
}, protected=TRUE)


setMethodS3("getCallStatisticsByArms", "CBS", function(fit, genomeData, ...) {
  # To please/trick R CMD check
  chromosome <- x <- NULL; rm(list=c("chromosome", "x"));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'genomeData':
  genomeData <- as.data.frame(genomeData);



  # Subset 'regions' by chromosomes segmented
  keep <- is.element(genomeData$chromosome, getChromosomes(fit));
  genomeData <- genomeData[keep,];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # p-arm
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  regions <- getChromosomeRanges(fit);
  regions$end <- genomeData$centroStart;
  regions$start <- pmin(regions$start, regions$end);

  # Shrink regions
  for (rr in seq(length=nrow(regions))) {
    chr <- regions[rr,"chromosome"];
    x0 <- regions[rr,"start"];
    x1 <- regions[rr,"end"];
    xs <- subset(fit$data, chromosome == chr & x0 <= x & x <= x1)$x;
    if (length(xs) > 0) {
      range <- range(xs, na.rm=TRUE);
      x0 <- max(c(x0, range[1]), na.rm=TRUE);
      x1 <- min(c(x1, range[2]), na.rm=TRUE);
      regions[rr,"start"] <- x0;
      regions[rr,"end"] <- x1;
    }
  } # for (rr ...)
  regions[,"length"] <- regions[,"end"] - regions[,"start"]; ## + 1L;
  callStats <- getCallStatistics(fit, regions=regions);
  callStats$label <- sprintf("%sp", callStats$label);
  callStatsP <- callStats;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # q-arm
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  regions <- getChromosomeRanges(fit);
  regions$start <- genomeData$centroEnd;
  regions$end <- pmax(regions$end, regions$start);

  # Shrink regions
  for (rr in seq(length=nrow(regions))) {
    chr <- regions[rr,"chromosome"];
    x0 <- regions[rr,"start"];
    x1 <- regions[rr,"end"];
    xs <- subset(fit$data, chromosome == chr & x0 <= x & x <= x1)$x;
    if (length(xs) > 0) {
      range <- range(xs, na.rm=TRUE);
      x0 <- max(c(x0, range[1]), na.rm=TRUE);
      x1 <- min(c(x1, range[2]), na.rm=TRUE);
      regions[rr,"start"] <- x0;
      regions[rr,"end"] <- x1;
    }
  } # for (rr ...)
  regions[,"length"] <- regions[,"end"] - regions[,"start"]; ## + 1L;

  callStats <- getCallStatistics(fit, regions=regions);
  callStats$label <- sprintf("%sq", callStats$label);
  callStatsQ <- callStats;



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Merge
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  callStats <- rbind(callStatsP, callStatsQ);

  # Not needed anymore
  regions <- callStatsP <- callStatsQ <- NULL;

  # Reorder
  o <- order(callStats$chromosome, callStats$start);
  callStats <- callStats[o,];

  callStats;
}, protected=TRUE); # getCallStatisticsByArms()


setMethodS3("callArms", "CBS", function(fit, genomeData, minFraction=0.95, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'minFraction':
  minFraction <- Arguments$getDouble(minFraction, range=c(0,1));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # p-arm
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  callStats <- getCallStatisticsByArms(fit, genomeData=genomeData);

  callTypes <- grep("Fraction", colnames(callStats), value=TRUE);
  callTypes <- gsub("Fraction", "", callTypes);

  keys <- sprintf("%sFraction", callTypes);
  rhos <- callStats[,keys];
  calls <- (rhos >= minFraction);
  colnames(calls) <- sprintf("%sCall", callTypes);

  callStats <- cbind(callStats, calls);

  callStats;
}, protected=TRUE); # callArms()




###########################################################################/**
# @RdocMethod mergeNonCalledSegments
#
# @title "Merge neighboring segments that are not called"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
#  \item{verbose}{@see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an object of the same class
#   with the same of fewer number of segments.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("mergeNonCalledSegments", "CBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Merging neighboring segments that are not called");

  # Identify call columns
  segs <- getSegments(fit, splitters=TRUE);
  keep <- grep("Call$", colnames(segs));
  nbrOfCalls <- length(keep);

  # Sanity check
  stopifnot(nbrOfCalls > 0);

  chromosomes <- getChromosomes(fit);
  fitList <- list();
  for (cc in seq(along=chromosomes)) {
    chromosome <- chromosomes[cc];
    verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", cc, chromosome, length(chromosomes)));


    fitCC <- extractChromosome(fit, chromosome=chromosome);
    n0 <- nbrOfSegments(fitCC);

    # Until no more neighboring non-called segments exists
    while (TRUE) {
      segs <- getSegments(fitCC, splitters=TRUE);
      calls <- as.matrix(segs[,keep]);

      # Find two neighboring segments that are not called
      isCalled <- rowAnys(calls, na.rm=TRUE);
      verbose && printf(verbose, "Number of segments not called: %d of %d\n", sum(!isCalled, na.rm=TRUE), length(isCalled));

      notCalled <- which(!isCalled);
      delta <- diff(notCalled);
      left <- notCalled[which(delta == 1)[1]];

      # No more segments to merge?
      if (is.na(left)) {
        break;
      }

      fitCC <- mergeTwoSegments(fitCC, left=left);
    } # while (...)

    n1 <- nbrOfSegments(fitCC);
    verbose && printf(verbose, "Number of segments merged: %d of %d\n", n0-n1, n0);

    fitList[[cc]] <- fitCC;
    verbose && exit(verbose);
  } # for (cc ...)

  verbose && enter(verbose, "Building result");
  res <- Reduce(append, fitList);
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}, protected=TRUE); # mergeNonCalledSegments()


setMethodS3("estimateDeltaCN", "CBS", function(fit, flavor=c("density(TCN)", "density(dTCN)", "dTCN"), adjust=0.3, ..., verbose=FALSE) {
  # This will load the 'aroma.light' namespace, if not already done.
  findPeaksAndValleys <- .use("findPeaksAndValleys", package="aroma.light");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'adjust':
  adjust <- Arguments$getDouble(adjust, range=c(0,10));


  if (flavor == "density(TCN)") {
    # Get segment mean levels
    segs <- getSegments(fit, splitters=FALSE);
    x <- segs$mean;
    w <- segs$nbrOfLoci;

    # Drop missing values
    keep <- is.finite(x) & is.finite(w);
    x <- x[keep];
    w <- w[keep];
    keep <- NULL; # Not needed anymore

    # Normalize weights
    w <- w / sum(w, na.rm=TRUE);

    # Estimate density
    d <- density(x, weights=w, adjust=adjust);

    w <- NULL; # Not needed anymore

    # Find peaks
    pv <- findPeaksAndValleys(d, ...);
    type <- NULL; rm(list="type"); # To please R CMD check
    p <- subset(pv, type == "peak");
    px <- p$x;
    pw <- p$density;

    # Distance between peaks
    dx <- diff(px);
    # Weights "between" peaks (AD HOC: sum up peak weights)
    dw <- pw[-length(pw)] + pw[-1L];

    deltaCN <- weighted.mean(dx, w=dw);
  } else if (flavor == "density(dTCN)") {
    # Get change-point magnitudes
    x <- getChangePoints(fit)[[1L]];
    x <- abs(x);

    # Drop missing values
    keep <- is.finite(x);
    x <- x[keep];
    keep <- NULL; # Not needed anymore


    # Estimate density
    d <- density(x, adjust=adjust);

    # Find peaks
    pv <- findPeaksAndValleys(d, ...);
    type <- NULL; rm(list="type"); # To please R CMD check
    p <- subset(pv, type == "peak");
    px <- p$x;
    pw <- p$density;

    # Distance between peaks
    dx <- diff(px);
    # Weights "between" peaks (AD HOC: sum up peak weights)
    dw <- pw[-length(pw)] + pw[-1L];

    throw("Still not implemented.");
  } else if (flavor == "dTCN") {
    # Get change-point magnitudes
    x <- getChangePoints(fit)[[1L]];
    x <- abs(x);

    # Drop missing values
    keep <- is.finite(x);
    x <- x[keep];
    keep <- NULL; # Not needed anymore

    deltaCN <- median(x);
  }

  # Sanity check
  deltaCN <- Arguments$getDouble(deltaCN, range=c(0, Inf));

  deltaCN;
}, protected=TRUE)



setMethodS3("encodeCalls", "data.frame", function(calls, flavor="UCSF", ...) {
  # Argument 'calls':
  stopifnot(all(is.element(c("chromosome", "x"), colnames(calls))));
  stopifnot(all(is.element(c("lossCall", "gainCall"), colnames(calls))));

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  calls0 <- calls;

  # Allocate
  calls <- rep(NA_real_, times=nrow(calls0));

  # Encode loss, neutral and gain (required)
  calls[!calls0$gainCall & !calls0$lossCall] <- 0;
  calls[calls0$gainCall] <- +1;
  calls[calls0$lossCall] <- -1;

  # Encode amplifications, if any/called.
  idxs <- which(calls0$amplificationCall);
  calls[idxs] <- +9;

  # Encode negative and positive outliers, if any/called.
  idxs <- which(calls0$negOutlierCall);
  calls[idxs] <- calls[idxs] - 0.1;

  idxs <- which(calls0$posOutlierCall);
  calls[idxs] <- calls[idxs] + 0.1;

  calls;
}, protected=TRUE) # encodeCalls()


setMethodS3("callGLAO", "CBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Call gains, losses, amplifications and (negative and positive) outliers");
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments(fit));

  # Call segments
  fitC <- callGainsAndLosses(fit, ..., verbose=verbose);
  fitC <- callAmplifications(fitC, ..., verbose=verbose);

  # Call loci, i.e. locus-level negative and positive outliers
  fitC <- callOutliers(fitC, ..., verbose=verbose);
  verbose && print(verbose, fitC);

  verbose && exit(verbose);

  fitC;
}, protected=TRUE) # callGLAO()


############################################################################
# HISTORY:
# 2013-12-17
# o Added argument 'flavor' to estimateDeltaCN() for CBS, which specifies
#   the type of estimator to use.
# 2013-11-27
# o Added callGLAO() for CBS.
# o Added encodeCalls() for 'data.frame' object returned by
#   getLocusData(..., addCalls=TRUE).
# 2013-11-23
# o BUG FIX: estimateDeltaCN() assumed aroma.light was loaded.
# 2013-11-14
# o Added estimateDeltaCN() for CBS.
# o BUG FIX: callGainsAndLosses() for CBS would not estimate the median
#   median CN level correctly if there were "empty" segments (e.g. gaps).
#   This was/is due to a bug in segments.summary() of the DNAcopy package.
#   Instead, we are now calculating the segment median levels ourselves.
# 2012-01-24
# o ROBUSTNESS: Now getCallStatistics() for CBS asserts that calls have
#   been made.  If not, an exception is thrown.
# 2011-12-13
# o Added "ucsf-dmad" to argument 'method' for callGainsAndLosses() of CBS.
# 2011-12-12
# o Now extractCallsByLocus() for CBS passes arguments
#   '...' to getLocusData().
# 2011-10-23
# o BUG FIX: getCallStatisticsByArms() for CBS would thrown a error if
#   argument 'genomeData' did not contain exactly the same chromosomes
#   as in the CBS object.
# o BUG FIX: The length of a segment must be defined as 'end-start' and
#   not 'end-start+1' so that the the total length of all segments
#   adds up correctly.
# o Added verbose output to callGainsAndLosses() and callAmplifications().
# o BUG FIX: callAmplifications() for CBS generated an error, if
#   more than one chromosome were called.
# 2011-10-08
# o Added mergeNonCalledSegments() for CBS.
# 2011-10-07
# o Now getCallStatistics() for CBS always return statistics for
#   all regions requested, even empty ones.
# o Now getCallStatistics() for CBS also returns a 'label' column.
# o Added getCallStatisticsByArms() and callArms() for CBS.
# 2011-10-06
# o Added optional argument 'regions' to getCallStatistics() for CBS.
# o Now getCallStatistics() for CBS also returns 'start' and 'end'
#   position of each chromosome.
# 2011-10-03
# o DOCUMENTATION: Added more help pages.
# 2011-10-02
# o DOCUMENTATION: Added an Rdoc help page for getFractionOfGenomeLost(),
#   getFractionOfGenomeGained(), getFractionOfGenomeAltered(), getFGL(),
#   getFGG() and getFGA().
# 2011-09-05
# o Added getCallStatistics() for CBS.
# 2011-09-04
# o Added extractCallsByLocus() for CBS.
# o Adopted the calling methods from ditto of the DNAcopy class.
# 2011-09-01
# o Now callGainsAndLosses() returns a DNAcopy where the segmentation
#   table has the new column 'tcnCall'.
# 2011-08-19
# o Added argument 'callParams' to plotTracks() for DNAcopy.
# 2011-07-24
# o Added callOutliers().
# 2011-07-21
# o Now amplified segments are also highlighted.
# 2011-07-20
# o Added callAmplifications().
# 2011-07-20
# o Now callGainsAndLosses() estimates the noise level on autosomes only.
# o Now callGainsAndLosses() returns parameters used.
# o Updated callGainsAndLosses() to estimate the std. dev. as the
#   MAD of the *residuals* (not the absolute) values.
# o Added support for estimateStandardDeviation(..., method="res").
# o Added extractSegmentMeansByLocus().
# o Added drawCentromeres().
# 2011-07-18
# o Added getSampleNames().
# o Added plotTracks() for DNAcopy.
# o Added callGainsAndLosses() to DNAcopy objects.
# o Added nbrOfSegments(), nbrOfLoci() and nbrOfSamples().
# 2011-07-17
# o Added estimateStandardDeviation() to DNAcopy objects.
############################################################################
