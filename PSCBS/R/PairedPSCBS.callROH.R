##########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod callROH
# @alias callROH.NonPairedPSCBS
#
# @title "Calls segments that are in ROH"
#
# \description{
#  @get "title", i.e. that have no (true) heterozygous genotypes.
#  Run of homozygosity (ROH) is a property of the normal (germline) sample.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments passed to @see "testROH".}
#   \item{updateMeans}{If @TRUE, DH and (C1,C2) mean levels are set
#    to @NA for segments called ROH, otherwise not.}
#   \item{force}{If @FALSE, and ROH calls already exits,
#    then nothing is done, otherwise the calls are done.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" object with ROH calls.
# }
#
# @author "PN, HB"
#
# \seealso{
#   Internally, @see "testROH" is used.
#   To call allelic balance (AB) see @seemethod "callAB".
#   To call loss of heterozygosity (LOH) see @seemethod "callLOH".
# }
#*/###########################################################################
setMethodS3("callROH", "PairedPSCBS", function(fit, ..., updateMeans=TRUE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling ROH");

  # Already done?
  segs <- getSegments(fit);
  calls <- segs$rohCall;
  if (!force && !is.null(calls)) {
    return(invisible(fit));
  }

  nbrOfSegments <- nrow(segs);
  calls <- rep(NA, times=nbrOfSegments);
  if (is.null(calls)) {
    segs <- cbind(segs, rohCall=calls);
  }
  delta <- NA_real_;

  # For each segment...
  for (ss in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment #%d of %d", ss, nbrOfSegments));

    fitT <- extractSegment(fit, ss);

    # Call only "non-splitter" segments
    if (nbrOfSegments(fitT) > 0L) {
      callSS <- callROHOneSegment(fitT, ..., verbose=less(verbose, 1));
      calls[ss] <- callSS;
      if (is.na(delta) && !is.na(callSS)) {
        delta <- attr(callSS, "delta");
      }
    }

    verbose && exit(verbose);
  } # for (ss ...)

  verbose && cat(verbose, "ROH calls:");
  verbose && str(verbose, calls);
  verbose && print(verbose, summary(calls));

  segs$rohCall <- calls;

  fit$output <- segs;

  # Append parameters
  params <- fit$params;
  params$deltaROH <- delta;
  fit$params <- params;

  # Set DH and (C1,C2) mean levels to NA?
  if (updateMeans) {
    fit <- updateMeans(fit, from="segments", adjustFor="roh",
                                             verbose=less(verbose, 20));
  }

  verbose && exit(verbose);

  invisible(fit);
}) # callROH()


# This method calls ROH for a single-segment PairedPSCBS object
setMethodS3("callROHOneSegment", "PairedPSCBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling ROH for a single segment");

  # Make sure that there is only a single segment in this object
  stopifnot(nbrOfSegments(fit, splitters=TRUE) == 1L);


  # Extract the locus-level data for the segment tested
  data <- getLocusData(fit);


  # Keep only SNPs:
  # SNPs are identifies as those loci that have non-missing
  # 'betaTN' & 'muN', cf. segmentByPairedPSCBS().
  isSnp <- (!is.na(data$betaTN) & !is.na(data$muN));
  nbrOfSnps <- sum(isSnp);
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSnps);
  data <- data[isSnp,];

  # Extract that SNP signals used for calling ROH
  betaN <- data$betaN;
  muN <- data$muN;
  csN <- data$csN;  # Genotyping confidence scores, if available

  # Test for ROH
  fit <- testROH(muN=muN, csN=csN, betaN=betaN, ..., verbose=less(verbose, 10));

  # Get the ROH call (TRUE, FALSE, or NA)
  call <- fit;

  verbose && exit(verbose);

  call;
}, private=TRUE)


##############################################################################
# HISTORY
# 2014-03-29 [HB]
# o BUG FIX: In rare cases, callROH() could throw "Error in if (is.na(delta))
#   { : argument is of length zero".
# 2012-05-30 [HB]
# o Now callROH() records paramter 'deltaROH' in the results.
# 2011-11-26 [HB]
# o Added argument 'updateMeans=TRUE' to callROH() for PairedPSCBS.
# 2011-11-12 [HB]
# o BUG FIX: ROH calls should be stored in column 'rohCall' (not 'rohCalls').
# 2011-11-04 [HB]
# o Added callROH() for PairedPSCBS.
# o Created.
##############################################################################
