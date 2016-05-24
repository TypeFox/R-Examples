###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod callAllelicBalanceByBAFs
#
# @title "Calls regions that are in allelic balance"
#
# \description{
#  @get "title" from the allele B fractions (BAF).
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A @see "PSCBS::PairedPSCBS" fit object as returned by
#     @see "PSCBS::segmentByPairedPSCBS".}
#   \item{maxScore}{A positive @double threshold.
#     If \code{"auto"}, the threshold is estimated empirically.}
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, an already called object is skipped, otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "PSCBS::PairedPSCBS" fit object
#   where columns for allelic imbalance scores and
#   p-values as well as allelic balance calls are added.
# }
#
# @examples "../incl/callAllelicBalanceByBAFs.PairedPSCBS.Rex"
#
# @author "HB, PN"
#
# \seealso{
#   Internally, @see "testAllelicBalanceByBAFs" is used.
#
#   Note that this AB caller differs from the default one in the
#   \pkg{PSCBS} package, cf. @see "PSCBS::callAB.PairedPSCBS".
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("callAllelicBalanceByBAFs", "PairedPSCBS", function(fit, maxScore="auto", ..., force=FALSE, cache=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'maxScore':
  if (is.character(maxScore)) {
    maxScore <- match.arg(maxScore);
  } else {
    maxScore <- Arguments$getDouble(maxScore, range=c(0,Inf));
  }

  # Extract segments
  segs <- as.data.frame(fit);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'cache':
  cache <- Arguments$getLogical(cache);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Nothing to do?
  if (!force && !is.null(segs$abCall)) {
    # Allelic balance segments are already called
    return(fit);
  }

  verbose && enter(verbose, "Calling allelic balance by BAFs");

  # Extract data
  data <- getLocusData(fit);
  betaTN <- data$betaTN;
  muN <- data$muN;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="callAllelicBalanceByBAFs", class=class(fit)[1],
    data=list(betaTN=betaTN, muN=muN, segments=as.data.frame(fit)),
    maxScore = maxScore,
    version="2011-11-02"
  );
  dirs <- c("aroma.cn", "ortho");
  if (!force) {
    res <- loadCache(key=key, dirs=dirs);
    if (!is.null(res)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(res);
    }
  }



  nbrOfSegments <- nrow(segs);
  naValue <- as.double(NA);
  df <- NULL;
  for (kk in seq_len(nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment #%d of %d", kk, nbrOfSegments));

    fitS <- extractDhSegment(fit, idx=kk, what="SNPs");
    if (is.null(fitS)) {
      verbose && cat(verbose, "A divider. Skipping.");
      dfKK <- data.frame(
        statistic=as.double(NA),
        p.value=as.double(NA)
      );
      df <- rbind(df, dfKK);
      verbose && exit(verbose);
      next;
    }

    verbose && print(verbose, fitS);

    dataS <- getLocusData(fitS);
    betaTN <- dataS$betaTN;
    muN <- dataS$muN;

    verbose && summary(verbose, betaTN);
    verbose && summary(verbose, muN);

    # AD HOC: For some unknown reason does resample() introduce NAs.
    # /HB 2010-09-15
    keep <- is.finite(betaTN) & is.finite(muN);
    keep <- which(keep);
    betaTN <- betaTN[keep];
    muN <- muN[keep];

    fitKK <- testAllelicBalanceByBAFs(betaTN, muN=muN);

    dfKK <- data.frame(
      statistic=fitKK$statistic,
      p.value=fitKK$p.value
    );

    df <- rbind(df, dfKK);

    verbose && exit(verbose);
  } # for (kk ...)
  rownames(df) <- NULL;

  colnames(df) <- c("ai", "ai.p.value");

  if (maxScore == "auto") {
    verbose && enter(verbose, "Estimating 'maxScore' cutoff empirically");
    d <- density(na.omit(df$ai), from=0, to=10, adjust=0.1);
    pvs <- .findPeaksAndValleys(d);
    verbose && print(verbose, pvs);

    type <- NULL; rm(list="type"); # To please R CMD check
    vs <- subset(pvs, type == "valley");
    v <- vs[1,];

    maxScore <- v$x;

    # Sanity check
    maxScore <- Arguments$getDouble(maxScore, range=c(0,Inf));

    attr(maxScore, "modelFit") <- list(density=d, pvs=pvs, v=v);

    verbose && cat(verbose, "maxScore: ", maxScore);
    verbose && exit(verbose);
  }

  params <- fit$params;
  paramsT <- list(maxScore=maxScore);
  params <- c(params, paramsT);

  df$abCall <- (df$ai <= maxScore);

  segs <- cbind(segs, df);

  fitC <- fit;
  fitC$output <- segs;
  fitC$params <- params;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    saveCache(key=key, dirs=dirs, fitC);
  }

  verbose && exit(verbose);

  fitC;
})



###########################################################################/**
# @RdocMethod callCopyNeutralRegions
#
# @title "Calls regions that are copy neutral"
#
# \description{
#  @get "title" from the allele B fractions (BAF).
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A @see "PSCBS::PairedPSCBS" fit object as returned by
#     @see "PSCBS::segmentByPairedPSCBS".}
#   \item{...}{Additional arguments passed to
#     @see "PSCBS::callCopyNeutral.PairedPSCBS".}
# }
#
# \value{
#   Returns a @see "PSCBS::PairedPSCBS" fit object
#   where a column with the copy-neutral call.
# }
#
# @examples "../incl/callCopyNeutralRegions.PairedPSCBS.Rex"
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("callCopyNeutralRegions", "PairedPSCBS", function(fit, ...) {
  # Call allelic balance or not, unless already done
  fit <- callAllelicBalanceByBAFs(fit, ...);
  callCopyNeutral(fit, flavor="TCN|AB", ...);
})



##############################################################################
# HISTORY
# 2012-02-24
# o CLEANUP: Now callCopyNeutralRegions() calls callCopyNeutral(), which
#   is now available in the PSCBS package.  The callCopyNeutralRegions()
#   method will eventually be deprecated for the favor of the PSCBS one.
# o Moved extractDhSegment() for PairedPSCBS to the PSCBS package.
# 2012-02-23
# o Moved drawC1C2Density() and plotC1C2Grid() to PairedPSCBS.PLOT3.R.
# o Made extractDhSegment() protected.
# 2011-12-15
# o Turned off default memoization for callAllelicBalanceByBAFs().
# o Now callAllelicBalanceByBAFs() for PairedPSCBS appends its parameter
#   settings to the ones in the PairedPSCBS object.  Before it dropped
#   them.
# 2011-10-16 [HB]
# o Now using getLocusData(fit) and getSegments(fit) where applicable.
# 2011-07-10 [HB]
# o Updated code to work with the new column names in PSCBS v0.11.0.
# 2010-10-26 [HB]
# o Now argument 'maxScore' for callAllelicBalanceByBAFs() defaults
#   to "auto", which corresponds to estimating the score empirically
#   from the allelic imbalance scores ('ai').
# o Added extractDhSegment() for PairedPSCBS.
# 2010-10-10 [HB]
# o Added memoization to callAllelicBalanceByBAFs().
# 2010-09-15 [HB]
# o Added Rdocs for callCopyNeutralRegions().
# 2010-09-09 [HB]
# o Added callCopyNeutralRegions() for PairedPSCBS.
# 2010-09-08 [HB]
# o Added subsetBySegments() for PairedPSCBS.
# o Added Rdocs with an example.
# o Added normalizeBAFsByRegions() for PairedPCSBS.
# o Created.
##############################################################################
