# @title "Re-calculates the segmented profile using non-TumorBoost BAFs"
setMethodS3("unTumorBoost", "PairedPSCBS", function(fit, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enterf(verbose, "Relcalculating %s profile without TumorBoost", class(fit)[1L]);

  # Nothing to do?
  if (!fit$params$tbn) {
    verbose && cat(verbose, "Profile is not from TumorBoost signals. Skipping.");
    verbose && exit(verbose);
    return(fit);
  }

  verbose && enter(verbose, "Updating locus-level data");
  data <- getLocusData(fit);
  data$betaTN <- data$betaT;
  isHet <- !is.na(data$rho);
  data$rho[isHet] <- 2*abs(data$betaTN[isHet]-1/2);
  fit$data <- data;
  verbose && exit(verbose);

  verbose && enter(verbose, "Updating segments");
  segs <- getSegments(fit);
  segs$abCall <- NULL;
  segs$lohCall <- NULL;
  fit$output <- segs;
  verbose && exit(verbose);

  verbose && enter(verbose, "Updating parameters");
  params <- fit$params;
  params$tbn <- FALSE;
  params$deltaAB <- params$alphaAB <- NULL;
  params$deltaLowC1 <- params$alphaLowC1 <- NULL;
  fit$params <- params;
  verbose && exit(verbose);

  verbose && enter(verbose, "Resetting miscellaneous parameters and estimates");
  fit$changepoints <- NULL;
  fit$deshearC1C2 <- NULL;
  fit$cScaled <- NULL;
  fit$kappa <- NULL;
  fit$scale <- NULL;
  fit <- clearBootstrapSummaries(fit, verbose=less(verbose, 50));
  verbose && exit(verbose);

  verbose && enter(verbose, "Update segment levels");
  fit <- updateMeans(fit, verbose=less(verbose, 50));
  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;
}, protected=TRUE) # unTumorBoost()


##############################################################################
# HISTORY
# 2014-03-28
# o Added unTumorBoost() for PairedPSCBS.
##############################################################################
