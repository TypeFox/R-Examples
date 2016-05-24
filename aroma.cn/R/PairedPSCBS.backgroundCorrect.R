setMethodS3("translateC1C2", "PairedPSCBS", function(fit, dC1=0, dC2=0, sC1=1, sC2=1, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dC1 <- Arguments$getNumeric(dC1);
  dC2 <- Arguments$getNumeric(dC2);
  sC1 <- Arguments$getNumeric(sC1);
  sC2 <- Arguments$getNumeric(sC2);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Scaling and translating by (C1,C2)");
  verbose && cat(verbose, "sC1: ", sC1);
  verbose && cat(verbose, "sC2: ", sC2);
  verbose && cat(verbose, "dC1: ", dC1);
  verbose && cat(verbose, "dC2: ", dC2);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  segs <- getSegments(fit, splitters=FALSE);
  stopifnot(!is.null(segs));

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # (C1,C2,...)
  X <- extractC1C2(fit, splitters=FALSE);

  # (C1,C2)
  C1C2 <- X[,1:2, drop=FALSE];

  C1C2[,1] <- dC1 + sC1*C1C2[,1];
  C1C2[,2] <- dC2 + sC2*C1C2[,2];

  verbose && enter(verbose, "(C1,C2) to (TCN,DH)");
  # (C1,C2) -> (TCN,DH)
  gamma <- rowSums(C1C2, na.rm=TRUE);
  dh <- 2*(C1C2[,2]/gamma - 1/2);
  verbose && exit(verbose);

  # Update segmentation means
  segs[,"tcnMean"] <- gamma;
  segs[,"dhMean"] <- dh;
  segs[,"c1Mean"] <- C1C2[,1];
  segs[,"c2Mean"] <- C1C2[,2];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitO <- fit;
  fitO$output <- segs;

  verbose && exit(verbose);

  fitO;
})



setMethodS3("transformC1C2", "PairedPSCBS", function(fit, fcn, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fcn':
  stopifnot(is.function(fcn));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Transform (C1,C2) by a function");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  segs <- getSegments(fit, splitters=FALSE);
  stopifnot(!is.null(segs));

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # (C1,C2,...)
  X <- extractC1C2(fit, splitters=FALSE);

  # (C1,C2)
  C1C2 <- X[,1:2, drop=FALSE];

  C1C2 <- fcn(C1C2, ...);

  verbose && enter(verbose, "(C1,C2) to (TCN,DH)");
  # (C1,C2) -> (TCN,DH)
  gamma <- rowSums(C1C2, na.rm=TRUE);
  dh <- 2*(C1C2[,2]/gamma - 1/2);
  verbose && exit(verbose);

  # Update segmentation means
  segs[,"tcnMean"] <- gamma;
  segs[,"dhMean"] <- dh;
  segs[,"c1Mean"] <- C1C2[,1];
  segs[,"c2Mean"] <- C1C2[,2];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitO <- fit;
  fitO$output <- segs;

  verbose && exit(verbose);

  fitO;
})



setMethodS3("estimateC2Bias", "PairedPSCBS", function(fit, ...) {
  # Identify region in allelic balance
  segs <- getSegments(fit, splitters=FALSE);
  abCall <- segs$abCall;
  if (is.null(abCall)) {
    throw("Allelic balance has not been called.");
  }
  idxs <- which(segs$abCall);
  segs <- segs[idxs,,drop=FALSE];

  # Extract (TCN,DH)
  gamma <- segs[, "tcnMean"];
  rho <- segs[, "dhMean"];

  # Calculate (C1,C2)
  C1 <- 1/2 * (1 - rho) * gamma;
  C2 <- gamma - C1;

  # Calculate bias in C2
  dC2 <- C2 - C1;

  # Calculate weighted average of all C2 biases
  w <- segs[,"dhNbrOfLoci"];
  w <- w / sum(w, na.rm=TRUE);
  dC2 <- weightedMedian(dC2, w=w);

  dC2;
}) # estimateC2Bias()



setMethodS3("backgroundCorrect", "PairedPSCBS", function(fit, targetMedian=c("same", "none"), ...) {
  # Argument 'targetMedian':
  if (is.character(targetMedian)) {
    targetMedian <- match.arg(targetMedian);
  } else {
    targetMedian <- Arguments$getDouble(targetMedian, range=c(0,Inf));
  }

  modelFit <- list();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (1) Estimate background (e.g. normal contamination and more)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  kappa <- estimateKappa(fit);
  modelFit$kappa <- kappa;
  fitBG <- translateC1C2(fit, dC1=-kappa, dC2=-kappa);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (2) Rescale
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate current median
  modelFit$targetMedian <- targetMedian;
  if (!identical(targetMedian, "none")) {
    segsBG <- getSegments(fitBG);
    count <- segsBG[["tcnNbrOfLoci"]];
    # chrs <- 1:22;
    # if (!is.null(chrs)) segsBG <- subset(segsBG, chromosome %in% chrs);
    CTBG <- segsBG[["tcnMean"]];
    muBG <- weightedMedian(CTBG, w=count);
    muBG <- Arguments$getDouble(muBG, range=c(0,Inf));
    modelFit$muBG <- muBG;

    if (identical(targetMedian, "same")) {
      segs <- getSegments(fit);
      # if (!is.null(chrs)) segs <- subset(segs, chromosome %in% chrs);
      CT <- segs[["tcnMean"]];
      mu <- weightedMedian(CT, w=count);
      mu <- Arguments$getDouble(mu, range=c(0,Inf));
    } else if (is.numeric(targetMedian)) {
      mu <- targetMedian;
    }
    modelFit$mu <- mu;

    scale <- mu/muBG;
    modelFit$scale <- scale;

    fitBG <- translateC1C2(fitBG, sC1=scale, sC2=scale);
  } # if (!identical(targetMedian, "none"))


  # Store model fit
  postMethods <- fit$postMethods;
  if (is.null(postMethods)) postMethods <- list();
  postMethods$backgroundCorrection <- list(modelFit=modelFit);
  fitBG$postMethods <- postMethods;

  fitBG;
}) # backgroundCorrect()



##############################################################################
# HISTORY
# 2013-08-21
# o Retired this version of deShearC1C2() by renaming it.
# 2012-09-22
# o Added backgroundCorrect() for PairedPSCBS.
# o Now fitDeltaC1C2ShearModel() "handles" also single elements.
# o BUG FIX: deShearC1C2() would incorrectly set the TCN level to zero
#   if (C1,C2) = (NA,NA).
# o Now deShearC1C2() stores the modelFit in postMethods$deShearC1C2$modelFit.
# o CORRECTION: The peak caller of fitDeltaXYShearModel() swapped the labels
#   of the horizontal and vertical peaks. Labels are only used for display.
# 2012-09-21 [HB]
# o BUG FIX: fitDeltaC1C2ShearModel(... dropABChangePoints=TRUE) could
#   generate an incorrect number change-point weights if dropping.
# o Added argument 'onError' to fitDeltaXYShearModel() for matrix.
# o Now attribute 'modelFit' fitDeltaXYShearModel() also contains
#   debug$cpAngleDensity.
# 2012-09-21 [PN]
# o Now fitDeltaC1C2ShearModel() will apply a fast naive AB caller,
#   iff AB calls are not available.
# o BUG FIX: horizontal/vertical axes were swapped for direction '|_'
# 2012-09-19 [HB]
# o MEMORY/"BUG FIX": fitDeltaXYShearModel() would return "H" deshearing
#   functions that each would hold the very large fitDeltaXYShearModel()
#   call environment.  Now we make sure each of those objects only
#   contains the few scalar parameter estimates needed.
# 2012-09-19 [PN]
# o Added direction '|_': similar to '|-' but enforces equal amount of
#   de-shearing in both dimensions and does not change TCN.
# 2011-10-17 [HB]
# o Added argument 'dropABChangePoints' to fitDeltaC1C2ShearModel().
# 2011-10-16 [HB]
# o Now deShearC1C2(), translateC1C2() and transformC1C2() also update
#   C1 and C2 mean levels.
# 2011-07-10 [HB]
# o Updated code to work with the new column names in PSCBS v0.11.0.
# 2010-10-20 [HB]
# o Now fitDeltaXYShearModel() uses both -pi/4 and +pi/4 to estimate
#   the diagonal parameters.
# 2010-10-08 [HB]
# o Added estimateC2Bias().
# o Added fitDeltaXYShearModel().
# o Now deShearC1C2() uses fitDeltaC1C2ShearModel().
# o Added fitDeltaC1C2ShearModel().
# o Now deShearC1C2() returns the 'modelFit'.
# o Now deShearC1C2() calls peaks using callPeaks() for PeaksAndValleys.
# 2010-09-26 [HB]
# o Added argument 'adjust' to deShearC1C2() with new default.
# o Added sanity checks to deShearC1C2().
# o Now normalizeBAFsByRegions() for PairedPSCBS handles multiple chromosomes.
# 2010-09-22 [PN]
# o Added deShearC1C2() for PairedPSCBS.
# 2010-09-19 [HB+PN]
# o Added orthogonalizeC1C2() for PairedPSCBS.
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
