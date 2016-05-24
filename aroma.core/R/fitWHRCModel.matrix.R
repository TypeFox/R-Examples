# Fit the log-additive model assuming *homoscedastic* error terms.
# Each data element can be given a weight.  Moreover, if there
# are missing values, these will be given zero weights.
setMethodS3("fitWRCModel", "matrix", function(y, ...) {
  # Explicit call to avoid method dispatching overheads.
  fitWHRCModel.matrix(y, ..., maxIterations=1);
})




# Fit the log-additive model assuming *heteroscedastic* error terms
# such that *each probe* may have a different standard deviation.
# Each data element can be given a weight.  Moreover, if there
# are missing values, these will be given zero weights.
# 'tau' is an optional penalty term to avoid zero variance estimates.
setMethodS3("fitWHRCModel", "matrix", function(y, w=NULL, hasNAs=TRUE, tau=1e-3, eps=1e-3, maxIterations=100, .checkArgs=TRUE, ..., .loadDeps=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Constants
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  PACKAGE <- "preprocessCore";
  if (.loadDeps) {
    require(PACKAGE, character.only=TRUE) || throw("Package not loaded: ", PACKAGE);
  }
  rcModelWPLM <- NULL; rm(list="rcModelWPLM"); # To please R CMD check

  I <- ncol(y);  # Number of arrays
  K <- nrow(y);  # Number of probes

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.checkArgs) {
    # Argument 'w':
    if (!is.null(w)) {
      if (!identical(dim(w), c(I,K)))
        throw(sprintf("The dimension of argument 'w' does not match that of 'y': %dx%d != %dx%d", nrow(w), ncol(w), I, K));
    }
  }


  # Scan for NAs?
  if (hasNAs) {
    isNA <- is.na(y);
  } else if (is.na(hasNAs)) {
    isNA <- is.na(y);
    hasNAs <- any(isNA);
  }

  if (hasNAs) {
    # Weight matrix where missing values gets zero weight.
    wNA <- matrix(1, nrow=K, ncol=I);
    wNA[isNA] <- 0;
    y[isNA] <- 0; # preprocessCore() only takes real values.

    badCells <- colAlls(isNA);
    if (any(badCells)) {
      return(list(Estimates=rep(NA_real_, times=I+K),
                  StdErrors=rep(NA_real_, times=I+K),
                  Weights=matrix(NA_real_, nrow=K, ncol=I),
                  Residuals=matrix(NA_real_, nrow=K, ncol=I)
                 )
            );
    }

    # Update prior weights?
    if (is.null(w)) {
      w <- wNA;
    } else {
      w <- w * wNA;
    }
  } else {
    w <- matrix(1, nrow=K, ncol=I);
  }

  thetaIdxs <- seq_len(I);
  phiIdxs <- I+seq_len(K);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Iterative re-weighted fit
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  w0 <- w;
  dParams <- NA_real_;
  lastParams <- NULL;

  nbrOfIterations <- 0;
  ready <- FALSE;
  while (!ready) {
    fit <- rcModelWPLM(y, w=w);

    nbrOfIterations <- nbrOfIterations + as.integer(1);

    params <- fit$Estimates;

    if (is.null(lastParams)) {
      hasConverged <- FALSE;
    } else {
      dParams <- mean(abs(params-lastParams), na.rm=TRUE);
      hasConverged <- (dParams < eps);
    }
    ready <- (hasConverged || nbrOfIterations >= maxIterations);

    if (!ready) {
      # Estimate probe standard deviations robustly
      phi <- params[phiIdxs];
      res <- fit$Residuals;
      res[isNA] <- NA;  # Missing values were set to zero.
      phiSd <- rowMads(res, center=0, na.rm=TRUE);

      # Downweights by variance (with a small penalty term)
      wV <- 1/(phiSd^2 + tau);
#      wV <- wV / sum(wV, na.rm=TRUE);

      # Update the weights
      w <- wV * w0;
    }

    lastParams <- params;
  } # while(!ready)

  if (maxIterations > 1) {
    fit$hasConverged <- hasConverged;
    fit$dParams <- dParams;
    fit$nbrOfIterations <- nbrOfIterations;
  }

  fit;
}) # fitWHRCModel()


############################################################################
# HISTORY:
# 2012-08-21
# o ROBUSTNESS: fitWHRCModel() now calls rcModelWPLM() of preprocessCore
#   instead of internal .Call(..., PACKAGE="preprocessCore") calls.
# o CLEANUP: Dropped arguments 'psiCode' and 'psiK' from fitWHRCModel().
# 2012-08-08
# o Now fitWHRCModel() allocates numerical NAs (instead of logical ones).
# o Added argument '.loadDeps' to fitWRMA().
# o Now R CMD check is no longer complaining about .Call(..., PACKAGE).
# 2007-10-05
# o Created.
############################################################################
