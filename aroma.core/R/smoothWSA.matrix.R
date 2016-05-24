setMethodS3("smoothWSA", "matrix", function(Y, x, w=NULL, kernel=gaussKernel, sd=100e3, na.rm=TRUE, ..., progress=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'Y'
  K <- nrow(Y);  # Number of positions
  I <- ncol(Y);  # Number of samples

  # Argument 'x'
  if (length(x) != K) {
    throw("Argument 'x' has different number of values that rows in 'Y': ",
                                                     length(x), " != ", K);
  }

  # Argument 'w'
  if (is.null(w)) {
    w <- 1; # Uniform prior weights.
  } else if (is.matrix(w)) {
    if (nrow(w) != K) {
      throw("Argument 'w' has different number of rows than 'Y': ",
                                                     nrow(w), " != ", K);
    }
    if (ncol(w) != I) {
      throw("Argument 'w' has different number of columns than 'Y': ",
                                                     ncol(w), " != ", I);
    }
  } else if (is.vector(w)) {
    if (length(w) != K) {
      throw("Argument 'w' has different number of values that rows in 'Y': ",
                                                     length(w), " != ", K);
    }
  }
  if (any(w < 0))
    throw("Argument 'w' contains negative weights.");
  if (any(!is.finite(w)))
    throw("Argument 'w' contains non-finite weights.");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup up
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit everything on the log-ratio scale
  YR <- rowMedians(Y, na.rm=TRUE);
  M <- log2(Y) - log2(YR);

  # Identify missing values?
  if (na.rm) {
    nas <- is.na(M);
    dim(nas) <- c(K,I);
  }

  # Allocate vector of smoothed signals
  naValue <- as.double(NA);
  theta <- matrix(naValue, nrow=K, ncol=I);
  phi <- rep(naValue, times=K);

  # At each position, calculate the weighed average using a
  # Gaussian kernel.
  cat("Progress: ");
  for (kk in seq_len(K)) {
    if (progress && kk %% 100 == 0)
      cat(kk, ", ", sep="");

    # Weights centered around x[kk]
    wK <- kernel(x, mean=x[kk], sd=sd);

    # Weight matrix
    wM <- matrix(wK, nrow=K, ncol=I);
    # Not needed anymore
    wK <- NULL;

    # Multiple with prior (row) weights
    wM <- w*wM;

    # Give missing values zero weight?
    if (na.rm)
      wM[nas] <- 0;

    wMR <- rowSums(wM);
    keep <- which(wMR > 0);
    # Not needed anymore
    wMR <- NULL;
    if (length(keep) > 0) {
      wM <- wM[keep,,drop=FALSE];
      m <- M[keep,,drop=FALSE];
      verbose && print(verbose, list(m=m, wM=wM));

      # Standardize wM such that each column sum to one.
      wMs <- colSums(wM, na.rm=TRUE);
      wM <- wM/wMs;

      # Weighted average
      theta[kk,] <- colSums(wM*m);

      # Not needed anymore
      m <- wMs <- NULL;
    } else {
      # If no data, keep theta:s and phi:s as (allocated) missing values.
    }
    # Not needed anymore
    wM <- wMR <- keep <- NULL;
  }
  cat(kk, "\n", sep="");

  # Above, theta holds smoothed M = log2(Y) - log2(YR).
  theta <- theta + log2(YR);

  # Return everything on the intensity scale
  theta <- 2^theta;
  phi <- 2^phi;

  list(theta=theta, phi=phi);
}) # smoothWSA()


############################################################################
# HISTORY:
# 2012-08-08
# o Now smoothWSA() allocates with numerical NAs (instead of logical ones).
# 2007-09-26
# o Added support for (optional) prior weights (either as row weights or
#   full matrix weights).
# 2007-09-18
# o Created from gaussianSmoothing.R.
############################################################################
