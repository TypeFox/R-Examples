# Fit the log-additive model assuming *homoscedastic* error terms.
# Each data element can be given a weight.  Moreover, if there
# are missing values, these will be given zero weights.
setMethodS3("fitWLAPLM", "matrix", function(y, ...) {
  # Explicit call to avoid method dispatching overheads.
  fitWLAPLM.matrix(y, ..., maxIterations=1);
})




setMethodS3("fitWHLAPLM", "matrix", function(y, ...) {
  K <- nrow(y);
  I <- ncol(y);
  y <- log2(y);
  thetaIdxs <- seq_len(I);
  phiIdxs <- I+seq_len(K);

  fit <- fitWHRCModel.matrix(y, ...);

  est <- fit$Estimates;
  se <- fit$StdErrors;

  # Chip effects
  thetaIdxs <- 1:I;
  beta <- est[thetaIdxs];
  theta <- 2^beta;

  # Probe affinities
  phiIdxs <- (I+1):(I+K);
  alpha <- est[phiIdxs];
  alpha[K] <- -sum(alpha[1:(K-1)]);
  phi <- 2^alpha;

  # The RMA model is already fitted with constraint prod(phi) = 1.
  # No rescaling needed.

  seTheta <- 2^(se[thetaIdxs]);
  sePhi <- 2^(se[phiIdxs]);

  fit$theta <- theta;
  fit$seTheta <- seTheta;
  fit$phi <- phi;
  fit$sePhi <- sePhi;

  fit;
}) # fitWHLAPLM()


############################################################################
# HISTORY:
# 2007-10-05
# o Created.
############################################################################
