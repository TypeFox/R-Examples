setConstructorS3("CartesianSnpData", function(data=NULL, ...) {
  if (!is.null(data))
    colnames(data) <- c("A", "B");
  extend(SnpData(data, ...), "CartesianSnpData");
})


setMethodS3("plot", "CartesianSnpData", function(x, xlim=range(x, na.rm=TRUE), ylim=xlim, ...) {
  NextMethod("plot", xlim=xlim, ylim=ylim);
})


setMethodS3("asCartesianSnpData", "CartesianSnpData", function(this, ...) {
  this;
})

setMethodS3("asPolarSnpData", "CartesianSnpData", function(this, ...) {
  theta <- this;

  data <- this;
  # Radius
  data[,1] <- sqrt(theta[,1]^2 + theta[,2]^2);
  # Azimuth
  data[,2] <- asin(theta[,2]/data[,1]);

  PolarSnpData(data, ...);
})

setMethodS3("asTotalFracBSnpData", "CartesianSnpData", function(this, ...) {
  theta <- this;

  data <- theta;
  # Total
  data[,1] <- theta[,"A"]+theta[,"B"];
  # FracB
  data[,2] <- theta[,"B"]/data[,1];

  TotalFracBSnpData(data, ...);
})


setMethodS3("pairedBoost", "CartesianSnpData", function(this, dataN, scaleByCN=TRUE, flavor=c("diagonal", "manhattan"), ...) {
##setMethodS3("pairedBoost", "CartesianSnpData", function(this, dataN, scaleByCN=TRUE, flavor=c("manhattan", "diagonal"), ...) {
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Coerce normal SNP signals
  thetaN <- asCartesianSnpData(dataN);

  # Assert compatibility
  if (!identical(dim(this), dim(thetaN))) {
    throw("Argument 'dataN' is of a non-compatible dimension.");
  }

  # Call normal genotypes
  thetaNC <- callGenotypes(thetaN, ...);

  isAA <- (thetaNC[,1] != 0 & thetaNC[,2] == 0);
  isAB <- (thetaNC[,1] != 0 & thetaNC[,2] != 0);
  isBB <- (thetaNC[,1] == 0 & thetaNC[,2] != 0);

  # Estimate Cartesian correction factors
  delta <- array(0, dim=dim(thetaN));

  # AA:s
  idxs <- which(isAA);
  delta[idxs,2] <- thetaN[idxs,2]-thetaNC[idxs,2];
  if (flavor == "diagonal") {
    delta[idxs,1] <- -delta[idxs,2];
  }

  # BB:s
  idxs <- which(isBB);
  delta[idxs,1] <- thetaN[idxs,1]-thetaNC[idxs,1];
  if (flavor == "diagonal") {
    delta[idxs,2] <- -delta[idxs,1];
  }

  # AB:s
  idxs <- which(isAB);
  delta[idxs,] <- thetaN[idxs,]-thetaNC[idxs,];

  # Scale by CN?
  if (scaleByCN) {
    # Estimate relative CNs
    thetaT <- this;
    totalT <- thetaT[,1]+thetaT[,2];
    totalN <- thetaN[,1]+thetaN[,2];
    C <- totalT/totalN;
    delta <- C*delta;
  }

  # Calibrate accordingly
  res <- this;
  res <- res - delta;

  res;
})


############################################################################
# HISTORY:
# 2009-03-31
# o Added pairedBoost() for CartesianSnpData.
# 2009-03-30
# o Created.
############################################################################
