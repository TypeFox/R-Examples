setConstructorS3("PolarSnpData", function(data=NULL, ...) {
  if (!is.null(data))
    colnames(data) <- c("radius", "azimuth");
  extend(SnpData(data, ...), "PolarSnpData");
})


setMethodS3("plot", "PolarSnpData", function(x, xlim=NULL, ylim=c(0,pi/2), ...) {
  NextMethod("plot", xlim=xlim, ylim=ylim);
})


setMethodS3("asPolarSnpData", "PolarSnpData", function(this, ...) {
  this;
})

setMethodS3("asCartesianSnpData", "PolarSnpData", function(this, ...) {
  polar <- this;

  theta <- this;
  theta[,1] <- polar[,1] * cos(polar[,2]);
  theta[,2] <- polar[,1] * sin(polar[,2]);

  CartesianSnpData(theta, ...);
})


setMethodS3("asTotalFracBSnpData", "PolarSnpData", function(this, ...) {
  theta <- asCartesianSnpData(this);
  asTotalFracBSnpData(theta, ...);
})


setMethodS3("pairedBoost", "PolarSnpData", function(this, dataN, ...) {
  # Coerce normal SNP signals
  polarN <- asPolarSnpData(dataN);

  # Assert compatibility
  if (!identical(dim(this), dim(polarN))) {
    throw("Argument 'dataN' is of a non-compatible dimension.");
  }

  # Call normal genotypes
  polarNC <- callGenotypes(polarN, ...);

  # Estimate beta correction factors
  delta <- polarN[,2]-polarNC[,2];

  # Calibrate accordingly
  res <- this;
  res[,2] <- res[,2] - delta;
  
  res;
})



############################################################################
# HISTORY:
# 2009-03-31
# o Added pairedBoost() for PolarSnpData.
# 2009-03-30
# o Created.
############################################################################
