setConstructorS3("TotalFracBSnpData", function(data=NULL, ...) {
  if (!is.null(data))
    colnames(data) <- c("total", "fracB");
  extend(SnpData(data, ...), "TotalFracBSnpData");
})

setMethodS3("plot", "TotalFracBSnpData", function(x, xlim=NULL, ylim=c(0,1), ...) {
  NextMethod("plot", xlim=xlim, ylim=ylim);
})


setMethodS3("asTotalFracBSnpData", "TotalFracBSnpData", function(this, ...) {
  this;
})


setMethodS3("asPolarSnpData", "TotalFracBSnpData", function(this, ...) {
  theta <- asCartesianSnpData(this);
  asPolarSnpData(theta, ...);
})

setMethodS3("asCartesianSnpData", "TotalFracBSnpData", function(this,
...) {
  totalFracB <- this;

  data <- this;
  # thetaB
  data[,2] <- totalFracB[,2]*totalFracB[,1];
  # thetaA
  data[,1] <- totalFracB[,1]-data[,2];

  CartesianSnpData(data, ...);
})


setMethodS3("callGenotypes", "TotalFracBSnpData", function(this, adjust=1.5, meanCN=2, ...) {
  data <- this;

  # Argument 'meanCN':
  meanCN <- Arguments$getInteger(meanCN, range=c(1,2));

  nbrOfUnits <- nrow(data);

  beta <- data[,"fracB"];
  fit <- .findPeaksAndValleys(beta, adjust=adjust);
  type <- NULL; rm(list="type"); # To please R CMD check
  fit <- subset(fit, type == "valley");
  nbrOfGenotypeGroups <- nrow(fit)+1;
  if (nbrOfGenotypeGroups == 1) {
    warning("PRECISION ERROR: Only one genotype group was detected.");
    if (meanCN == 1) {
      mu <- rep(1/2, times=nbrOfUnits);
      a <- 1/2;
      mu[beta < a] <- 0;
      mu[beta > a] <- 1;
    } else if (meanCN == 2) {
      a <- 1/3;
      b <- 2/3;
      mu[beta < a] <- 0;
      mu[beta > b] <- 1;
    }
  } else if (nbrOfGenotypeGroups == 2) {
    a <- fit$x[1L];
    mu <- rep(0, times=nbrOfUnits);
    mu[beta > a] <- 1;
  } else if (nbrOfGenotypeGroups == 3) {
    a <- fit$x[1L];
    b <- fit$x[2L];
    mu <- rep(1/2, times=nbrOfUnits);
    mu[beta < a] <- 0;
    mu[beta > b] <- 1;
  } else {
    throw("Unexpected number of genotype groups: ", nbrOfGenotypeGroups);
  }
  mu[is.na(beta)] <- NA;

  res <- this;
  res[,2] <- mu;
  attr(res, "fit") <- fit;

  res;
})



setMethodS3("pairedBoost", "TotalFracBSnpData", function(this, dataN, ...) {
  # Coerce normal SNP signals
  tfN <- asTotalFracBSnpData(dataN);

  # Assert compatibility
  if (!identical(dim(this), dim(tfN))) {
    throw("Argument 'dataN' is of a non-compatible dimension.");
  }

  # Call normal genotypes
  tfNC <- callGenotypes(tfN, ...);

  # Estimate beta correction factors
  delta <- tfN[,2L]-tfNC[,2L];

  # Calibrate accordingly
  res <- this;
  res[,2] <- res[,2L] - delta;

  res;
})


############################################################################
# HISTORY:
# 2009-04-28
# o Added argument 'meanCN=2' to callGenotypes().  It used in order to
#   identify the modes in case the expect CN == 1, e.g. Chr X or Chr Y.
# 2009-03-31
# o Added pairedBoost() for TotalFracBSnpData.
# 2009-03-30
# o Created.
############################################################################
