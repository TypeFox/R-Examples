fitCalMaTeMedians <- function(dataT, references, fB1=1/3, fB2=2/3, ...) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.
  nbrOfSNPs <- nrow(dataT);
  nbrOfReferences <- length(references);

  # Adding a small value so there are "non" 0 values
  eps <- 1e-6;
  dataT[dataT < eps] <- eps;

  eps2 <- 1e-4;  
  a <- max(max(dataT[2,references] / (pmax(dataT[1,references],0) + eps2)), max(dataT[1,references] / (pmax(dataT[2,references],0) + eps2)));
  Giro <- matrix(c(1, 1/a, 1/a, 1), nrow=2, ncol=2, byrow=FALSE);
  Giro <- solve(Giro);
  dataT <- Giro %*% dataT;

  # Extract the signals for the reference set
  TR <- dataT[,references, drop=FALSE];
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Genotyping
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fracB <- TR[2,] / (TR[1,] + TR[2,]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Group the three possibilities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  wgths <- c(sum(naiveGenoDiff == +2),
             sum(naiveGenoDiff ==  0),
             sum(naiveGenoDiff == -2))
  
  TR <- cbind(rowMedians(TR[,naiveGenoDiff ==+2, drop=FALSE]),
              rowMedians(TR[,naiveGenoDiff == 0, drop=FALSE]),
              rowMedians(TR[,naiveGenoDiff ==-2, drop=FALSE]));
  # Remove possible missing values
  TR[is.nan(TR)] <- 0;
  
  # Handle cases a single genotype is present in the samples

  # A very unlikely case: all the samples are heterozygous
  if ((wgths[1] == 0) && (wgths[3] == 0)) {
    TR[1,1] <- 2.0*TR[1,2];
    TR[2,1] <- 0.1*TR[2,2];
    wgths[1] <- wgths[2];
  }
  # Only BB
  if ((wgths[1] == 0) && (wgths[2] == 0)) {
    TR[,1] <- TR[2:1,3];
    wgths[1] <- wgths[3];
  }
  if ((wgths[3] == 0) && (wgths[2] == 0)) {
    TR[,3] <- TR[2:1,1];
    wgths[3] <- wgths[1];
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibration
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Genotypes <- cbind(c(2,0), c(1,1), c(0,2));
  P <- qr.solve(t(TR) * wgths^2, t(Genotypes) * wgths^2);
  res <- t(P) %*% dataT;
  
  dataT <- res;

  TR <- dataT[,references, drop=FALSE];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Have the genotypes changed?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fracB <- TR[2,] / (TR[1,] + TR[2,]);
  naiveGenoDiff_1 <- 2*(fracB < fB1) - 2*(fracB > fB2);
  
  if (!identical(naiveGenoDiff_1, naiveGenoDiff)) {
    # Repeat the process once
    naiveGenoDiff <- naiveGenoDiff_1;
    wgths <- c(sum(naiveGenoDiff == +2),
               sum(naiveGenoDiff ==  0),
               sum(naiveGenoDiff == -2))

    TR <- cbind(rowMedians(TR[,naiveGenoDiff == +2, drop=FALSE]),
                rowMedians(TR[,naiveGenoDiff ==  0, drop=FALSE]),
                rowMedians(TR[,naiveGenoDiff == -2, drop=FALSE]));
    # Remove possible missing values
    TR[is.nan(TR)] <- 0;

    # Handle cases in which not all the genotypes are present in the samples

    # A very unlikely case: all the samples are heterozygous
    if ((wgths[1] == 0) && (wgths[3] == 0)) {
      TR[1,1] <- 2.0*TR[1,2];
      TR[2,1] <- 0.1*TR[2,2];
      wgths[1] <- wgths[2];
    }
    # Only BB
    if ((wgths[1] == 0) && (wgths[2] == 0)) {
      TR[,1] <- TR[2:1,3];
      wgths[1] <- wgths[3];
    }
    if ((wgths[3] == 0) && (wgths[2] == 0)) {
      TR[,3] <- TR[2:1,1];
      wgths[3] <- wgths[1];
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Callibration
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    Genotypes <- cbind(c(2,0), c(1,1), c(0,2));
    P <- qr.solve(t(TR) * wgths^2, t(Genotypes) * wgths^2);
    res <- t(P) %*% dataT;
  }

  res;
} # fitCalMaTeMedians()

###########################################################################
# HISTORY:
# 2012-02-21 [HB]
# o CODE HARMONIZATION: Renamed argument 'T' to 'dataT', cf.
#   fitCalMaTeV1() and fitCalMaTeV2().
# o SPEEDUP: Made fitCalMaTeMedians() a plain function to avoid the
#   method-dispatch overhead, cf. fitCalMaTeV1() and fitCalMaTeV2().
#   DOCUMENTATION: Moved documentation to fitCalMaTeInternal.R.
# o CODE CLEANUP: Cleaned up spacing etc. Dropped unused 'data' variable.
# o CODE ROBUSTNESS: Always use 'if (x && y)' and never 'if (x & y)'.
# 2011-04-12 [AR] [Is date correct? /HB 2012-02-21]
# o Created from fitCalMaTe().
#   - The previous version used rlm(). This one does not.
#   - Robustness is obtained by using medians among different genotypes
#     and weighting according to the number or samples in each genotype.
###########################################################################
