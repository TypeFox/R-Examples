# \references{
#   [1] P. Neuvial, \emph{Improved CN estimation and detection for SNP arrays}, private notes, 2009-12-01.\cr
# }
##
## NOMENCLATURE?!?
##
## ascnT [Jx2], ascnN [Jx2]  => (thetaA, thetaB)
## 
## pscnT [Jx2], pscnN [Jx2]  => (theta, beta)
## 
## pscnT [Jx2], pscnN [Jx2]  => (theta1, theta2)
##
## data <- calculatePairedPSCNByGenotype(pscnT, pscnN);
## fit <- segmentByPairedPSCBS(data, tbn=FALSE);
## fit <- callROH(fit);
##
setMethodS3("calculatePairedPSCNByGenotype", "numeric", function(thetaT, betaT, thetaN, betaN, muN=callNaiveGenotypes(betaN, ...), flavor=c("v1"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'thetaT':
  thetaT <- as.numeric(thetaT);

  # Argument 'betaT':
  betaT <- as.numeric(betaT);
  J <- length(thetaT);
  if (length(betaT) != J) {
    stop("The length of arguments 'betaT' and 'thetaT' differ: ", 
                                                   length(betaT), " != ", J);
  }

  # Argument 'thetaN':
  thetaN <- as.numeric(thetaN);
  if (length(thetaN) != J) {
    stop("The length of arguments 'thetaN' and 'thetaT' differ: ", 
                                                   length(thetaN), " != ", J);
  }

  # Argument 'betaN':
  betaN <- as.numeric(betaN);
  if (length(betaN) != J) {
    stop("The length of arguments 'betaN' and 'thetaT' differ: ", 
                                                   length(betaN), " != ", J);
  }
  
  # Argument 'muN':
  if (is.null(muN)) {
    muN <- callNaiveGenotypes(betaN, ...);
  }
  if (length(muN) != J) {
    stop("The length of arguments 'muN' and 'thetaT' differ: ", 
                                                   length(muN), " != ", J);
  }
  knownGenotypes <- c(0, 1/2, 1, NA);
  unknown <- which(!is.element(muN, knownGenotypes));
  n <- length(unknown);
  if (n > 0) {
    unknown <- unique(muN[unknown]);
    unknownStr <- paste(unknown, collapse=", ");
    stop("Argument 'muN' contains unknown values: ", unknownStr);
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate CT for all loci [and record (betaT,muN)]
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  CTx <- 2 * thetaT/thetaN;
  betaTx <- betaT;
  muNx <- muN;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate (CT,betaT) for SNPs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify set to be updated (our best what are SNPs)
  isSNP <- which(is.finite(betaN) & is.finite(muN));

  # Subset
  thetaT <- thetaT[isSNP];
  betaT <- betaT[isSNP];
  thetaN <- thetaN[isSNP];
  betaN <- betaN[isSNP];
  muN <- muN[isSNP];

  # Transform to (thetaA, thetaB) for tumors and normals (only for SNPs)
  thetaTB <- thetaT * betaT;
  thetaTA <- thetaT - thetaTB;

  thetaNB <- thetaN * betaN;
  thetaNA <- thetaN - thetaNB;

  # Identify genotypes
  isAA <- (muN == 0);
  isAB <- (muN == 1/2);
  isBB <- (muN == 1);

  # The true allele-specific CNs in the normal (CNA,CNB)
  CNA <- 2*isAA + 1*isAB;
  CNB <- 2*isBB + 1*isAB;

  # Calculate CT stratified by genotype  
  CTA <- thetaTA/thetaNA*CNA;
  CTB <- thetaTB/thetaNB*CNB;

  CT <- CTA + CTB;
  betaT <- CTB / CT;

  # Sanity checks
  stopifnot(all(betaT[isAA] == 0));
  stopifnot(all(betaT[isBB] == 1));
  stopifnot(all(CTA[isBB] == 0));
  stopifnot(all(CTB[isAA] == 0));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Expand to all loci
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update SNPs
  CTx[isSNP] <- CT;
  betaTx[isSNP] <- betaT;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return calculate data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  stopifnot(length(CTx) == J);
  stopifnot(length(betaTx) == J);
  stopifnot(length(muNx) == J);

  list(CT=CTx, betaT=betaTx, muN=muNx);
}) # calculatePairedPSCNByGenotype()


setMethodS3("calculatePairedPSCNByGenotype", "data.frame", function(data, ...) {
  res <- calculatePairedPSCNByGenotype(thetaT=data$thetaT, betaT=data$betaT, thetaN=data$thetaN, betaN=data$betaN, muN=data$muN, ...);
  data$CT <- res$CT;
  data$betaT <- res$betaT;
  data$muN <- res$muN;
  data;
})



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# TO BE RETIRED
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# \references{
#   [1] P. Neuvial, \emph{Improved CN estimation and detection for SNP arrays}, private notes, 2009-12-01.\cr
# }
setMethodS3("calculateTumorPSCNByGenotypeUsingCTCN", "numeric", function(CT, betaT, betaN, muN=callNaiveGenotypes(betaN, ...), flavor=c("v1"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'CT':
  CT <- as.numeric(CT);

  # Argument 'betaT':
  betaT <- as.numeric(betaT);
  J <- length(CT);
  if (length(betaT) != J) {
    stop("The length of arguments 'betaT' and 'CT' differ: ", 
                                                   length(betaT), " != ", J);
  }

  # Argument 'betaN':
  betaN <- as.numeric(betaN);
  if (length(betaN) != J) {
    stop("The length of arguments 'betaN' and 'CT' differ: ", 
                                                   length(betaN), " != ", J);
  }
  
  # Argument 'muN':
  if (is.null(muN)) {
    muN <- callNaiveGenotypes(betaN, ...);
  }
  if (length(muN) != J) {
    stop("The length of arguments 'muN' and 'CT' differ: ", 
                                                   length(muN), " != ", J);
  }
  knownGenotypes <- c(0, 1/2, 1, NA);
  unknown <- which(!is.element(muN, knownGenotypes));
  n <- length(unknown);
  if (n > 0) {
    unknown <- unique(muN[unknown]);
    unknownStr <- paste(unknown, collapse=", ");
    stop("Argument 'muN' contains unknown values: ", unknownStr);
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate CT for all loci [and record (betaT,muN)]
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  CTx <- CT;
  betaTx <- betaT;
  muNx <- muN;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate (CT,betaT) for SNPs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify set to be updated (our best what are SNPs)
  isSNP <- which(is.finite(betaN) & is.finite(muN));

  # Subset
  CT <- CT[isSNP];
  betaT <- betaT[isSNP];
  betaN <- betaN[isSNP];
  muN <- muN[isSNP];

  # Identify genotypes
  isAA <- (muN == 0);
  isAB <- (muN == 1/2);
  isBB <- (muN == 1);

  # True (CNA,CNB)
  CNA <- 2*isAA + 1*isAB;
  CNB <- 2*isBB + 1*isAB;

  # Calculate (CTA,CTB) stratified by genotype
  alphaT <- 1 - betaT;
  alphaN <- 1 - betaN;
  kappaA <- (alphaT/alphaN) * CT/2;
  kappaB <- ( betaT/ betaN) * CT/2;
  CTA <- kappaA * CNA;
  CTB <- kappaB * CNB;

  CT <- CTA + CTB;
  betaT <- CTB / CT;

  # Sanity checks
  stopifnot(all(betaT[isAA] == 0, na.rm=TRUE));
  stopifnot(all(betaT[isBB] == 1, na.rm=TRUE));
  stopifnot(all(CTA[isBB] == 0, na.rm=TRUE));
  stopifnot(all(CTB[isAA] == 0, na.rm=TRUE));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Expand to all loci
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update SNPs
  CTx[isSNP] <- CT;
  betaTx[isSNP] <- betaT;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return calculate data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  stopifnot(length(CTx) == J);
  stopifnot(length(betaTx) == J);
  stopifnot(length(muNx) == J);

  list(CT=CTx, betaT=betaTx, muN=muNx);
}) # calculateTumorPSCNByGenotypeUsingCTCN()



############################################################################
# HISTORY:
# 2012-09-16
# o Renamed to calculatePairedPSCNByGenotype().
# o Workname: "SACKS"? (from Sack's coffee in Berkeley)
# o Added calculateTumorPSCNByGenotypeUsingTheta().
# o Added calculateTumorPSCNByGenotype().
# o Created.
# 2012-09-14 [two strange looking scientist hanging out at the Cal Fountain
#   across from Caffe Strada just before sunset]
# o Decided to go back to this tumor PSCN definition from 2009, which was 
#   discussed while trying to come up with an alternative interpretation
#   of that TumorBoost normalization is doing.  See PN's extensive
#   "private" report [1] on this.
############################################################################ 
