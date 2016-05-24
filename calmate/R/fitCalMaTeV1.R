fitCalMaTeV1 <- function(dataT, references, fB1=1/3, fB2=2/3, maxIter=50, ...) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.
  nbrOfSNPs <- nrow(dataT);
  nbrOfReferences <- length(references);
  
  # Adding a small value so there are "non" 0 values
  eps <- 1e-6;
  dataT[dataT < eps] <- eps;
  
  eps2 <- 1e-4;
  a <- max(max(dataT[2,] / (pmax(dataT[1,],0) + eps2)), max(dataT[1,] / (pmax(dataT[2,],0) + eps2)));
  Giro <- matrix(c(1, 1/a, 1/a, 1), nrow=2, ncol=2, byrow=FALSE);
  Giro <- solve(Giro);
  dataT <- Giro %*% dataT;

  # Extract the signals for the reference set
  TR <- dataT[,references, drop=FALSE];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Checking if all the samples are homozygous
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fracB <- TR[2,] / (TR[1,] + TR[2,]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Twist half of the reference samples in case there is only one allele?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  onlyOneAllele <- (abs(sum(naiveGenoDiff)/2) == length(naiveGenoDiff));
  if (onlyOneAllele) {
    idxsSwap <- references[seq(length=ncol(TR)/2)];
    dataT[1:2,idxsSwap] <- dataT[2:1,idxsSwap, drop=FALSE];

    # Update precalcalculated signals
    TR <- dataT[,references, drop=FALSE];
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Total copy numbers must be close to 2 for the reference samples or
  # (if there are not control samples) for most of the samples
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  H <- matrix(2, nrow=nbrOfReferences, ncol=1, byrow=FALSE);
  fit <- rlm(t(TR), H, maxit=maxIter);
  matSum <- fit$coefficients;
  coeffs <- fit$w;
  dataT <- diag(matSum) %*% dataT;

  # Reextract the signals for the reference set
  TR <- dataT[,references, drop=FALSE];

  # The difference of the copy numbers must be 2, 0 or -2 depending genotyping
  fracB <- TR[2,] / (TR[1,] + TR[2,]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
  fit <- rlm(t(TR), naiveGenoDiff, maxit=maxIter, weights=coeffs);
  matDiff <- fit$coefficients;

  # T matrix is:
  #  [1  1] [   ] = [MatSum[1]   MatSum[2]] (We have already applied it) MatSum is 1,1
  #  [1 -1] [ T ]   [MatDiff[1] MatDiff[2]]
  U <- matrix(c(0.5, 0.5, 0.5, -0.5), nrow=2, ncol=2, byrow=FALSE);
  V <- matrix(c(c(1, 1), matDiff), nrow=2, ncol=2, byrow=TRUE);
  T <- U %*% V;
  
  res <- T %*% dataT;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Undo the previous change applied to the data in case there is 
  # only one allele    
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (onlyOneAllele) {
    res[1:2,idxsSwap] <- res[2:1,idxsSwap, drop=FALSE];
  }

  # Return parameter estimates(?)
  ## attr(res, "modelFit") <- list(fit=fit);

  res;
} # fitCalMaTeV1()

###########################################################################
# HISTORY:
# 2012-02-19 [HB]
# o Clarified in the source code comments that it is only the reference
#   samples that are "twisted".
# o Created internal fit functions for the different versions of CalMaTe.
# 2012-01-31 [MO]
# o BUG FIX: the index "idxs" was recalculated to undo the change when 
#   there is only one allele, and it was done as the previous version,
#   taking into account all the samples, not only the references.     
# 2011-11-29 [MO]
# o Change matrix "T" by "dataT" and "P" by "T"
# 2010-08-02 [HB]
# o ROBUSTNESS: Now fitCalMaTe() also works (technically) when there is
#   only one reference.
# o Made into an S3 method for matrix:es.
# 2010-06-18 [HB]
# o Created from refineCN.list().
###########################################################################
