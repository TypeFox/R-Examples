setMethodS3("extractTheta", "SnpCnvQSet", function(this, ..., transform=function(y, ...) { 2^y }, addNames=TRUE, verbose=FALSE) {
  eSet <- this;

  # To please R CMD check
  ns <- loadNamespace("oligo");
  thetaA <- get("thetaA", mode="function", envir=ns);
  thetaB <- get("thetaB", mode="function", envir=ns);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  nbrOfUnits <- nrow(eSet);
  nbrOfSamples <- ncol(eSet);
  nbrOfGroups <- 2L;  # (thetaA, thetaB)

  # Extract sample names
  sampleNames <- .sampleNames(eSet);
  sampleNames <- gsub("[.](cel|CEL)$", "", sampleNames);
  sampleNames <- gsub(",.*$", "", sampleNames);

  # Extract unit names
  snpNames <- .featureNames(eSet);

  # Allocate result object
  naValue <- as.double(NA);
  theta <- array(naValue, dim=c(nbrOfUnits, nbrOfGroups, nbrOfSamples));
  dimnames(theta)[[3L]] <- sampleNames;
  if (addNames)
    dimnames(theta)[[1L]] <- snpNames;

  # Populate with estimates
  theta[,1L,] <- transform(thetaA(eSet));
  theta[,2L,] <- transform(thetaB(eSet));

  theta;
})


############################################################################
# HISTORY:
# 2008-12-05
# o Created.
############################################################################
