# extractTheta() for SnpQSet and SnpCnvQSet. Currently,
# it is not possible to subset by unit indices.
setMethodS3("extractTheta", "SnpQSet", function(this, ..., transform=function(y, ...) { 2^y }, addNames=TRUE, verbose=FALSE) {
  eSet <- this;

  # To please R CMD check
  ns <- getNamespace("oligo");
  if (!exists("senseThetaA", mode="function", envir=ns)) {
    throw("This methods is only supported for older versions of the 'oligo' package: oligo v", packageVersion("oligo"));
  }
  senseThetaA <- get("senseThetaA", mode="function", envir=ns);
  senseThetaB <- get("senseThetaB", mode="function", envir=ns);
  antisenseThetaA <- get("antisenseThetaA", mode="function", envir=ns);
  antisenseThetaB <- get("antisenseThetaB", mode="function", envir=ns);


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
  nbrOfGroups <- 4L;  # (thetaA+, thetaB+, thetaA-, thetaB-);

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
  theta[,1L,] <- transform(senseThetaA(eSet));
  theta[,2L,] <- transform(senseThetaB(eSet));
  theta[,3L,] <- transform(antisenseThetaA(eSet));
  theta[,4L,] <- transform(antisenseThetaB(eSet));

  theta;
})


############################################################################
# HISTORY:
# 2013-10-07
# o ROBUSTNESS: Now extractTheta() for SnpQSet gives an informative error
#   message that it only applies with older version of 'oligo'.
# 2008-12-05
# o Created.
############################################################################
