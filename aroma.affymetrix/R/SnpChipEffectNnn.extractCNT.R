setMethodS3("extractCNT", "SnpChipEffectFile", function(this, reference, units=NULL, chromosomes=NULL, addNames=TRUE, fields=c("Log2Ratio", "FreqB"), ..., order=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cdf <- getCdf(this);
  # Argument 'reference':
  if (inherits(this, "SnpChipEffectSet")) {
    className <- class(getOneFile(this))[1L];
  } else {
    className <- class(this)[1L];
  }
  reference <- Arguments$getInstanceOf(reference, className);

  # Argument 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  } else {
    units <- 1:nbrOfUnits(cdf);
  }

  # Argument 'chromosomes':
  if (!is.null(chromosomes)) {
    chromosomes <- Arguments$getIndices(chromosomes, max=999);
  }

  # Argument 'addNames':
  addNames <- Arguments$getLogical(addNames);
  
  # Argument 'fields':
  fields <- match.arg(fields, several.ok=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  gi <- getGenomeInformation(cdf);
  gp <- readDataFrame(gi, units=units, verbose=less(verbose, 25));

  # Subset by chromosome?
  if (!is.null(chromosomes)) {
    keep <- which(gp[,1] %in% chromosomes);
    units <- units[keep];
    gp <- gp[keep,,drop=FALSE];
    # Not needed anymore
    keep <- NULL;
  }

  if (order) {
    o <- order(gp[,1], gp[,2]);
    units <- units[o];
    gp <- gp[o,,drop=FALSE];
    # Not needed anymore
    o <- NULL;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  dataHead <- data.frame(
    Chromosome=gp[,1],
    Position=gp[,2]
  );
  # Not needed anymore
  gp <- NULL;

  if (addNames) {
    unitNames <- getUnitNames(cdf, units=units);
    dataHead <- cbind(ProbeSet=unitNames, dataHead);
    # Not needed anymore
    unitNames <- NULL;
  }
  # Not needed anymore
  cdf <- NULL;


  # Special case/Nothing (more) to do?
  if (length(fields) == 0) {
    return(dataHead);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract (theta, freqB) of the sample
  data <- extractTotalAndFreqB(this, units=units);

  # Rename
  colnames(data) <- gsub("^total", "Log2Ratio", colnames(data));
  colnames(data) <- gsub("^freqB", "FreqB", colnames(data));

  # Drop dimensions? [extract() to deal w/ either a matrix or an array]
  keep <- match(fields, colnames(data));
  data <- extract(data, "2"=keep, drop=FALSE);
  # Not needed anymore
  keep <- NULL;

  # Special case/ad hoc: If 'this' is a data set...
  if (length(dim(data)) == 3) {
    data <- wrap(data, map=list(1, 2:3));
  }

  # Calculate log ratios
  cc <- grep("^Log2Ratio", colnames(data));
  if (length(cc) > 0) {
    # Extract (total) theta:s of the reference
    thetaR <- extractTotalAndFreqB(reference, units=units)[,"total"];
    data[,cc] <- log2(data[,cc]/thetaR);
    # Not needed anymore
    thetaR <- NULL;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Build return structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- cbind(dataHead, data);
  # Not needed anymore
  dataHead <- NULL;

  data;  
}, protected=TRUE)

setMethodS3("extractCNT", "SnpChipEffectSet", extractCNT.SnpChipEffectFile, protected=TRUE)



############################################################################
# HISTORY:
# 2008-12-29
# o Added extractCNT() and writeCNT().
# o Created.
############################################################################
