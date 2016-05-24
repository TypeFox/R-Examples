setMethodS3("getAM", "AromaUnitTotalCnBinaryFile", function(this, other, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ugp <- getAromaUgpFile(this);
  # Argument 'other':
  if (is.null(other)) {
    # Do not calculate ratios relative to a reference
    other <- "none";
  }
  if (is.character(other)) {
    choices <- c("none", "constant(1)", "constant(2)");
    if (!is.element(other, choices)) {
      throw(sprintf("Argument 'other' should be one of %s: %s", paste(dQuote(choices), collapse=", "), other[1]));
    }
    other <- match.arg(other, choices=choices);
  } else {
    other <- Arguments$getInstanceOf(other, "CopyNumberDataFile");
  }
  
  # Argument 'units':
  if (is.null(units)) {
    units <- seq_len(nbrOfUnits(ugp));
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(ugp));
  }
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  


  verbose && enter(verbose, "Getting (A,M)-transformed chip effects");

  nunits <- length(units);

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get thetas from the sample
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving sample thetas");

  theta <- this[units,1, drop=TRUE];
  nTheta <- length(theta);
  if (!identical(nTheta, nunits)) {
    verbose && str(verbose, theta);
    verbose && print(verbose, nunits);
    throw("Internal error: The number of chip-effect values is not equal to the number of units requested: ", nTheta, " != ", nunits);
  }

  # Unlog?
  logBase0 <- NULL;
  if (hasTag(this, "log2ratio")) {
    logBase0 <- 2;
  } else if (hasTag(this, "log10ratio")) {
    logBase0 <- 10;
  } else if (hasTag(this, "logRatio")) {
    logBase0 <- 10;
  }
  if (!is.null(logBase0)) {
    theta <- logBase0^theta;
    verbose && printf(verbose, "Transformed theta = %f^theta\n", logBase0);
  }

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get thetas from the other
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.character(other)) {
    verbose && enter(verbose, "Calculating special other thetas");
    verbose && cat(verbose, "Argument 'other': ", other);

    if (other == "none") {
      thetaR <- 1;
    } else if (other == "constant(1)") {
      thetaR <- 1;
    } else if (other == "constant(2)") {
      thetaR <- 2;
    } else {
      throw("Unknown value of 'other': ", other);
    }

    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "Retrieving other thetas");
  
    # Get the other theta estimates
    thetaR <- other[units,1, drop=TRUE];
    stopifnot(identical(length(thetaR), nTheta));

    # Unlog?
    logBase0 <- NULL;
    if (hasTag(other, "log2ratio")) {
      logBase0 <- 2;
    } else if (hasTag(other, "log10ratio")) {
      logBase0 <- 10;
    } else if (hasTag(other, "logRatio")) {
      logBase0 <- 10;
    }

    if (!is.null(logBase0)) {
      thetaR <- logBase0^thetaR;
      verbose && printf(verbose, "Transformed thetaR = %f^thetaR\n", logBase0);
    }

    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate raw copy numbers relative to the other
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  M <- theta/thetaR;
  A <- theta*thetaR;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # On the log2 scale?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  M <- log(M, base=2);
  A <- log(A, base=2)/2;
  stopifnot(identical(length(M), nTheta));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the return matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  am <- matrix(c(A,M), ncol=2)
  colnames(am) <- c("A", "M");
  rownames(am) <- as.character(units);

  verbose && exit(verbose);

  am;
}) # getAM()



setMethodS3("getXAM", "AromaUnitTotalCnBinaryFile", function(this, other, chromosome, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'other':
  # Validated by getAM() below.

  # Argument 'chromosome':
  chromosome <- Arguments$getCharacter(chromosome);  # integer? /HB 2009-05-18

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting (X,A,M)-transformed chip effects");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve genome information, i.e. chromosome positions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving genome information");
  ugp <- getAromaUgpFile(this);

  # Extract the units from the given chromosome.  Requested units not on
  # chromosome are ignored.
  if (!is.null(units)) {
    verbose && str(verbose, units);
  }

  units0 <- getUnitsOnChromosome(ugp, chromosome=chromosome, verbose=less(verbose));
  if (!is.null(units)) {
    units0 <- intersect(units0, units);
  }
  units <- units0;

  nunits <- length(units);
  if (nunits == 0) {
    throw("No units found on requested chromosome: ", chromosome);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the relative copy-number estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  am <- getAM(this, other=other, units=units, verbose=less(verbose));

  # Get the unit indices for all unit groups
  units <- rownames(am);
  # Sanity check
  if (is.null(units)) {
    throw("Internal error: getAM() did not return unit names: NULL");
  }
  units <- as.integer(units);
  verbose && cat(verbose, "Units read:");
  verbose && str(verbose, units);
  
  # Get the positions of all unit groups
  x <- getPositions(ugp, units=units, verbose=less(verbose));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove units for which we have no position information
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  keep <- which(!is.na(x));
  nexcl <- length(x) - length(keep);
  if (nexcl > 0) {
    msg <- sprintf("Could not find position information on %d unit groups: ", nexcl);
    verbose && cat(verbose, msg);
    warning(msg);
    x <- x[keep];
    units <- units[keep];
  }

  am <- am[,c("M","A"), drop=FALSE]; # Ad hoc /HB 2007-02-19
  xam <- cbind(x=x, am);

  verbose && cat(verbose, "(X,A,M):");
  verbose && str(verbose, xam);

  verbose && exit(verbose);

  xam;
}, protected=TRUE) # getXAM()


############################################################################
# HISTORY:
# 2012-10-21
# o Now getAM() accepts values "none", "constant(1)", and "constant(2)"
#   for argument 'other'.
# o CLEANUP: Now argument 'other' of getXAM() is validated by getAM().
# 2009-11-19
# o Created.
############################################################################
