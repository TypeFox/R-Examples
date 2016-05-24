#
# @keyword internal
setMethodS3("segmentByPruneCBS", "RawGenomicSignals", function(this, ...,       flavor=c("v1", "v2"), strict=FALSE, normalMean=0.0, debugPlot=FALSE, ylim=c(0,6), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'normalMean':
  normalMean <- Arguments$getDouble(normalMean);

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Use an old version?
  if (flavor == "v1") {
    return(segmentByPruneCBSv1(this, ..., strict=strict, normalMean=normalMean, debugPlot=debugPlot, ylim=ylim, verbose=verbose));
  }


  verbose && enter(verbose, "Segment using PruneCBS");

  if (debugPlot) {
    cols <- c("atomic island"="red",
              "ambigous atomic region"="orange",
              "extreme region"="purple");
  }

  res <- list();

  depth <- 1;
  cnr <- NULL;
  cn <- this;
  while(TRUE) {
    verbose && enter(verbose, sprintf("Segmentation depth %d", depth));
    cnrPrev <- cnr;

    fit <- segmentByCBS(cn);
    cnr <- extractCopyNumberRegions(fit);

    # Converged?
    if (equals(cnr, cnrPrev)) {
      cnr$type <- "constant";
      res[[depth]] <- cnr;
      verbose && cat(verbose, "Quitting, because nothing has changed since last round.");
      verbose && exit(verbose);
      break;
    }


    # Decrease the weights close to change points.
    # This will lower the risk for false change points
    # in succeeding segmentation iterations.
    # The weighting should be on the same "x scale" as the
    # segmentation operates on, that is, if it segments with
    # physical distances, then the weighting function should
    # be on that too, and if it segments with index distances,
    # then the weighting function should to that too.
    verbose && enter(verbose, "Updating locus-specific weights");
    w <- cn$w;
    if (is.null(w)) {
      w <- rep(1, times=nbrOfLoci(cn));
    }
    dx <- 100e3;
    regions <- as.data.frame(cnr)[,c("start","stop")];
    regions <- as.matrix(regions);
    regions2 <- regions;
    regions2[,"start"] <- regions2[,"start"] + dx;
    regions2[,"stop"] <- regions2[,"stop"] - dx;
    nok <- (regions2[,"start"] > regions[,"stop"]);
    regions2[nok,"start"] <- regions[nok,"stop"];
    nok <- (regions2[,"stop"] < regions[,"start"]);
    regions2[nok,"stop"] <- regions[nok,"start"];
    regions3 <- regions;
    regions3[,"stop"] <- regions2[,"start"];
    regions4 <- regions;
    regions4[,"start"] <- regions2[,"stop"];
    regions <- rbind(regions3, regions4);
    x <- getPositions(cn);
    for (kk in seq_len(nrow(regions))) {
      region <- regions[kk,,drop=TRUE];
      keep <- (region[1] <= x & x <= region[2]);
      w[keep] <- 0.2*w[keep];
    }
    cn$w <- w;
    # Not needed anymore
    w <- NULL;
    verbose && exit(verbose);

    if (debugPlot) {
#      if (depth > 1 && (depth-1) %% np == 0) {
#        readline("wait...");
#      }
      plot(cn, col="#aaaaaa", ylim=ylim, axes=FALSE);
      axis(side=2);
      drawLevels(cnr, col="black", lwd=3);
      verbose && print(verbose, cnr);
    }

    # Done?
    if (nbrOfRegions(cnr) == 1) {
      cnr$type <- "constant";
      res[[depth]] <- cnr;
      verbose && exit(verbose);
      break;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Find atomic aberrations of minimal length (H) to be pruned
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fitAR <- findAtomicAberrations(cnr, data=cn, H=1, ..., verbose=less(verbose,10));
    verbose && print(verbose, fitAR);

    # (a) Any atomic islands?
    minimalRegions <- fitAR$atomicIslands;
    n <- length(minimalRegions);
    if (n > 0) {
      type <- "atomic island";
    } else {
      # (b) Any ambigous atomic islands?
      minimalRegions <- fitAR$ambigousRegions;
      n <- length(minimalRegions);
      if (n > 0) {
        type <- "ambigous atomic region";
        cnrT <- subset(cnr, fitAR$ambigousRegions);
        verbose && cat(verbose, "atomicSiblings:");
        verbose && print(verbose, cnrT);

        # Find the furthest away from the normal state
        data <- as.data.frame(cnrT);
        normalMean <- 0.0;
        dMu <- (data$mean - normalMean);
        idx <- which.max(abs(dMu));
        minimalRegions <- fitAR$atomicRegions[idx];
      } else {
        type <- "extreme region";
        # (c) Extreme regions...
        # Find the furthest away from the normal state
        data <- as.data.frame(cnr);
        dMu <- (data$mean - normalMean);
        minimalRegions <- which.max(abs(dMu));
      }
    }

    n <- length(minimalRegions);
    verbose && printf(verbose, "Selected %d minimal regions that are %ss:\n", n, type);
    verbose && print(verbose, minimalRegions);

    if (n > 1 && strict) {
      verbose && cat(verbose, "Keeping only one minimal region (argument strict=TRUE)");
      # TO DO: For now, just the first one. /HB 2010-07-24.
      minimalRegions <- minimalRegions[1];
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Record the identified minimal regions
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    cnrX <- subset(cnr, minimalRegions);
    cnrX$type <- type;
    res[[depth]] <- cnrX;

    if (debugPlot) {
      drawLevels(cnrX, col=cols[type]);
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Drop minimal region(s)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    cnrT <- subset(cnr, -minimalRegions);
    cn <- extractRegions(cn, regions=cnrT);

    depth <- depth + 1;

    verbose && exit(verbose);
  } # while(...)

  types <- sapply(res, FUN=function(cnr) cnr$type);
  names(res) <- types;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reverse history
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- rev(res);

  verbose && exit(verbose);

  res;
}) # segmentByPruneCBS()





setMethodS3("segmentByPruneCBSv1", "RawGenomicSignals", function(this, ...,       flavor=c("v1", "v2"), strict=FALSE, normalMean=0.0, debugPlot=FALSE, ylim=c(0,6), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'normalMean':
  normalMean <- Arguments$getDouble(normalMean);

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Segment using PruneCBS");

  if (debugPlot) {
    cols <- c("atomic island"="red",
              "ambigous atomic region"="orange",
              "extreme region"="purple");
  }

  res <- list();

  depth <- 1;
  cnr <- NULL;
  cn <- this;
  while(TRUE) {
    verbose && enter(verbose, sprintf("Segmentation depth %d", depth));
    cnrPrev <- cnr;

    fit <- segmentByCBS(cn);
    cnr <- extractCopyNumberRegions(fit);

    # Converged?
    if (equals(cnr, cnrPrev)) {
      cnr$type <- "constant";
      res[[depth]] <- cnr;
      verbose && cat(verbose, "Quitting, because nothing has changed since last round.");
      verbose && exit(verbose);
      break;
    }


    # Decrease the weights close to change points.
    # This will lower the risk for false change points
    # in succeeding segmentation iterations.
    # The weighting should be on the same "x scale" as the
    # segmentation operates on, that is, if it segments with
    # physical distances, then the weighting function should
    # be on that too, and if it segments with index distances,
    # then the weighting function should to that too.
    verbose && enter(verbose, "Updating locus-specific weights");
    w <- cn$w;
    if (is.null(w)) {
      w <- rep(1, times=nbrOfLoci(cn));
    }
    dx <- 100e3;
    regions <- as.data.frame(cnr)[,c("start","stop")];
    regions <- as.matrix(regions);
    regions2 <- regions;
    regions2[,"start"] <- regions2[,"start"] + dx;
    regions2[,"stop"] <- regions2[,"stop"] - dx;
    nok <- (regions2[,"start"] > regions[,"stop"]);
    regions2[nok,"start"] <- regions[nok,"stop"];
    nok <- (regions2[,"stop"] < regions[,"start"]);
    regions2[nok,"stop"] <- regions[nok,"start"];
    regions3 <- regions;
    regions3[,"stop"] <- regions2[,"start"];
    regions4 <- regions;
    regions4[,"start"] <- regions2[,"stop"];
    regions <- rbind(regions3, regions4);
    x <- getPositions(cn);
    for (kk in seq_len(nrow(regions))) {
      region <- regions[kk,,drop=TRUE];
      keep <- (region[1] <= x & x <= region[2]);
      w[keep] <- 0.2*w[keep];
    }
    cn$w <- w;
    # Not needed anymore
    w <- NULL;
    verbose && exit(verbose);

    if (debugPlot) {
#      if (depth > 1 && (depth-1) %% np == 0) {
#        readline("wait...");
#      }
      plot(cn, col="#aaaaaa", ylim=ylim, axes=FALSE);
      axis(side=2);
      drawLevels(cnr, col="black", lwd=3);
      verbose && print(verbose, cnr);
    }

    # Done?
    if (nbrOfRegions(cnr) == 1) {
      cnr$type <- "constant";
      res[[depth]] <- cnr;
      verbose && exit(verbose);
      break;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Find minimal region(s) to be pruned
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fitAR <- findAtomicAberrations(cnr, data=cn, H=1, ..., verbose=less(verbose,10));
    verbose && print(verbose, fitAR);

    # (a) Any atomic islands?
    minimalRegions <- fitAR$atomicIslands;
    n <- length(minimalRegions);
    if (n > 0) {
      type <- "atomic island";
    } else {
      # (b) Any ambigous atomic islands?
      minimalRegions <- fitAR$ambigousRegions;
      n <- length(minimalRegions);
      if (n > 0) {
        type <- "ambigous atomic region";
        cnrT <- subset(cnr, fitAR$ambigousRegions);
        verbose && cat(verbose, "atomicSiblings:");
        verbose && print(verbose, cnrT);

        # Find the furthest away from the normal state
        data <- as.data.frame(cnrT);
        normalMean <- 0.0;
        dMu <- (data$mean - normalMean);
        idx <- which.max(abs(dMu));
        minimalRegions <- fitAR$atomicRegions[idx];
      } else {
        type <- "extreme region";
        # (c) Extreme regions...
        # Find the furthest away from the normal state
        data <- as.data.frame(cnr);
        dMu <- (data$mean - normalMean);
        minimalRegions <- which.max(abs(dMu));
      }
    }

    n <- length(minimalRegions);
    verbose && printf(verbose, "Selected %d minimal regions that are %ss:\n", n, type);
    verbose && print(verbose, minimalRegions);

    if (n > 1 && strict) {
      verbose && cat(verbose, "Keeping only one minimal region (argument strict=TRUE)");
      # TO DO: For now, just the first one. /HB 2010-07-24.
      minimalRegions <- minimalRegions[1];
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Record the identified minimal regions
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    cnrX <- subset(cnr, minimalRegions);
    cnrX$type <- type;
    res[[depth]] <- cnrX;

    if (debugPlot) {
      drawLevels(cnrX, col=cols[type]);
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Drop minimal region(s)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    cnrT <- subset(cnr, -minimalRegions);
    cn <- extractRegions(cn, regions=cnrT);

    depth <- depth + 1;

    verbose && exit(verbose);
  } # while(...)

  types <- sapply(res, FUN=function(cnr) cnr$type);
  names(res) <- types;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reverse history
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- rev(res);

  verbose && exit(verbose);

  res;
}, protected=TRUE) # segmentByPruneCBSv1()



############################################################################
# HISTORY:
# 2010-09-08
# o Updated to use findAtomicAberrations().
# 2010-07-24
# o CLEAN UP: Now the notation of the code better reflect the algorithm.
# o Now findAtomicRegions() returns ambigous atomic regions too.
# o Added argument 'ylim'.
# 2010-07-20
# o Added argument 'debugPlot'.
# 2010-07-19
# o Added trial version of segmentByPruneCBS().
# o TO DO: Down-weight loci that were close to earlier
#   change points in the succeeding segmentations.
# o Added prototype version of findAtomicRegions().
# o Added prototype version of callByPruning().
# o Created.
############################################################################
