setMethodS3("updateMeansTogether", "PairedPSCBS", function(fit, idxList, ..., avgTCN=c("mean", "median"), avgDH=c("mean", "median"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nbrOfSegments(fit, splitters=TRUE);

  # Argument 'idxList':
  if (!is.list(idxList)) {
    idxList <- list(idxList);
  }
  idxList <- lapply(idxList, FUN=function(idxs) {
    idxs <- Arguments$getIndices(idxs, max=nbrOfSegments);
    sort(unique(idxs));
  });

  # Argument 'avgTCN' & 'avgDH':
  avgTCN <- match.arg(avgTCN);
  avgDH <- match.arg(avgDH);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Updating mean level estimates of multiple segments");

  verbose && cat(verbose, "Segments:");
  verbose && str(verbose, idxList);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up averaging functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  avgList <- list(
    tcn = get(avgTCN, mode="function"),
    dh = get(avgDH, mode="function")
  );


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit);

  segs <- getSegments(fit, splitters=TRUE);

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Total number of segments: ", nbrOfSegments);

  for (ss in seq(along=idxList)) {
    idxs <- idxList[[ss]];

    fitT <- extractSegments(fit, idxs);
    verbose && cat(verbose, "Number of segments: ", nbrOfSegments(fitT));

    dataT <- getLocusData(fitT);
    segsT <- getSegments(fitT);

    CT <- dataT$CT;
    rho <- dataT$rho;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Update the TCN segments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Recalculate (TCN,DH,C1,C2) means");
    naValue <- NA_real_;
    mus <- c(tcn=naValue, dh=naValue, c1=naValue, c2=naValue);
    for (key in c("tcn", "dh")) {
      avgFUN <- avgList[[key]];

      # (c) Adjust for missing values
      if (key == "tcn") {
        value <- CT;
      } else if (key == "dh") {
        value <- rho;
      }
      keep <- which(!is.na(value));

      # (d) Update mean
      gamma <- avgFUN(value[keep]);

      # Sanity check
      stopifnot(length(gamma) == 0 || !is.na(gamma));

      mus[key] <- gamma;
    } # for (what ...)

    mus["c1"] <- 1/2*(1-mus["dh"])*mus["tcn"];
    mus["c2"] <- mus["tcn"] - mus["c1"];
    names(mus) <- sprintf("%sMean", names(mus));
    verbose && print(verbose, mus);
    verbose && exit(verbose);

    for (key in names(mus)) {
      segs[idxs,key] <- mus[key];
    }
  } # for (ss ...)

  # Return results
  res <- fit;
  res$output <- segs;
  res <- setMeanEstimators(res, tcn=avgTCN, dh=avgDH);

  verbose && exit(verbose);

  res;
}, private=TRUE) # updateMeansTogether()



############################################################################
# HISTORY:
# 2011-11-28
# o Dropped kmeansCNs() stub.
# o Added Rdoc comments.
# o Now hclustCNs() also handles segments with missing (C1,C2) levels,
#   which for instance can happen after calling ROH.
# 2011-10-14
# o Implemented hclustCNs() and pruneByHClust() for AbstractCBS.
# o Implemented extractCNs() for PairedPSCBS.
# o Created.
############################################################################
