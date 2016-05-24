setMethodS3("applyByRegion", "PairedPSCBS", function(fit, FUN, ..., subset=NULL, append=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'FUN':
  stopifnot(is.function(FUN));

  # Argument 'subset':
  if (!is.null(subset)) {
    subset <- Arguments$getIndices(subset, range=c(1, nbrOfSegments(fit)));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Apply function region by region");
  verbose && cat(verbose, "Segments:");
  verbose && str(verbose, subset);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit);
  segs <- getSegments(fit);
  tcnSegRows <- fit$tcnSegRows;
  dhSegRows <- fit$dhSegRows;
  params <- fit$params;

  # Sanity checks
  if (!params$joinSegments) {
    throw("Cannot applyByRegion() unless PSCNs are segmented using joinSegments=TRUE.");
  }
  dataRows <- tcnSegRows;

  # Sanity checks
  stopifnot(all(!is.na(data$chromosome) & !is.na(data$x)));
  stopifnot(length(tcnSegRows) == length(dhSegRows));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For each segment...
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # Allocate result objects?
  if (append) {
    dataN <- outputN <- dataRowsN <- NULL;
  }

  if (is.null(subset)) {
    subset <- seq_len(nbrOfSegments);
  }

  for (rr in subset) {
    verbose && enter(verbose, sprintf("Segment #%d of %d", rr, nbrOfSegments));

    # Extract segment
    segRR <- segs[rr,,drop=FALSE];

    # Nothing todo?
    if (is.na(segRR[["tcnId"]]) && is.na(segRR[["dhId"]])) {
      verbose && cat(verbose, "A divider. Nothing to do.");
      outputN <- rbind(outputN, NA);
      dataRowsN <- rbind(dataRowsN, NA);
      verbose && exit(verbose);
      next;
    }

    verbose && str(verbose, segRR, level=-20);

    # Extract data
    dataRowsRR <- dataRows[rr,,drop=FALSE];
    from <- dataRowsRR[[1]];
    to <- dataRowsRR[[2]];
    ok <- (!is.na(from) & !is.na(to));
    from <- from[ok];
    to <- to[ok];
    keep <- logical(nrow(data));
    for (kk in seq(along=from)) {
      keep[from[kk]:to[kk]] <- TRUE;
    }
    dataRowsRR <- which(keep);
    verbose && printf(verbose, "Identified %d (%.2f%%) of %d data rows:\n", length(dataRowsRR), 100*length(dataRowsRR)/nrow(data), nrow(data));
    verbose && str(verbose, dataRowsRR);
    dataRR <- data[dataRowsRR,,drop=FALSE];
    verbose && str(verbose, dataRR, level=-20);

    verbose && enter(verbose, "Applying function 'FUN' to segment");
    resRR <- FUN(rr, segRR, dataRR, ...);
    verbose && cat(verbose, "Returned result:");
    verbose && str(verbose, resRR, level=-20);
    verbose && exit(verbose);

    # Nothing to update/store?
    if (!is.list(resRR)) {
      verbose && cat(verbose, "Nothing more to do for this segment since nothing was returned: ", rr);
      verbose && exit(verbose);
      next;
    }

    # Extract return data
    dataRRN <- resRR$data;
    segRRN <- resRR$output;
    # Sanity check
    stopifnot(!is.null(dataRRN));
    stopifnot(is.data.frame(dataRRN));
    stopifnot(!is.null(segRRN));
    stopifnot(is.data.frame(segRRN));

    if (append) {
      # Modified locus-level data
      dataRowsRRN <- c(1L, nrow(dataRRN));
      if (!is.null(dataN)) {
        dataRowsRRN <- dataRowsRRN + nrow(dataN);
      }
      dataN <- rbind(dataN, dataRRN);
      # Sanity checks
      stopifnot(nrow(dataN) == max(dataRowsN, na.rm=TRUE));

      # Update segment table?
      outputN <- rbind(outputN, segRRN);
      dataRowsN <- rbind(dataRowsN, dataRowsRRN);
      # Sanity check
      stopifnot(nrow(outputN) == nrow(dataRowsN));
      # Sanity checks
      stopifnot(nrow(dataN) == max(dataRowsN, na.rm=TRUE));
    } else {
      # Modified locus-level data
      verbose && enter(verbose, "Updating locus-level data");
      # Sanity check
      stopifnot(dim(dataRRN) == dim(dataRR));
      stopifnot(length(dataRowsRR) == nrow(dataRRN));
      data[dataRowsRR,] <- dataRRN;
      str(data[dataRowsRR,]);
      verbose && exit(verbose);

      # Modified segment data
      verbose && enter(verbose, "Updating segment data");
      # Sanity check
      stopifnot(dim(segRRN) == dim(segRR));
      segs[rr,] <- segRRN;
      verbose && exit(verbose);
    }

    # Not needed anymore
    dataRRN <- segRRN <- NULL;
    dataRR <- segRR <- NULL;
    resRR <- NULL;

    verbose && exit(verbose);
  } # for (rr ...)

  if (append) {
    if (!is.null(dataRowsN)) {
      rownames(dataRowsN) <- NULL;
      colnames(dataRowsN) <- colnames(dataRows);
      dataRowsN <- as.data.frame(dataRowsN);
      # Sanity checks
      stopifnot(!is.null(dataN));
      stopifnot(!is.null(outputN));
      stopifnot(!is.null(dataRowsN));

      data <- dataN;
      segs <- outputN;

      # Not needed anymore
      dataN <- outputN <- NULL;
    }
  }

  # Return result
  res <- fit; # "clone"
  res$data <- data;
  res$output <- segs;

  # Not needed anymore
  data <- segs <- NULL;

  # Update segment-to-locus index tables
  if (append && !is.null(dataRowsN)) {
    res$tcnSegRows <- dataRowsN;
    res$dhSegRows <- dataRowsN;  # Is this really the case? /HB 2011-01-17
  }

  verbose && exit(verbose);

  res;
}, private=TRUE)


.addC1C2WithStatitics <- function(rr, output, data, robust=TRUE, ...) {
  # Calculate locus-level (C1,C2)
  C <- data$CT;
  rho <- data$rho;
  C1 <- 1/2 * (1 - rho) * C;
  C2 <- C - C1;
  CC <- data.frame(C1=C1, C2=C2);

  if (robust) {
    meanFcn <- function(x, ...) median(x, na.rm=TRUE);
    sdFcn <- function(x, ...) mad(x, na.rm=TRUE);
  } else {
    meanFcn <- function(x, ...) mean(x, na.rm=TRUE);
    sdFcn <- function(x, ...) sd(x, na.rm=TRUE);
  }

  # Calculate region-level (C1,C2) means and std devs.
  muCC <- apply(CC, MARGIN=2, FUN=meanFcn);
  sigmaCC <- apply(CC, MARGIN=2, FUN=sdFcn);
  rhoCC <- cor(CC[,1], CC[,2], use="pairwise.complete.obs");

  names(muCC) <- c("c1Avg", "c2Avg");
  names(sigmaCC) <- c("c1Sd", "c2Sd");

  # Update data
  data <- cbind(data, CC);


  # Update segment table
  outputT <- c(muCC, sigmaCC, c1c2.cor=rhoCC);
  outputT <- as.list(outputT);
  outputT <- as.data.frame(outputT);

  output <- cbind(output, outputT);

  list(data=data, output=output);
} # .addC1C2WithStatitics()


.addCACBWithStatitics <- function(rr, output, data, beta=c("betaTN", "betaT"), stratifyBy=c("all", "hets", "homs"), robust=TRUE, ...) {
  # Argument 'beta':
  beta <- match.arg(beta);

  # Argument 'stratifyBy':
  stratifyBy <- match.arg(stratifyBy);

  # Calculate locus-level (CA,CB)
  C <- data$CT;
  beta <- data[[beta]];
  CB <- beta * C;
  CA <- C - CB;
  CC <- data.frame(CA=CA, CB=CB);

  # Update data
  data <- cbind(data, CC);


  if (robust) {
    meanFcn <- function(x, ...) median(x, na.rm=TRUE);
    sdFcn <- function(x, ...) mad(x, na.rm=TRUE);
  } else {
    meanFcn <- function(x, ...) mean(x, na.rm=TRUE);
    sdFcn <- function(x, ...) sd(x, na.rm=TRUE);
  }

  if (stratifyBy == "hets") {
    muN <- data$muN;
    keep <- (muN == 1/2);
    CC <- CC[keep,,drop=FALSE];
  } else if (stratifyBy == "homs") {
    muN <- data$muN;
    keep <- (muN == 0 | muN == 1);
    CC <- CC[keep,,drop=FALSE];
  }

  # Calculate region-level (CA,CB) means and std devs.
  muCC <- apply(CC, MARGIN=2, FUN=meanFcn);
  sigmaCC <- apply(CC, MARGIN=2, FUN=sdFcn);
  if (nrow(CC) < 3) {
    rhoCC <- NA_real_;
  } else {
    rhoCC <- cor(CC[,1], CC[,2], use="pairwise.complete.obs");
  }
  names(muCC) <- c("caAvg", "cbAvg");
  names(sigmaCC) <- c("caSd", "cbSd");

  # Update segment table
  outputT <- c(muCC, sigmaCC, cacbCor=rhoCC);
  outputT <- as.list(outputT);
  outputT <- as.data.frame(outputT);

  output <- cbind(output, outputT);

  list(data=data, output=output);
} # .addCACBWithStatitics()



#############################################################################
# HISTORY:
# 2013-10-21
# o Added argument 'subset' to applyByRegion() for PairedPSCBS.
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-01-27
# o Added .addCACBWithStatitics().
# o Added .addC1C2WithStatitics().
# o Added applyByRegion().
# o Created.
#############################################################################
