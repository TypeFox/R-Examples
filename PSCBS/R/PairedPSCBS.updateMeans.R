setMethodS3("updateMeans", "PairedPSCBS", function(fit, from=c("loci", "segments"), adjustFor=NULL, ..., avgTCN=c("asis", "mean", "median"), avgDH=c("asis", "mean", "median"), clear=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'from':
  from <- match.arg(from);

  # Argument 'adjustFor':
  if (!is.null(adjustFor)) {
    adjustFor <- Arguments$getCharacters(adjustFor);
    adjustFor <- tolower(adjustFor);
    knownValues <- c("ab", "loh", "roh");
    adjustFor <- match.arg(adjustFor, choices=knownValues, several.ok=TRUE);
  }

  # Argument 'avgTCN' & 'avgDH':
  avgTCN <- match.arg(avgTCN);
  avgDH <- match.arg(avgDH);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Updating mean level estimates");
  verbose && cat(verbose, "Adjusting for:");
  verbose && print(verbose, adjustFor);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up averaging functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (avgTCN == "asis" || avgDH == "asis") {
    est <- fit$params$meanEstimators;
    if (avgTCN == "asis") {
      avgTCN <- est$tcn;
      if (is.null(avgTCN)) avgTCN <- "mean";
      avgTCN <- match.arg(avgTCN);
    }
    if (avgDH == "asis") {
      avgDH <- est$dh;
      if (is.null(avgDH)) avgDH <- "mean";
      avgDH <- match.arg(avgDH);
    }
  }

  avgList <- list(
    tcn = get(avgTCN, mode="function"),
    dh = get(avgDH, mode="function")
  );


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  segs <- getSegments(fit, splitters=TRUE);
  segRows <- list(tcn=fit$tcnSegRows, dh=fit$dhSegRows);
  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert that adjustments can be made
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.element("ab", adjustFor)) {
    if (!is.element("abCall", names(segs))) {
      adjustFor <- setdiff(adjustFor, "ab");
      throw("Cannot adjust for AB, because they haven't been called.");
    }
  }

  if (is.element("loh", adjustFor)) {
    if (!is.element("lohCall", names(segs))) {
      adjustFor <- setdiff(adjustFor, "loh");
      throw("Cannot adjust for LOH, because they haven't been called.");
    }
  }

  if (is.element("roh", adjustFor)) {
    if (!is.element("rohCall", names(segs))) {
      adjustFor <- setdiff(adjustFor, "roh");
      throw("Cannot adjust for ROH, because they haven't been called.");
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update the (TCN,DH) mean levels from locus-level data?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (from == "loci") {
    data <- getLocusData(fit);
    chromosome <- data$chromosome;
    x <- data$x;
    CT <- data$CT;
    rho <- data$rho;

    isSplitter <- isSegmentSplitter(fit);
    for (ss in seq(length=nbrOfSegments)[!isSplitter]) {
      verbose && enter(verbose, sprintf("Segment %d of %d", ss, nbrOfSegments));
      seg <- segs[ss,];
      verbose && print(verbose, seg);

      chr <- seg[["chromosome"]];
      chrTag <- sprintf("chr%02d", chr);

      for (what in c("tcn", "dh")) {
        segRow <- segRows[[what]][ss,];

        # (a) A splitter - nothing todo?
        if (!is.finite(segRow[[1]]) || !is.finite(segRow[[2]])) {
          next;
        }

        # (b) Identify units (loci)
        units <- segRow[[1]]:segRow[[2]];

        # (c) Adjust for missing values
        if (what == "tcn") {
          value <- CT;
        } else if (what == "dh") {
          value <- rho;
        }
        keep <- which(!is.na(value[units]));
        units <- units[keep];

        # (d) Update mean
        avgFUN <- avgList[[what]];
        gamma <- avgFUN(value[units]);

        # Sanity check
        stopifnot(length(units) == 0 || !is.na(gamma));

        # Update the segment boundaries, estimates and counts
        key <- paste(what, "Mean", sep="");
        seg[[key]] <- gamma;
      }

      verbose && print(verbose, seg);

      segs[ss,] <- seg;

      verbose && exit(verbose);
    } # for (ss ...)
  } # if (from ...)



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Adjust segment means from various types of calls
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (length(adjustFor) > 0) {
    verbose && enter(verbose, "Adjusting segment means");
    verbose && cat(verbose, "Adjusting for:");
    verbose && print(verbose, adjustFor);

    if (is.element("ab", adjustFor)) {
      verbose && enter(verbose, "Adjusting for AB");
      calls <- segs$abCall;
      segs$dhMean[calls] <- 0;
      verbose && exit(verbose);
    }

    if (is.element("loh", adjustFor)) {
      verbose && enter(verbose, "Adjusting for LOH");
      calls <- segs$lohCall;
      segs$dhMean[calls] <- 1;
      verbose && exit(verbose);
    }

    if (is.element("roh", adjustFor)) {
      verbose && enter(verbose, "Adjusting for ROH");
      calls <- segs$rohCall;
      segs$dhMean[calls] <- NA_real_;
      verbose && exit(verbose);
    }

    verbose && exit(verbose);
  } # if (length(adjustFor) > 0)


  # Update
  fit$output <- segs;
  fit <- setMeanEstimators(fit, tcn=avgTCN, dh=avgDH);
  if (clear) {
    fit <- clearBootstrapSummaries(fit);
  }

  # Update (C1,C2) mean levels
  fit <- updateMeansC1C2(fit, verbose=verbose);

  verbose && exit(verbose);

  fit;
}, private=TRUE) # updateMeans()


setMethodS3("updateMeansC1C2", "PairedPSCBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Updating (C1,C2) segment mean levels");
  segs <- getSegments(fit);

  if (nrow(segs) > 0L) {
    tcn <- segs$tcnMean;
    dh <- segs$dhMean;

    C1 <- 1/2*(1-dh)*tcn;
    C2 <- tcn - C1;

    segs$c1Mean <- C1;
    segs$c2Mean <- C2;

    # Preserve (C1,C2) swaps / change-point flips?
    swap <- segs$c1c2Swap;
    if (!is.null(swap)) {
      swap <- which(swap);
      if (length(swap) > 0L) {
        segs[swap, c("c1Mean","c2Mean")] <- segs[swap, c("c2Mean","c1Mean")];
      }
    }

    fit$output <- segs;
  }

  verbose && exit(verbose);

  fit;
}, protected=TRUE) # updateMeansC1C2()



##############################################################################
# HISTORY
# 2014-03-26
# o BUG FIX: updateMeansC1C2() for PairedPSCBS did not handle missing
#   values (=splitters) in the 'c1c2Swap' field.
# 2014-03-25
# o BUG FIX: updateMeans() for PairedPSCBS and NonPairedPSCBS returned the
#   incorrect DH segment levels for region in AB if adjustFor="ab" and
#   likewise for segments in LOH if adjustFor="loh".
# 2013-11-23
# o Now updateMeans(..., clear=TRUE) clears bootstrap summaries, otherwise
#   not.
# 2013-10-26
# o Now updateMeans() for PairedPSCBS always clears the bootstrap summaries.
# o Added updateMeansC1C2().
# o updateMeans() for PairedPSCBS did not preserve the (C1,C2) swaps.
# 2012-04-21
# o CLEANUP: Removed unused objects in updateMeans().
# 2011-11-12
# o Added arguments 'from' and 'adjustFor' to updateMeans().
# 2011-01-16
# o BUG FIX: updateMeans() save to the incorrect column names.
# 2011-01-12
# o Added updateMeans() for PairedPSCBS.
##############################################################################
