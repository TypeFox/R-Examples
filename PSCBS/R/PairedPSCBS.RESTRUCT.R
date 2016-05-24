setMethodS3("extractSegments", "PairedPSCBS", function(this, idxs, ..., verbose=FALSE) {
  fit <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  updateSegRows <- function(segRows, idxs=NULL) {
    verbose && str(verbose, segRows);
    if (!is.null(idxs)) {
      segRows <- segRows[idxs,,drop=FALSE];
    }
#    verbose && cat(verbose, "Number of segments: ", nrow(segRows));
#    verbose && str(verbose, segRows);

    # Treat splitters separately
    isSplitter <- (is.na(segRows[,1]) & is.na(segRows[,2]));

    ns <- segRows[,2] - segRows[,1] + 1L;
#    verbose && cat(verbose, "Number of loci per segment:");
#    verbose && str(verbose, ns);

    ns <- ns[!isSplitter];
    from <- c(1L, cumsum(ns)[-length(ns)]+1L);
    to <- from + (ns - 1L);
    segRows[!isSplitter,1] <- from;
    segRows[!isSplitter,2] <- to;
    verbose && str(verbose, segRows);

    # Sanity check
    ns2 <- segRows[,2] - segRows[,1] + 1L;
    ns2 <- ns2[!isSplitter];
    stopifnot(all(ns2 == ns));

    segRows;
  } # updateSegRows()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'idxs':
  idxs <- Arguments$getIndices(idxs, max=nbrOfSegments(fit, splitters=TRUE));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extracting subset of segments");

  verbose && cat(verbose, "Number of segments: ", length(idxs));
  verbose && str(verbose, idxs);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit);
  tcnSegRows <- fit$tcnSegRows;
  dhSegRows <- fit$dhSegRows;
  segs <- getSegments(fit);
  params <- fit$params;

  # Sanity checks
  stopifnot(all(!is.na(data$chromosome) & !is.na(data$x)));
  stopifnot(length(tcnSegRows) == length(dhSegRows));

  # Sanity checks
  if (!params$joinSegments) {
    throw("Cannot extract subset of segments unless CNs are segmented using joinSegments=TRUE.");
  }

  if (params$flavor == "tcn,dh") {
    throw("NOT IMPLEMENTED: Extracting a subset of segments is not supported for flavor '", params$flavor, "'.");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Update table of segments");
  segsT <- segs[idxs,,drop=FALSE];
  verbose && str(verbose, segsT);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset data accordingly
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Update locus data");

  segRows <- tcnSegRows;
  segRows <- segRows[idxs,,drop=FALSE];
  from <- segRows[[1]];
  to <- segRows[[2]];
  ok <- (!is.na(from) & !is.na(to));
  from <- from[ok];
  to <- to[ok];
  keep <- logical(nrow(data));
  for (rr in seq(along=from)) {
    keep[from[rr]:to[rr]] <- TRUE;
  }
  keep <- which(keep);
  verbose && printf(verbose, "Identified %d (%.2f%%) of %d data rows:\n", length(keep), 100*length(keep)/nrow(data), nrow(data));
  verbose && str(verbose, keep);

  dataT <- data[keep,,drop=FALSE];
  verbose && str(verbose, dataT);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update 'segRows'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Update 'segRows'");

  segRows <- updateSegRows(tcnSegRows, idxs=idxs);
  d <- tcnSegRows[idxs,] - segRows;
  # Sanity check
  stopifnot(identical(d[,1], d[,2]));
  d <- d[,1];
  verbose && cat(verbose, "Row deltas:");
  verbose && str(verbose, d);

  tcnSegRows <- tcnSegRows[idxs,,drop=FALSE] - d;
  verbose && str(verbose, tcnSegRows);
  # Sanity checks
  segRows <- tcnSegRows;
  stopifnot(suppressWarnings(max(segRows, na.rm=TRUE)) <= nrow(dataT));
  drow <- segRows[-1,1] - segRows[-nrow(segRows),2];
  if (!all(is.na(drow) | (drow > 0))) {
    print(segRows);
    throw("INTERNAL ERROR: Generated 'tcnSegRows' is invalid, because it contains overlapping data chunks.");
  }

  dhSegRows <- dhSegRows[idxs,,drop=FALSE] - d;
  verbose && str(verbose, dhSegRows);
  # Sanity checks
  segRows <- dhSegRows;
  stopifnot(suppressWarnings(max(segRows, na.rm=TRUE)) <= nrow(dataT));
  drow <- segRows[-1,1] - segRows[-nrow(segRows),2];
  stopifnot(all(is.na(drow) | (drow > 0)));
  if (!all(is.na(drow) | (drow > 0))) {
    print(segRows);
    throw("INTERNAL ERROR: Generated 'dhSegRows' is invalid, because it contains overlapping data chunks.");
  }

  verbose && exit(verbose);


  # Create new object
  res <- fit;
  res$data <- dataT;
  res$output <- segsT;
  res$tcnSegRows <- tcnSegRows;
  res$dhSegRows <- dhSegRows;

  verbose && exit(verbose);

  res;
}, protected=TRUE) # extractSegments()



###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod mergeTwoSegments
#
# @title "Merge two neighboring segments"
#
# \description{
#   @get "title" by recalculating segment statistics.
# }
#
# @synopsis
#
# \arguments{
#  \item{left}{An @integer specifying the segments (left, left+1)
#    to be merged.}
#  \item{update}{If @TRUE, segment statistics are updated.}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" with one less segment.
# }
#
# @author "HB"
#
# \seealso{
#   To drop regions (a connected set of segments) see \code{dropRegions()}.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("mergeTwoSegments", "PairedPSCBS", function(this, left, update=TRUE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nbrOfSegments(this, splitters=TRUE);
  # Argument 'left':
  left <- Arguments$getIndex(left, max=nbrOfSegments-1L);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Merging two segments");
  verbose && printf(verbose, "Segments to be merged: %s & %s\n", left, left+1);
  verbose && cat(verbose, "Number of segments before merging: ", nbrOfSegments);
  verbose && cat(verbose, "Number of segments after merging: ", nbrOfSegments-1L);

  segs <- getSegments(this, splitters=TRUE);
  tcnSegRows <- this$tcnSegRows;
  dhSegRows <- this$dhSegRows;

  rows <- c(left,left+1);
  segsT <- segs[rows,,drop=FALSE];

  # Sanity check
  chrs <- segsT[["chromosome"]];
  if (chrs[1] != chrs[2]) {
    throw("Cannot merge segments that are on different chromosomes: ", chrs[1], " != ", chrs[2]);
  }

  # Merge segments
  segT <- segsT[1,];
  fields <- colnames(segsT);

  # (chromosome, tcnId, dhId)
  idxsUsed <- 1:3;

  # Starts
  idxs <- grep("Start$", fields);
  T <- as.matrix(segsT[,idxs,drop=FALSE]);
  segT[,idxs] <- colMins(T, na.rm=TRUE);
  idxsUsed <- c(idxsUsed, idxs);

  # Ends
  idxs <- grep("End$", fields);
  T <- as.matrix(segsT[,idxs,drop=FALSE]);
  segT[,idxs] <- colMaxs(T, na.rm=TRUE);
  idxsUsed <- c(idxsUsed, idxs);

  # Counts
  idxs <- grep("NbrOf", fields);
  segT[,idxs] <- colSums(segsT[,idxs,drop=FALSE]);
  idxsUsed <- c(idxsUsed, idxs);

  # "Invalidate" remaining entries
  if (update) {
    idxsTodo <- setdiff(seq(along=fields), idxsUsed);
    segT[,idxsTodo] <- NA;
  }

  # Update segment table
  segs[rows[1],] <- segT;
  segs <- segs[-rows[2],];

  # Update 'segRows' tables
  segRows <- tcnSegRows;
  segRows[rows[1],2] <- segRows[rows[2],2];
  segRows <- segRows[-rows[2],];
  tcnSegRows <- segRows;

  segRows <- dhSegRows;
  segRows[rows[1],2] <- segRows[rows[2],2];
  segRows <- segRows[-rows[2],];
  dhSegRows <- segRows;

  # Create results object
  res <- this;
  res$output <- segs;
  res$tcnSegRows <- tcnSegRows;
  res$dhSegRows <- dhSegRows;

  # Update the segment statistics?
  if (update) {
    res <- updateMeans(res);
  }

  verbose && exit(verbose);

  res;
}, private=TRUE)



############################################################################
# HISTORY:
# 2012-02-24
# o BUG FIX: The local updateSegRows() function inside extractSegments()
#   for PairedPSCBS would return incorrect and invalid row indices.
#   Copied ditto for CBS, which seems to work.
# o ROBUSTNESS: Added more sanity checks validating the correctness of
#   what is returned by extractSegments() for PairedPSCBS.
# 2012-01-09
# o ROBUSTNESS: Now extractSegments() for PairedPSCBS gives an informative
#   error message that it is not supported if CNs were segmented using
#   flavor "tcn,dh".
# 2011-10-16
# o Added argument 'update' to mergeTwoSegments().
# 2011-10-02
# o CLEANUP: Dropped empty callSegments() for PairedPSCBS.
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-04-08
# o BUG FIX: postsegmentTCN() for PairedPSCBS could generate an invalid
#   'tcnSegRows' matrix, where the indices for two consecutive segments
#   would overlap, which is invalid.
# 2011-04-05
# o BUG FIX: estimateHighDHQuantileAtAB() for PairedPSCBS would throw
#   an error on an undefined 'trim' if verbose output was used.
# 2011-02-17
# o Added arguments 'robust' and 'trim' to estimateMeanForDH().
# 2011-02-03
# o Added argument 'tauTCN' to estimateMeanForDH().
# 2011-01-27
# o Added flavor="DHskew" to estimateTauAB().
# o Added flavor="DH" to estimateTauAB() to estimate from DH instead
#   of hBAF.  As argued by the equations in the comments, these two
#   approaches gives virtually the same results.  The advantage with the
#   DH approach is that it requires one less degree of freedom.
# o Added estimateMeanForDH().
# 2011-01-18
# o BUG FIX: 'tcnSegRows' and 'dhSegRows' where not updated by
#   extractByRegions() for PairedPSCBS.
# 2011-01-14
# o Added estimateTauAB() for estimating the DeltaAB parameter.
# o Added estimateStdDevForHeterozygousBAF() for PairedPSCBS.
# o BUG FIX: extractByRegions() did not handle the case where multiple loci
#   at the same position are split up in two different segments.
# 2011-01-12
# o Added extractByRegions() and extractByRegion() for PairedPSCBS.
# o Now postsegmentTCN(..., force=TRUE) for PairedPSCBS also updates
#   the TCN estimates even for segments where the DH segmentation did
#   not find any additional change points.
# 2010-12-02
# o Now postsegmentTCN() assert that total number of TCN loci before
#   and after is the same.
# o Now postsegmentTCN() assert that joinSegment is TRUE.
# 2010-12-01
# o Now postsegmentTCN() checks if it is already postsegmented.
# 2010-11-30
# o TODO: postsegmentTCN() does not make sure of 'dhLociToExclude'. Why?
# o Now postsegmentTCN() recognizes the new 'tcnLociToExclude'.
# 2010-11-28
# o BUG FIX: postsegmentTCN() did not handle loci with the same positions
#   and that are split in two different segments.  It also did not exclude
#   loci with missing values.
# 2010-11-21
# o Adjusted postsegmentTCN() such that the updated TCN segment boundaries
#   are the maximum of the DH segment and the support by the loci.  This
#   means that postsegmentTCN() will work as expected both when signals
#   where segmented with 'joinSegments' being TRUE or FALSE.
# 2010-10-25
# o Now subsetByDhSegments() for PairedPSCBS handles the rare case when
#   markers with the same positions are split in two different segments.
# o Renamed subsetBySegments() for PairedPSCBS to subsetByDhSegments().
# 2010-09-26
# o Now subsetBySegments() for PairedPSCBS handles multiple chromosomes.
# o Now postsegmentTCN() PairedPSCBS handles multiple chromosomes.
# 2010-09-21
# o Added postsegmentTCN() for PairedPSCBS.
# 2010-09-19
# o BUG FIX: plot() used non-defined nbrOfLoci; now length(x).
# 2010-09-15
# o Added subsetBySegments().
# o Added linesC1C2() and arrowsC1C2().
# o Now the default 'cex' for pointsC1C2() corresponds to 'dh.num.mark'.
# o Now extractTotalAndDH() also returns 'dh.num.mark'.
# 2010-09-08
# o Added argument 'add=FALSE' to plot().
# o Added plotC1C2().
# o Added extractTotalAndDH() and extractMinorMajorCNs().
# 2010-09-04
# o Added drawLevels() for PairedPSCBS.
# o Added as.data.frame() and print() for PairedPSCBS.
# 2010-09-03
# o Added plot() for PairedPSCBS.
# o Created.
############################################################################
