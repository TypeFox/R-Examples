setMethodS3("shiftTCN", "CBS", function(fit, shift, update=TRUE, ...) {
  # Argument 'shift':
  shift <- Arguments$getDouble(shift, disallow=c("NA", "NaN", "Inf"));

  data <- getLocusData(fit);
  data$y <- data$y + shift;
  fit$data <- data;
  # Not needed anymore
  data <- NULL;

  if (update) {
    fit <- updateMeans(fit, ...);
  }

  fit;
}, protected=TRUE)


###########################################################################/**
# @set "class=CBS"
# @RdocMethod append
#
# @title "Appends one segmentation result to another"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{x, other}{The two @see "CBS" objects to be combined.}
#  \item{other}{A @see "PSCBS" object.}
#  \item{addSplit}{If @TRUE, a "divider" is added between chromosomes.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "CBS" object of the same class as argument \code{x}.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("append", "CBS", function(x, other, addSplit=TRUE, ...) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'other':
  other <- Arguments$getInstanceOf(other, class(this)[1]);
  for (field in c("data", "output")) {
    dataA <- this[[field]]
    dataB <- other[[field]]
    namesA <- colnames(dataA)
    namesB <- colnames(dataB)
    if (!all(namesA == namesB)) {
      throw(sprintf("Cannot merge %s objects. Arguments 'other' and 'this' has different sets of columns in field '%s': {%s} [n=%d] != {%s} [n=%d]", class(this)[1], field, paste(namesA, collapse=", "), length(namesA), paste(namesB, collapse=", "), length(namesB)))
    }
  }

  # Argument 'addSplit':
  addSplit <- Arguments$getLogical(addSplit);


  # Allocate results
  res <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locus data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(this);
  res$data <- rbind(data, getLocusData(other));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segmentation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  indexOffset <- nrow(data);
  fields <- c("output", "segRows");
  for (field in fields[-1]) {
    other[[field]] <- other[[field]] + indexOffset;
  }

  splitter <- if (addSplit) NA else NULL;
  for (field in fields) {
    res[[field]] <- rbind(this[[field]], splitter, other[[field]]);
    rownames(res[[field]]) <- NULL;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ksT <- this$params$knownSegments;
  ksT$length <- NULL;  # In case it's been added
  ksO <- other$params$knownSegments;
  ksO$length <- NULL;  # In case it's been added
  res$params$knownSegments <- rbind(ksT, ksO);


  # Sanity check
  ns <- sapply(res[fields], FUN=nrow);
  stopifnot(all(ns == ns[1]));

  res;
}) # append()



setMethodS3("extractSegments", "CBS", function(this, idxs, ..., verbose=FALSE) {
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
  segRows <- fit$segRows;
  segs <- getSegments(fit);
  params <- fit$params;

  # Sanity checks
  stopifnot(all(!is.na(data$chromosome) & !is.na(data$x)));

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

  segRowsT <- segRows[idxs,,drop=FALSE];
  from <- segRowsT[[1]];
  to <- segRowsT[[2]];
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
  segRowsT <- updateSegRows(segRowsT);
  d <- segRows[idxs,] - segRowsT;

  # Sanity check
  stopifnot(identical(d[,1], d[,2]));
  d <- d[,1];
  verbose && cat(verbose, "Row deltas:");
  verbose && str(verbose, d);

  segRows <- segRows[idxs,,drop=FALSE] - d;
  verbose && str(verbose, segRows);
  # Sanity checks
  stopifnot(suppressWarnings(max(segRows, na.rm=TRUE)) <= nrow(dataT));
  drow <- segRows[-1,1] - segRows[-nrow(segRows),2];
  stopifnot(all(is.na(drow) | (drow > 0)));
  if (!all(is.na(drow) | (drow > 0))) {
    print(segRows);
    throw("INTERNAL ERROR: Generated 'segRows' is invalid, because it contains overlapping data chunks.");
  }

  verbose && exit(verbose);


  # Create new object
  res <- fit;
  res$data <- dataT;
  res$output <- segsT;
  res$segRows <- segRows;

  verbose && exit(verbose);

  res;
}, protected=TRUE) # extractSegments()



setMethodS3("mergeTwoSegments", "CBS", function(this, left, update=TRUE, verbose=FALSE, ...) {
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

  segs <- getSegments(this);
  segRows <- this$segRows;

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
  idxsUsed <- c();

  # (id) [as in label]
  idxs <- grep("(I|i)d$", fields);
  idxsUsed <- c(idxsUsed, idxs);

  # (chromosome)
  idxs <- grep("chromosome$", fields);
  idxsUsed <- c(idxsUsed, idxs);

  # Starts
  idxs <- grep("(S|s)tart$", fields);
  T <- as.matrix(segsT[,idxs,drop=FALSE]);
  segT[,idxs] <- colMins(T, na.rm=TRUE);
  idxsUsed <- c(idxsUsed, idxs);

  # Ends
  idxs <- grep("(E|e)nd$", fields);
  T <- as.matrix(segsT[,idxs,drop=FALSE]);
  segT[,idxs] <- colMaxs(T, na.rm=TRUE);
  idxsUsed <- c(idxsUsed, idxs);

  # Counts
  idxs <- grep("(N|n)brOf", fields);
  segT[,idxs] <- colSums(segsT[,idxs,drop=FALSE]);
  idxsUsed <- c(idxsUsed, idxs);

  # "Invalidate" remaining entries
  idxsTodo <- setdiff(seq(along=fields), idxsUsed);
  segT[,idxsTodo] <- NA;

  # Update segment table
  segs[rows[1],] <- segT;
  segs <- segs[-rows[2],];

  # Update 'segRows' tables
  segRows[rows[1],2] <- segRows[rows[2],2];
  segRows <- segRows[-rows[2],];

  # Create results object
  res <- this;
  res$output <- segs;
  res$segRows <- segRows;

  # Update the segment statistics?
  if (update) {
    res <- updateMeans(res);
  }

  verbose && exit(verbose);

  res;
}, protected=TRUE) # mergeTwoSegments()



############################################################################
# HISTORY:
# 2012-09-13
# o Added shiftTCN() for CBS.
# 2012-02-24
# o ROBUSTNESS: Added more sanity checks validating the correctness of
#   what is returned by extractSegments() for CBS.
# 2011-11-17
# o BUG FIX: extractSegments() for CBS would throw an error when
#   there were multiple chromosomes.
# 2011-11-15
# o BUG FIX: extractSegments() for CBS would throw an error, because in
#   most cases it would created a corrupt internal 'segRows' field.
# 2011-10-20
# o Now append() for CBS also appends '...$params$knownSegments'.
# 2011-10-16
# o Added argument 'update' to mergeTwoSegments().
# 2011-10-10
# o Replaced extractRegions() with extractSegments() for CBS.
# o Added extractRegions() for CBS.
# 2011-10-08
# o Relabelled column 'id' to 'sampleName' returned by getSegments().
# o BUG FIX: getSegments() for CBS would not set 'id' for "splitter" rows.
# o Added mergeTwoSegments() for CBS.
# o Added updateMeans() for CBS.
# o Added all.equal() for CBS.
# 2011-10-02
# o CLEANUP: Moved getChromosomes(), nbrOfChromosomes(), nbrOfSegments(),
#   nbrOfLoci() and print() to AbstractCBS.
# o Now the CBS class extends the AbstractCBS class.
# 2011-09-04
# o Added writeSegments() for CBS.
# o Added writeLocusData() for CBS.
# o Added getSignalType() for CBS.
# o Added argument 'addCalls' to getLocusData().
# o Added getSampleName() for CBS.
# 2011-09-03
# o Added print() and as.character() for CBS.
# o Added CBS() constructor. Although it rairly will be used
#   it we be a place holder for the documentation.
# 2011-09-02
# o Added nbrOfLoci(), nbrOfSegments(), nbrOfChromosomes() and
#   getChromosomes() for CBS.
# 2010-11-19
# o Added append() for CBS objects.
############################################################################
