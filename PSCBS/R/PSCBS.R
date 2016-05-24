###########################################################################/**
# @RdocClass PSCBS
#
# @title "The PSCBS class"
#
# \description{
#  @classhierarchy
#
#  A PSCBS is an object containing results from parent-specific copy-number
#  (PSCN) segmentation.
# }
#
# \usage{PSCBS(fit=list(), ...)}
#
# \arguments{
#   \item{fit}{A @list structure containing the PSCN segmentation results.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#   @see "PairedPSCBS".
# }
#*/###########################################################################
setConstructorS3("PSCBS", function(fit=list(), ...) {
  # Argument 'fit':
  if (!is.list(fit)) {
    throw("Argument 'fit' is not a list: ", class(fit)[1]);
  }

  extend(AbstractCBS(fit, ...), "PSCBS");
})


setMethodS3("as.data.frame", "PSCBS", function(x, ...) {
  getSegments(x, splitter=TRUE, ...);
}, protected=TRUE)


setMethodS3("getLocusSignalNames", "PSCBS", function(fit, ...) {
  c("CT", "rho");
}, protected=TRUE)

setMethodS3("getSegmentTrackPrefixes", "PSCBS", function(fit, ...) {
  c("tcn", "dh");
}, protected=TRUE)


setMethodS3("getLocusData", "PSCBS", function(fit, indices=NULL, fields=c("asis"), ...) {
  # Argument 'indices':
  if (!is.null(indices)) {
    indices <- Arguments$getIndices(indices);
  }

  # Argument 'fields':
  fields <- match.arg(fields);

  data <- fit$data;

  # Return requested indices
  if (!is.null(indices)) {
    # Map of final indices to current indices
    map <- match(indices, data$index);

    # Extract/expand...
    data <- data[map,];

    # Sanity check
    stopifnot(nrow(data) == length(indices));
  }

  data;
}, protected=TRUE) # getLocusData()



setMethodS3("isSegmentSplitter", "PSCBS", function(fit, ...) {
  segs <- fit$output;

  isSplitter <- lapply(segs[-1], FUN=is.na);
  isSplitter <- Reduce("&", isSplitter);

  isSplitter;
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getSegments
#
# @title "Gets the segments"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{simplify}{If @TRUE, redundant and intermediate information is dropped.}#  \item{splitters}{If @TRUE, "splitters" between chromosomes are
#     preserved, otherwise dropped.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a SxK @data.frame, where S in the number of segments,
#   and K is the number of segment-specific fields.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getSegments", "PSCBS", function(fit, simplify=FALSE, splitters=TRUE, addGaps=FALSE, ...) {
  # Argument 'splitters':
  splitters <- Arguments$getLogical(splitters);

  segs <- fit$output;

  # Drop chromosome splitters?
  if (!splitters) {
    isSplitter <- isSegmentSplitter(fit);
    segs <- segs[!isSplitter,];
  }

  # Add splitters for "gaps"...
  if (splitters && addGaps) {
    # Chromosome gaps
    n <- nrow(segs);
    chrs <- segs$chromosome;
    gapsAfter <- which(diff(chrs) != 0L);
    gapsAfter <- gapsAfter[!is.na(chrs[gapsAfter])];
    nGaps <- length(gapsAfter);
    if (nGaps > 0L) {
      idxs <- seq(length=n);
      values <- rep(NA_integer_, times=nGaps);
      idxs <- insert(idxs, at=gapsAfter+1L, values=values);
      segs <- segs[idxs,];
    }

    # Other gaps
    n <- nrow(segs);
    chrs <- segs$chromosome;
    starts <- segs$tcnStart[-1L];
    ends <- segs$tcnEnd[-n];
    gapsAfter <- which(starts != ends);
    onSameChr <- (chrs[gapsAfter+1L] == chrs[gapsAfter] );
    gapsAfter <- gapsAfter[onSameChr];
    nGaps <- length(gapsAfter);
    if (nGaps > 0L) {
      idxs <- seq(length=n);
      values <- rep(NA_integer_, times=nGaps);
      idxs <- insert(idxs, at=gapsAfter+1L, values=values);
      segs <- segs[idxs,];
    }
  }

##  if (nrow(segs) > 0) {
##    segs$id <- getSampleName(fit);
##  }

  if (simplify) {
    # If joinSegments was used (i.e. (start,end) are equal for TCN and DH)...
    if (fit$params$joinSegments) {
      # Sanity check
      stopifnot(all(segs$tcnStart == segs$dhStart, na.rm=TRUE));
      stopifnot(all(segs$tcnEnd == segs$dhEnd, na.rm=TRUE));

      names <- colnames(segs);
      keep <- !is.element(names, c("dhStart", "dhEnd"));
      segs <- segs[,keep];
      names <- colnames(segs);
      names[names == "tcnStart"] <- "start";
      names[names == "tcnEnd"] <- "end";
      colnames(segs) <- names;
    }

    # Drop bootstrap columns, if any
    names <- colnames(segs);
    keep <- (regexpr("_[0-9]+(|[.][0-9]+)%$", names) == -1);
    segs <- segs[,keep];
  }

  segs;
}, private=TRUE)



setMethodS3("getChangePoints", "PSCBS", function(fit, ...) {
  # Already available?
  cps <- fit$changepoints;
  if (!is.null(cps)) return(cps);

  segs <- getSegments(fit, splitters=TRUE);
  tcn <- segs[["tcnMean"]];
  dh <- segs[["dhMean"]];
  C1 <- (1-dh) * tcn / 2;
  C2 <- tcn - C1;
  n <- length(tcn);

  # Calculate observed (alpha, radius, manhattan, dc1, dc2) data
  D1 <- C1[-n] - C1[-1L];
  D2 <- C2[-n] - C2[-1L];
  cps <- data.frame(
    alpha = atan2(D2, D1), # Changepoint angles in (0,2*pi)
    radius = sqrt(D2^2 + D1^2),
    manhattan = abs(D2) + abs(D1),
    d1 = D1,
    d2 = D2
  );

  cps;
}, private=TRUE) # getChangePoints()



############################################################################
# HISTORY:
# 2013-10-20
# o Added getChangePoints() for PSCBS.
# 2012-09-21
# o Now getSegments(..., splitters=TRUE) for CBS and PSCBS inserts NA
#   rows whereever there is a "gap" between segments.  A "gap" is when
#   two segments are not connected (zero distance).
# 2012-04-21
# o CLEANUP: Moved getSegmentSizes() from PSCBS to AbstractCBS.
# 2012-04-21
# o CLEANUP: Moved getSegmentSizes() from PairedPSCBS to PSCBS.
# 2012-02-27
# o Added argument 'fields' to getLocusData() for PairedPSCBS.
# 2011-12-12
# o Added optional argument 'indices' to getLocusData() to be able
#   to retrieve the locus-level data as indexed by input data.
# 2011-12-03
# o Added argument 'simplify' to getSegments().
# 2011-10-16
# o Added isSegmentSplitter().
# 2011-10-02
# o Now the CBS class extends the AbstractCBS class.
# o Added print() and as.data.frame() to PSCBS.
# o Added getSegments() to PSCBS.
# o DOCUMENTATION: Added Rdoc for several PSCBS methods.
# o Added a PSCBS constructor with documentation.
# 2010-12-01
# o Now also extractByChromosomes() and append() for PSCBS recognizes
#   fields 'tcnLociToExclude' and 'dhLociToExclude'.
# o BUG FIX: extractByChromosome() for PSCBS would call it self instead
#   of extractByChromosomes().
# 2010-11-26
# o Added extractByChromosomes() for PSCBS.
# 2010-09-26
# o getChromosomes() no longer returns NA divers.
# 2010-09-24
# o Added append() and more for PSCBS objects.
############################################################################
