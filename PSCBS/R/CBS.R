###########################################################################/**
# @RdocClass CBS
#
# @title "The CBS class"
#
# \description{
#   A CBS object holds results from the
#   Circular Binary Segmentation (CBS) method
#   for \emph{one} sample for one or more chromosomes.
#
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to the constructor of @see "AbstractCBS".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Difference to DNAcopy object}{
#   A CBS object is similar to DNAcopy objects with the major
#   difference that a CBS object holds only one sample, whereas
#   a DNAcopy object can hold more than one sample.
# }
#
# \section{See also}{
#  The @see "segmentByCBS" method returns an object of this class.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("CBS", function(...) {
  extend(AbstractCBS(list(data=NULL, output=NULL), ...), "CBS");
})


setMethodS3("all.equal", "CBS", function(target, current, check.attributes=FALSE, ...) {
  # Compare class attributes
  res <- all.equal(class(target), class(current));
  if (!isTRUE(res)) {
    return(res);
  }

  # WORKAROUND: segmentByCBS() return getSegments(fit)$id without NA:s for
  # splitters, unless append() is used.
  # TO DO: Fix segmentByCBS() /HB 2011-10-08
  segs <- getSegments(target);
  if (nrow(segs) > 0) {
    isSplitter <- isSegmentSplitter(target);
    segs[isSplitter, "sampleName"] <- NA;
    target$output <- segs;
  }

  segs <- getSegments(current);
  if (nrow(segs) > 0) {
    isSplitter <- isSegmentSplitter(current);
    segs[isSplitter, "sampleName"] <- NA;
    current$output <- segs;
  }

  # NOTE: Here arguments 'target' and 'current' are lists and does not
  # have to be passed explicitly (although they have been modified).
  # If passed explicity, note that they must be named *and* that the
  # first/dispatch argument have to be passed as 'object=target'
  # (and never as 'target=target').  /HB 2014-02-03
  NextMethod("all.equal", object=target, current=current, check.attributes=check.attributes);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod as.data.frame
#
# @title "Gets the table of segments"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @data.frame, where each row corresponds to
#   a unique segment.
# }
#
# @author
#
# \seealso{
#   Utilizes @seemethod "getSegments".
#   @seeclass.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("as.character", "CBS", function(x, ...) {
  # To please R CMD check
  fit <- x;

  s <- sprintf("%s:", class(fit)[1]);

  s <- c(s, sprintf("Sample name: %s", getSampleName(fit)));

  s <- c(s, sprintf("Signal type: %s", getSignalType(fit)));

  s <- c(s, sprintf("Number of segments: %d", nbrOfSegments(fit)));

  s <- c(s, sprintf("Number of loci: %d", nbrOfLoci(fit)));

  n <- getSegments(fit)$nbrOfLoci;
  q <- quantile(n, probs=c(0.00, 0.05, 0.25, 0.50, 0.75, 0.95, 1.00), na.rm=TRUE);
  qs <- sprintf("%g [%s]", q, names(q));
  s <- c(s, sprintf("Number of loci per segment: %s", paste(qs, collapse=", ")));

  chrs <- getChromosomes(fit);
  s <- c(s, sprintf("Chromosomes: [%d] %s", length(chrs), hpaste(chrs)));

  s <- c(s, sprintf("Standard deviation: %g", estimateStandardDeviation(fit)));

  tt <- grep("Call$", colnames(getLocusData(fit)), value=TRUE);
  s <- c(s, sprintf("Locus calls: [%d] %s", length(tt), hpaste(tt)));

  segs <- getSegments(fit);
  callCols <- grep("Call$", colnames(segs), value=TRUE);
  callTypes <- gsub("Call$", "", callCols);
  s <- c(s, sprintf("Types of segment calls: [%d] %s", length(callTypes), hpaste(callTypes)));
  for (kk in seq(along=callCols)) {
    key <- callCols[kk];
    type <- callTypes[kk];
    n <- sum(segs[,key], na.rm=TRUE);
    if (type == "loss") {
      nC <- sum(isWholeChromosomeLost(fit));
    } else if (type == "gain") {
      nC <- sum(isWholeChromosomeGained(fit));
    } else {
      nC <- NA;
    }
    s <- c(s, sprintf("Number of chromosomes (segments) called '%s': %d (%d)", type, nC, n));
  }

  GenericSummary(s);
}, protected=TRUE)


setMethodS3("as.data.frame", "CBS", function(x, ...) {
  getSegments(x, splitter=FALSE, ...);
}, protected=TRUE)


setMethodS3("getSignalType", "CBS", function(fit, ...) {
  type <- fit$signalType;
  if (is.null(type)) type <- as.character(NA);
  type;
}, protected=TRUE)


setMethodS3("signalType", "CBS", function(fit, ...) {
  getSignalType(fit);
}, protected=TRUE)


"signalType<-" <- function(x, value) {
  UseMethod("signalType<-");
}

setMethodS3("signalType<-", "CBS", function(x, value) {
  fit <- x;

  # Argument 'value':
  value <- Arguments$getCharacter(value);

  fit$signalType <- value;
  fit;
}, private=TRUE, addVarArgs=FALSE)



setMethodS3("getLocusSignalNames", "CBS", function(fit, ...) {
  data <- fit$data;
  names <- colnames(data);
  if (is.element("y", names)) {
    return("y");
  } else if (is.element("CT", names)) {
    return("CT");
  }

  throw("INTERNAL ERROR: Unknown locus signal names: ", paste(names, collapse=", "));
}, protected=TRUE)

setMethodS3("getSegmentTrackPrefixes", "CBS", function(fit, ...) {
  c("");
}, protected=TRUE)


setMethodS3("getLocusData", "CBS", function(fit, indices=NULL, addCalls=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'indices':
  if (!is.null(indices)) {
    indices <- Arguments$getIndices(indices);
  }

  # Argument 'addCalls':
  if (is.logical(addCalls)) {
    addCalls <- Arguments$getLogical(addCalls);
    if (!addCalls) {
      addCalls <- NULL;
    }
  } else {
    addCalls <- Arguments$getCharacters(addCalls);
  }

  data <- fit$data;

  # Append segment calls?
  if (length(addCalls) > 0) {
    callsL <- extractCallsByLocus(fit);
    if (is.character(addCalls)) {
      callsL <- callsL[,addCalls];
    }

    # Sanity check
    stopifnot(nrow(callsL) == nrow(data));

    data <- cbind(data, callsL);
  }

  # Return requested indices
  if (!is.null(indices)) {
    # Map of final indices to current indices
    map <- match(indices, data$index);

    # Extract/expand...
    data <- data[map,];
    rownames(data) <- NULL;

    # Sanity check
    stopifnot(nrow(data) == length(indices));
  }

  data;
}, private=TRUE) # getLocusData()


setMethodS3("isSegmentSplitter", "CBS", function(fit, ...) {
  segs <- fit$output;

  isSplitter <- lapply(segs[-1], FUN=is.na);
  isSplitter <- Reduce("&", isSplitter);

  isSplitter;
}, protected=TRUE)


setMethodS3("getSegments", "CBS", function(fit, simplify=FALSE, splitters=TRUE, addGaps=FALSE, ...) {
  # Argument 'splitters':
  splitters <- Arguments$getLogical(splitters);

  segs <- fit$output;

  isSplitter <- isSegmentSplitter(fit);

  # Add 'sampleName' column?
  if (nrow(segs) > 0) {
    sampleName <- rep(getSampleName(fit), times=nrow(segs));
    sampleName[isSplitter] <- as.character(NA);
    if (!is.element("sampleName", colnames(segs))) {
      segs <- cbind(sampleName=I(sampleName), segs);
    } else {
      segs[,"sampleName"] <- sampleName;
    }
  }

  # Drop chromosome splitters?
  if (!splitters) {
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

  segs;
}, private=TRUE)



setMethodS3("getChangePoints", "CBS", function(fit, ...) {
  # Already available?
  cps <- fit$changepoints;
  if (!is.null(cps)) return(cps);

  segs <- getSegments(fit, splitters=TRUE);
  tcn <- segs[["mean"]];
  n <- length(tcn);

  # Calculate observed (d) data
  D <- tcn[-n] - tcn[-1L];
  cps <- data.frame(
    d = D
  );

  cps;
}, private=TRUE) # getChangePoints()



setMethodS3("updateBoundaries", "CBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Updating boundaries");
  verbose && cat(verbose, "Number of segments: ",
                                  nbrOfSegments(fit, splitters=FALSE));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit);
  segs <- getSegments(fit, splitters=TRUE);
  segRows <- fit$segRows;

  nbrOfSegments <- nrow(segs);
  chromosome <- data$chromosome;
  x <- data$x;
  y <- data$y;
  w <- data$w;
  hasWeights <- !is.null(w);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (ss in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment %d of %d", ss, nbrOfSegments));
    segRow <- segRows[ss,];
    seg <- segs[ss,];

    # A splitter - nothing todo?
    if (is.na(segRow[[1]]) && is.na(segRow[[2]])) {
      next;
    }

    # (a) Identify units (loci)
    units <- segRow[[1]]:segRow[[2]];
    verbose && cat(verbose, "Loci:");
    verbose && str(verbose, units);

    # (b) Extract signals
    ySS <- y[units];
    xSS <- x[units];
    cSS <- chromosome[units];
    if (hasWeights) {
      wSS <- w[units];
    }

    # (c) Drop missing values
    keep <- (!is.na(ySS) & !is.na(xSS) & !is.na(cSS));
    if (hasWeights) {
      keep <- keep & (!is.na(wSS) & wSS > 0);
    }
    keep <- which(keep);
    ySS <- ySS[keep];
    xSS <- xSS[keep];
    cSS <- cSS[keep];
    if (hasWeights) {
      wSS <- wSS[keep];
    }
    units <- units[keep];
    verbose && cat(verbose, "Loci (non-missing):");
    verbose && str(verbose, units);

    # (d) Identify (chromosome, start, stop)
    stopifnot(all(cSS == cSS[1]));
    cSS <- cSS[1];
    xRange <- range(xSS, na.rm=TRUE);
    verbose && cat(verbose, "Range:");
    verbose && print(verbose, xRange);

    # (e) Update segment information
    seg$chromosome <- cSS;
    seg$start <- xRange[1];
    seg$end <- xRange[2];

    segs[ss,] <- seg;

    verbose && exit(verbose);
  } # for (ss ...)

  # Update results
  res <- fit;
  res$output <- segs;

  # Rejoin segments?
  if (isTRUE(res$params$joinSegments)) {
    res <- joinSegments(res, verbose=less(verbose,10));
  }

  verbose && exit(verbose);

  res;
}, protected=TRUE) # updateBoundaries()



setMethodS3("updateMeans", "CBS", function(fit, ..., avg=c("asis", "mean", "median"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'avg':
  avg <- match.arg(avg);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Updating mean level estimates");
  verbose && cat(verbose, "Number of segments: ",
                                  nbrOfSegments(fit, splitters=FALSE));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit);
  segs <- getSegments(fit, splitters=TRUE);
  segRows <- fit$segRows;

  nbrOfSegments <- nrow(segs);
  chromosome <- data$chromosome;
  x <- data$x;
  y <- data$y;
  w <- data$w;
  hasWeights <- !is.null(w);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up averaging functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (avg == "asis") {
    est <- fit$params$meanEstimators;
    avg <- est$y;
    if (is.null(avg)) avg <- "mean";
    avg <- match.arg(avg);
  }

  if (hasWeights) {
    if(avg == "mean") {
      avgFUN <- weighted.mean;
    } else if(avg == "median") {
      avgFUN <- weightedMedian;
    } else {
      throw("Value of argument 'avg' is not supported with weights: ", avg);
    }
  } else {
    avgFUN <- get(avg, mode="function");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (ss in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment %d of %d", ss, nbrOfSegments));
    segRow <- segRows[ss,];
    seg <- segs[ss,];

    # A splitter - nothing todo?
    if (!is.finite(segRow[[1]]) || !is.finite(segRow[[2]])) {
      next;
    }

    # (a) Identify units (loci)
    units <- segRow[[1]]:segRow[[2]];

    # (b) Extract signals
    ySS <- y[units];
    if (hasWeights) {
      wSS <- w[units];
    }

    # (c) Drop missing values
    keep <- (!is.na(ySS));
    if (hasWeights) {
      keep <- keep & (!is.na(wSS) & wSS > 0);
    }
    keep <- which(keep);
    ySS <- ySS[keep];
    if (hasWeights) {
      wSS <- wSS[keep];
    }
    units <- units[keep];
    nbrOfLoci <- length(units);

    # (d) Update mean
    if (hasWeights) {
      wSS <- wSS / sum(wSS);
      gamma <- avgFUN(ySS, w=wSS);
    } else {
      gamma <- avgFUN(ySS);
    }

    # Sanity check
    stopifnot(nbrOfLoci == 0 || !is.na(gamma));

    # (d) Update the segment statistics
    seg$mean <- gamma;
    seg$nbrOfLoci <- nbrOfLoci;

    segs[ss,] <- seg;

    verbose && exit(verbose);
  } # for (ss ...)

  # Return results
  res <- fit;
  res$output <- segs;
  res <- setMeanEstimators(res, y=avg);

  verbose && exit(verbose);

  res;
}, protected=TRUE) # updateMeans()


setMethodS3("resegment", "CBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Resegmenting a ", class(fit)[1], " object");
  segFcnName <- "segmentByCBS";
  segFcn <- getMethodS3(segFcnName, "default");

  # Use the locus-level data of the segmentation object
  data <- getLocusData(fit);
  class(data) <- "data.frame";
  drop <- c("index");
  keep <- !is.element(colnames(data), drop);
  data <- data[,keep];
  verbose && str(verbose, data);

  verbose && cat(verbose, "Number of loci: ", nrow(data));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup arguments to be passed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Overriding default arguments");

  # (a) The default arguments
  formals <- formals(segFcn);

  formals <- formals[!sapply(formals, FUN=is.language)];
  formals <- formals[!sapply(formals, FUN=is.name)];
  drop <- c("chromosome", "x", "y", "w", "...");
  keep <- !is.element(names(formals), drop);
  formals <- formals[keep];

  # (b) The arguments used in previous fit
  params <- fit$params;
  keep <- is.element(names(params), names(formals));
  params <- params[keep];

  # (c) The arguments in '...'
  userArgs <- list(..., verbose=verbose);

  # (d) Merge
  args <- formals;
  args2 <- append(params, userArgs);
  for (kk in seq(along=args2)) {
    value <- args2[[kk]];
    if (!is.null(value)) {
      key <- names(args2)[kk];
      if (!is.null(key)) {
        args[[key]] <- value;
      } else {
        args <- append(args, list(value));
      }
    }
  } # for (key ...)
  verbose && str(verbose, args[names(args) != "verbose"]);

  args <- append(list(data), args);
  verbose && cat(verbose, "Arguments with data:");
  verbose && str(verbose, args[names(args) != "verbose"]);
  verbose && exit(verbose);

  verbose && enter(verbose, sprintf("Calling %s()", segFcnName));
  fit <- do.call(segFcnName, args);
  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;
}, protected=TRUE) # resegment()


############################################################################
# HISTORY:
# 2014-02-03
# o BUG FIX: all.equal() for CBS would pass the first/dispatch argument
#   to NextMethod() as 'target=target' and not as 'object=target', which
#   would result in it being passed it twice both named and non-named
#   where the latter would become argument 'tolerance=target' in an
#   internal call to all.equal() for numerics.  In recent R-devel version
#   this would generate "Error in all.equal.numeric(target[[i]],
#   current[[i]], check.attributes = check.attributes, : 'tolerance'
#   should be numeric  Calls: stopifnot ... all.equal.default ->
#   all.equal.list -> all.equal -> all.equal.numeric".
# 2013-12-17
# o BUG FIX: getChangePoints() for CBS returned empty results.
# 2013-10-20
# o Added getChangePoints() for CBS.
# 2012-09-21
# o Now getSegments(..., splitters=TRUE) for CBS and PSCBS inserts NA
#   rows whereever there is a "gap" between segments.  A "gap" is when
#   two segments are not connected (zero distance).
# 2012-06-03
# o BUG FIX: all.equal(target, current) for CBS objects would give an
#   error if either 'target' or 'current' had zero segments.
# 2011-12-12
# o Added optional argument 'indices' to getLocusData() to be able
#   to retrieve the locus-level data as indexed by input data.
# 2011-11-17
# o Added resegment() for CBS for easy resegmentation.
# 2011-11-15
# o Now updateMeans() uses locus-specific weights, iff available.
# o Added updateBoundaries() for CBS to update (start,stop) per segment.
# o CORRECTNESS: Now updateMeans() for CBS identify loci by the internal
#   'segRows' field and no longer by locations of segment boundaries,
#   which gave slightly incorrect estimates for "tied" loci.
# 2011-10-16
# o Added isSegmentSplitter().
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
