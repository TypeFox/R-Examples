###########################################################################/**
# @RdocClass NonPairedPSCBS
#
# @title "The NonPairedPSCBS class"
#
# \description{
#  @classhierarchy
#
#  A NonPairedPSCBS is an object containing the results from the
#  Non-paired PSCBS method.
# }
#
# \usage{NonPairedPSCBS(fit=list(), ...)}
#
# \arguments{
#   \item{fit}{A @list structure containing the Non-paired PSCBS results.}
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
#   The @see "segmentByNonPairedPSCBS" method returns an object of this class.
# }
#*/###########################################################################
setConstructorS3("NonPairedPSCBS", function(fit=list(), ...) {
  # Argument 'fit':
  if (!is.list(fit)) {
    throw("Argument 'fit' is not a list: ", class(fit)[1]);
  }

  extend(PSCBS(fit=fit, ...), "NonPairedPSCBS");
})


setMethodS3("getLocusData", "NonPairedPSCBS", function(fit, ..., fields=c("asis", "full")) {
  # Argument 'fields':
  fields <- match.arg(fields);


  data <- NextMethod("getLocusData", fields="asis");


  if (fields == "full") {
    names <- colnames(data);

    data$isHet <- (data$muN == 1/2);

    # BACKWARD COMPATIBILITY: If 'rho' does not exists, calculate
    # it on the fly from 'betaT'.
    # NOTE: This should give an error in the future. /HB 2013-10-25
    if (is.null(data$rho)) {
      data$rho <- 2*abs(data$betaT-1/2);
      data$rho[!data$isHet] <- NA_real_;
      warning("Locus-level DH signals ('rho') did not exist and were calculated from tumor BAFs ('betaT')");
    }

    data$c1 <- 1/2*(1-data$rho)*data$CT;
    data$c2 <- data$CT - data$c1;

    data$isSNP <- (!is.na(data$betaT) | !is.na(data$muN));
    data$type <- ifelse(data$isSNP, "SNP", "non-polymorphic locus");

    # Labels
    data$muNx <- c("AA", "AB", "BB")[2*data$muN + 1L];
    data$isHetx <- c("AA|BB", "AB")[data$isHet + 1L];
  }

  data;
}, protected=TRUE) # getLocusData()


setMethodS3("callROH", "NonPairedPSCBS", function(fit, ...) {
  throw(sprintf("Cannot call ROH from '%s' data.", class(fit)[1L]));
}, private=TRUE) # callROH()


setMethodS3("updateMeans", "NonPairedPSCBS", function(fit, from=c("loci", "segments"), adjustFor=NULL, ..., avgTCN=c("asis", "mean", "median"), avgDH=c("asis", "mean", "median"), verbose=FALSE) {
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
      segs$dhMean[calls] <- 1/2;
      verbose && exit(verbose);
    }

    if (is.element("loh", adjustFor)) {
      verbose && enter(verbose, "Adjusting for LOH");
      calls <- segs$lohCall;
      segs$dhMean[calls] <- 0;
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



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update (C1,C2) mean levels
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Update (C1,C2) per segment");
  # Append (C1,C2) estimates
  tcn <- segs$tcnMean;
  dh <- segs$dhMean;
  C1 <- 1/2*(1-dh)*tcn;
  C2 <- tcn - C1;
  segs$c1Mean <- C1;
  segs$c2Mean <- C2;
  verbose && exit(verbose);


  # Return results
  res <- fit;
  res$output <- segs;
  res <- setMeanEstimators(res, tcn=avgTCN, dh=avgDH);

  verbose && exit(verbose);

  res;
}, private=TRUE) # updateMeans()



setMethodS3("resegment", "NonPairedPSCBS", function(fit, ..., verbose=FALSE) {
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

  # Use the locus-level data of the PairedPSCBS object
  data <- getLocusData(fit);
  class(data) <- "data.frame";
  drop <- c("rho", "betaTN", "index");
  keep <- !is.element(colnames(data), drop);
  data <- data[,keep];
  verbose && str(verbose, data);

  verbose && cat(verbose, "Number of loci: ", nrow(data));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup arguments to be passed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Overriding default arguments");
  segFcnName <- "segmentByPairedNonPSCBS";
  segFcn <- getMethodS3(segFcnName, "default");

  # (a) The default arguments
  formals <- formals(segFcn);

  formals <- formals[!sapply(formals, FUN=is.language)];
  formals <- formals[!sapply(formals, FUN=is.name)];
  drop <- c("chromosome", "x", "w", "CT", "betaT", "betaN", "muN", "...");
  keep <- !is.element(names(formals), drop);
  formals <- formals[keep];

  # (b) The arguments used in previous fit
  params <- fit$params;
  keep <- is.element(names(params), names(formals));
  params <- params[keep];
  # Don't trust 'tbn'!  TODO. /HB 20111117
  params$tbn <- NULL;

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

  verbose && enter(verbose, sprintf("Calling %s()", segFcnName));
  args <- append(list(data), args);
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, args[names(args) != "verbose"]);
  verbose && exit(verbose);

  fit <- do.call(segFcnName, args);
  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;
}, protected=TRUE) # resegment()



##############################################################################
# HISTORY
# 2013-04-23
# o BUG FIX: updateMeans() for PairedPSCBS and NonPairedPSCBS could
#   include a signal from a neighboring segment when averaging, iff
#   that signal was located at the exact locus of the change point.
# 2013-04-09
# o Added callROH() for NonPairedPSCBS that throws an exception.
# 2013-03-08
# o Added getLocusData() for NonPairedPSCBS.
# 2013-01-15
# o regsegment() was defined for PairedPSCBS instead of NonPairedPSCBS.
# 2012-04-21
# o Created from PairedPSCBS.R.
##############################################################################
