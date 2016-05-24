###########################################################################/**
# @set class=DNAcopy
# @RdocMethod as.CBS
#
# @title "Coerces a DNAcopy object to a CBS object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A @see "DNAcopy" object (of the \pkg{DNAcopy} package.)}
#   \item{sample}{An index specifying which sample to extract,
#     if more than one exists.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "CBS" object.
# }
#
# @author "HB"
#
# \seealso{
#   \code{\link[PSCBS:as.DNAcopy.CBS]{as.DNAcopy()}}.
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("as.CBS", "DNAcopy", function(fit, sample=1L, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- fit$data;
  sample <- Arguments$getIndex(sample, max=ncol(data)-2L);

  sampleName <- colnames(data)[sample+2L];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the 'data' field
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- fit$data;
  rownames <- rownames(data);

  data <- data.frame(
    chromosome = data$chrom,
    x          = data$maploc,
    y          = data[,sample+2L,drop=TRUE],
    stringsAsFactors=FALSE
  );
  rownames(data) <- rownames;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the 'output' field
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  output <- fit$output;
  ID <- NULL; rm(list="ID"); # To please R CMD check
  output <- subset(output, ID == sampleName);
  rownames <- rownames(output);

  output <- data.frame(
    chromosome = output$chrom,
    start      = output$loc.start,
    end        = output$loc.end,
    nbrOfLoci  = as.integer(output$num.mark),
    mean       = output$seg.mean,
    stringsAsFactors=FALSE
  );
  rownames(output) <- rownames;

  # Add chromosome splitter
  ats <- which(diff(output$chromosome) != 0) + 1L;
  if (length(ats) > 0) {
    idxs <- seq(length=nrow(output));
    values <- rep(NA_integer_, times=length(ats));
    expand <- insert(idxs, ats=ats, values=values); # R.utils::insert()
    output <- output[expand,];
    rownames(output) <- NULL;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup up 'CBS' object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (sampleName == "<NA>") sampleName <- as.character(NA);
  res <- list();
  res$sampleName <- sampleName;
  res$data <- data;
  res$output <- output;
  res$params <- list();
  class(res) <- class(CBS());

  res;
}) # as.CBS()


setMethodS3("extractTotalCNs", "CBS", function(fit, ...) {
  data <- getSegments(fit, ...);
  data[,c("mean", "nbrOfLoci"), drop=FALSE];
}, protected=TRUE)


setMethodS3("extractCNs", "CBS", function(fit, ...) {
  data <- extractTotalCNs(fit, ...);
  data <- data[,c("mean"), drop=FALSE];
  data <- as.matrix(data);
  data;
}, protected=TRUE)



setMethodS3("extractChromosomes", "CBS", function(x, chromosomes, ...) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chromosomes':
  disallow <- c("NaN", "Inf");
  chromosomes <- Arguments$getIntegers(chromosomes, range=c(0,Inf), disallow=disallow);
  stopifnot(all(is.element(chromosomes, getChromosomes(this))));

  # Always extract in order
  chromosomes <- unique(chromosomes);
  chromosomes <- sort(chromosomes);

  # Allocate results
  res <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locus data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chromosome <- NULL; rm(list="chromosome"); # To please R CMD check
  data <- getLocusData(this);
  class <- class(data);
  class(data) <- "data.frame";
  data <- subset(data, chromosome %in% chromosomes);
  class(data) <- class;
  res$data <- data;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segmentation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify rows to subset
  rows <- which(is.element(res$output$chromosome, chromosomes));
  for (field in c("output", "segRows")) {
    res[[field]] <- res[[field]][rows,,drop=FALSE];
  }

  # Identify chromosome offsets
  data <- getLocusData(this);
  chrStarts <- match(getChromosomes(this), data$chromosome);
  chrEnds <- c(chrStarts[-1]-1L, nrow(data));
  chrLengths <- chrEnds - chrStarts + 1L;

  chrLengthsExcl <- chrLengths;

  keep <- match(chromosomes, getChromosomes(this));
  chrLengthsExcl[keep] <- 0L;
  cumChrLengthsExcl <- cumsum(chrLengthsExcl);

  shifts <- cumChrLengthsExcl[keep];
  stopifnot(all(is.finite(shifts)));

  # Adjust indices
  for (cc in seq(along=chromosomes)) {
    chromosome <- chromosomes[cc];
    shift <- shifts[cc];
    # Nothing to do?
    if (shift == 0) next;
    for (field in c("segRows")) {
      segRows <- res[[field]];
      rows <- which(res$output$chromosome == chromosome);
      segRows[rows,] <- segRows[rows,] - shift;
      res[[field]] <- segRows;
    }
  }

  res;
}, protected=TRUE)


setMethodS3("subset", "CBS", function(x, chromlist=NULL, ...) {
  extractChromosomes(x, chromosomes=chromlist, ...);
}, private=TRUE)



###########################################################################/**
# @set "class=CBS"
# @RdocMethod extractSegmentMeansByLocus
#
# @title "Extracts segments means at each locus"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getLocusData".}
# }
#
# \value{
#  Returns a @numeric @vector of length \code{nbrOfLoci()}.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("extractSegmentMeansByLocus", "CBS", function(fit, ...) {
  data <- getLocusData(fit, ...);
  chromosome <- data$chromosome;
  x <- data$x;
  y <- data[,3];

  segs <- getSegments(fit);
  nbrOfSegments <- nrow(segs);
  nbrOfLoci <- nrow(data);

  # Get mean estimators
  estList <- getMeanEstimators(fit, "y");
  avgY <- estList$y;

  yS <- y;
  for (ss in seq(length=nbrOfSegments)) {
    seg <- segs[ss,];
    idxs <- which(seg$chromosome == chromosome &
                  seg$start <= x & x <= seg$end);
    idxs <- Arguments$getIndices(idxs, max=nbrOfLoci);

    ySS <- y[idxs];
    ok <- is.finite(ySS);

    # Sanity check
    ## stopifnot(sum(ok) == seg$nbrOfLoci); # Not dealing with ties

    mu <- avgY(ySS[ok]);
    yS[idxs] <- mu;
  } # for (ss ...)

  yS;
}, private=TRUE) # extractSegmentMeansByLocus()



###########################################################################/**
# @set "class=CBS"
# @RdocMethod estimateStandardDeviation
#
# @title "Estimates the whole-genome standard deviation of the signals"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{chromosomes}{An optional @vector specifying the subset of
#    chromosomes used for the estimate.  If @NULL, all chromosomes are used.}
#  \item{method}{A @character string specifying the method used.}
#  \item{estimator}{A @character string or a @function specifying the
#    internal estimator.}
#  \item{na.rm}{If @TRUE, missing values are dropped, otherwise not.}
#  \item{weights}{An optional @double @vector of \code{nbrOfLoci()}
#    non-negative weights.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a non-negative @numeric scale.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("estimateStandardDeviation", "CBS", function(fit, chromosomes=NULL, method=c("diff", "res", "abs", "DNAcopy"), estimator=c("mad", "sd"), na.rm=TRUE, weights=NULL, ...) {
  # Local copies of DNAcopy functions
  DNAcopy_trimmed.variance <- .use("trimmed.variance", package="DNAcopy");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chromosomes':
  if (!is.null(chromosomes)) {
  }

  # Argument 'method':
  method <- match.arg(method);

  # Argument 'estimator':
  estimator <- match.arg(estimator);


  # Argument 'weights':
  if (!is.null(weights)) {
    nbrOfLoci <- nbrOfLoci(fit);
    weights <- Arguments$getNumerics(weights, range=c(0,Inf),
                                     length=rep(nbrOfLoci, times=2));
  }


  # Get the estimator function
  if (!is.null(weights)) {
    estimator <- sprintf("weighted %s", estimator);
    estimator <- R.utils::toCamelCase(estimator);
  }

  if (method == "DNAcopy") {
    estimatorFcn <- function(y, trim=0.025, ...) {
      sigma2 <- DNAcopy_trimmed.variance(y, trim=trim);
      sqrt(sigma2);
    }
  } else {
    estimatorFcn <- get(estimator, mode="function");
  }


  # Subset by chromosomes?
  if (!is.null(chromosomes)) {
    fit <- extractChromosomes(fit, chromosomes=chromosomes);
  }

  nbrOfLoci <- nbrOfLoci(fit);
  # Nothing to do?
  if (nbrOfLoci <= 1) {
    sigma <- NA_real_;
    attr(sigma, "nbrOfLoci") <- nbrOfLoci;
    attr(sigma, "df") <- NA_integer_;
    return(sigma);
  }

  data <- getLocusData(fit);
  y <- data[,3];

  if (method == "diff") {
    dy <- diff(y);
    # Weighted estimator?
    if (!is.null(weights)) {
      # Calculate weights per pair
      weights <- (weights[1:(nbrOfLoci-1)]+weights[2:nbrOfLoci])/2;
      sigma <- estimatorFcn(dy, w=weights, na.rm=na.rm)/sqrt(2);
    } else {
      sigma <- estimatorFcn(dy, na.rm=na.rm)/sqrt(2);
    }
    df <- length(dy);
  } else if (method == "res") {
    yS <- extractSegmentMeansByLocus(fit);
    dy <- y - yS;
    if (!is.null(weights)) {
      sigma <- estimatorFcn(dy, w=weights, na.rm=na.rm);
    } else {
      sigma <- estimatorFcn(dy, na.rm=na.rm);
    }
    df <- length(dy);
  } else if (method == "abs") {
    if (!is.null(weights)) {
      sigma <- estimatorFcn(y, w=weights, na.rm=na.rm);
    } else {
      sigma <- estimatorFcn(y, na.rm=na.rm);
    }
    df <- length(y);
  } else if (method == "DNAcopy") {
    if (na.rm) {
      y <- y[!is.na(y)];
    }
    sigma <- estimatorFcn(y, ...);
    df <- length(y);
  } else {
    throw("Method no implemented: ", method);
  }

  attr(sigma, "nbrOfLoci") <- nbrOfLoci;
  attr(sigma, "df") <- df;

  sigma;
}) # estimateStandardDeviation()



setMethodS3("getChromosomeRanges", "CBS", function(fit, ...) {
  # To please R CMD check, cf. subset()
  chromosome <- NULL; rm(list="chromosome");

  segs <- getSegments(fit, splitter=FALSE);
  chromosomes <- sort(unique(segs$chromosome));

  # Allocate
  naValue <- NA_real_;
  res <- matrix(naValue, nrow=length(chromosomes), ncol=3);
  rownames(res) <- chromosomes;
  colnames(res) <- c("start", "end", "length");

  # Get start and end of each chromosome.
  for (ii in seq(length=nrow(res))) {
    chr <- chromosomes[ii];
    segsII <- subset(segs, chromosome == chr);
    res[ii,"start"] <- min(segsII$start, na.rm=TRUE);
    res[ii,"end"] <- max(segsII$end, na.rm=TRUE);
  } # for (ii ...)

  res[,"length"] <- res[,"end"] - res[,"start"] + 1L;

  # Sanity check
  stopifnot(nrow(res) == length(chromosomes));

  res <- as.data.frame(res);
  res <- cbind(chromosome=chromosomes, res);

  res;
}, protected=TRUE) # getChromosomeRanges()



############################################################################
# HISTORY:
# 2012-06-03
# o BUG FIT: The recent updates of as.CBS() for DNAcopy would not work
#   for samples with name '<NA>'.
# 2012-05-30
# o BUG FIX: as.CNA() for DNAcopy added incorrect chromosome splitters.
# o BUG FIX: as.CNA() for DNAcopy would ignore argument 'sample' and
#   always return the first sample.
# 2011-12-12
# o Now extractSegmentMeansByLocus() for CBS passes arguments
#   '...' to getLocusData().
# 2011-11-28
# o extractCNs() for CBS would not return a matrix but a data.frame.
# o BUG FIX: extractTotalCNs() for CBS would give an error.
# 2011-11-15
# o Added method="DNAcopy" to estimateStandardDeviation() for CBS, which
#   estimates the std. dev. using DNAcopy:::trimmed.variance().
# 2011-10-16
# o Added extractTotalCNs() for CBS.
# o Implemented extractCNs() for CBS.
# 2011-10-08
# o BUG FIX: The object returned by as.CBS() of DNAcopy did not have the
#   correct class hierarchy.
# 2011-10-06
# o Now getChromosomeRanges() of CBS returns a data.frame instead of
#   a matrix, and first column is now 'chromosome'.
# 2011-09-05
# o Added getChromosomeRanges() for CBS.
# 2011-09-04
# o Added estimateStandardDeviation() for CBS.
# o Added extractSegmentMeansByLocus() for CBS.
# 2011-09-03
# o Added as.CBS() for DNAcopy to coerce a DNAcopy object to a CBS object.
# 2011-09-02
# o Added extractByChromosomes() for CBS.
# o Added subset() for CBS for backward compatibility.
# o Added nbrOfLoci(), getChromosomes() and getSampleNames() for CBS.
# 2010-11-19
# o Added append() for CBS objects.
############################################################################
