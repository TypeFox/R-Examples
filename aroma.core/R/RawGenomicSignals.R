###########################################################################/**
# @RdocClass RawGenomicSignals
#
# @title "The RawGenomicSignals class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{y}{A @numeric @vector of length J specifying the signal
#     at each locus.}
#   \item{x}{A (optional) @numeric @vector of length J specifying the
#     position of each locus.}
#   \item{w}{A (optional) non-negative @numeric @vector of length J
#     specifying a weight of each locus.}
#   \item{chromosome}{An (optional) @integer specifying the chromosome for
#     these genomic signals.}
#   \item{name}{An (optional) @character string specifying the sample name.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("RawGenomicSignals", function(y=NULL, x=NULL, w=NULL, chromosome=0L, name=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  object <- NULL;
  if (!is.null(y)) {
    if (inherits(y, "RawGenomicSignals")) {
      object <- y;
      y <- getSignals(object);
      x <- object$x;
      w <- object$w;
      chromosome <- object$chromosome;
      name <- getName(object);
    }

    if (!is.vector(y)) {
      throw("Argument 'y' must be a vector: ", mode(y)[1]);
    }

    if (!is.numeric(y)) {
      throw("Argument 'y' must be a numeric: ", class(y)[1]);
    }
  }
  n <- length(y);

  # Argument 'chromosome':
  if (is.null(chromosome)) {
    chromosome <- as.integer(NA);
  }
  if (length(chromosome) == 1) {
    chromosome <- rep(chromosome, times=n);
  }
  chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), length=c(n,n), disallow=c("NaN", "Inf"));

  # Argument 'x':
  if (!is.null(x)) {
    x <- Arguments$getNumerics(x, length=c(n,n));
  }

  # Argument 'w':
  if (!is.null(w)) {
    w <- Arguments$getNumerics(w, range=c(0,Inf), length=c(n,n));
  }

  # Arguments '...':
  args <- list(...);
  if (length(args) > 0) {
    argsStr <- paste(names(args), collapse=", ");
    throw("Unknown arguments: ", argsStr);
  }

  # Setup data frame with optional fields
  data <- RichDataFrame(chromosome=chromosome, .name=name);
  data$x <- x;  # optional
  data$y <- y;
  data$w <- w;  # optional

  this <- extend(data, "RawGenomicSignals");

  # Append other locus fields?
  if (!is.null(object)) {
    fields <- setdiff(colnames(object), colnames(this));
    for (field in fields) {
      this[[field]] <- object[[field]];
    }
  }

  this;
})


setMethodS3("as.character", "RawGenomicSignals", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- getGenericSummary(this);
  chrs <- getChromosomes(this);
  nbrOfChrs <- length(chrs);
  s <- c(s, sprintf("Chromosomes: %s [%d]", seqToHumanReadable(chrs), nbrOfChrs));
  n <- nbrOfLoci(this);
  s <- c(s, sprintf("Number of loci: %d", n));
  if (nbrOfChrs == 1) {
    xRange <- xRange(this);
    s <- c(s, sprintf("Position range: [%g,%g]", xRange[1], xRange[2]));
    dAvg <- if (n >= 2) diff(xRange)/(n-1) else as.double(NA);
    s <- c(s, sprintf("Mean distance between loci: %g", dAvg));
  }

  s;
}, protected=TRUE)



setMethodS3("print", "RawGenomicSignals", function(x, ...) {
  print(as.character(x));
}, protected=TRUE)


setMethodS3("getSignalColumnNames", "RawGenomicSignals", function(this, translate=TRUE, ...) {
  names <- "y";
  if (translate) {
    names <- translateColumnNames(this, names);
  }
  names;
}, protected=TRUE)

setMethodS3("getSignalColumnName", "RawGenomicSignals", function(this, ...) {
  getSignalColumnNames(this, ...)[1L];
}, protected=TRUE)


setMethodS3("nbrOfLoci", "RawGenomicSignals", function(this, na.rm=FALSE, ...) {
  if (!na.rm) {
    return(nrow(this));
  }

  y <- getSignals(this, ...);
  y <- y[is.finite(y)];
  length(y);
})


setMethodS3("getChromosomes", "RawGenomicSignals", function(this, ...) {
  chrs <- this$chromosome;
  if (is.null(chrs)) {
    chrs <- as.integer(NA);
  }
  chrs <- as.integer(chrs);
  chrs <- unique(chrs);
  chrs <- sort(chrs);
  chrs;
})

setMethodS3("nbrOfChromosomes", "RawGenomicSignals", function(this, ...) {
  chrs <- getChromosomes(this, ...);
  length(chrs);
})


setMethodS3("assertOneChromosome", "RawGenomicSignals", function(this, ...) {
  if (nbrOfChromosomes(this) > 1) {
    throw(sprintf("Cannot perform operation. %s has more than one chromosome: %s", class(this)[1], nbrOfChromosomes(this)));
  }
}, protected=TRUE)




# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# EXTRACT METHODS
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
setMethodS3("extractSubset", "RawGenomicSignals", function(this, subset, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'subset':
  n <- nbrOfLoci(this);
  if (is.logical(subset)) {
    subset <- Arguments$getLogicals(subset, length=c(n,n));
    subset <- which(subset);
  } else {
    subset <- Arguments$getIndices(subset, max=n);
  }

  res <- this[subset,,drop=FALSE];

  # Sanity check
  stopifnot(nbrOfLoci(res) == length(subset));

  res;
}, protected=TRUE) # extractSubset()



setMethodS3("extractRegion", "RawGenomicSignals", function(this, region, chromosome=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chromosome':

  # Argument 'region':
  region <- Arguments$getNumerics(region, length=c(2,2));
  stopifnot(region[1] <= region[2]);


  # Subset by chromosome
  if (!is.null(chromosome)) {
    rgs <- extractChromosome(this, chromosome=chromosome);
  } else {
    rgs <- this;
  }

  # This is a single-chromosome method. Assert that's the case.
  assertOneChromosome(rgs);

  x <- getPositions(rgs);
  keep <- which(region[1] <= x & x <= region[2]);

  extractSubset(rgs, keep);
}, protected=TRUE) # extractRegion()



setMethodS3("extractRegions", "RawGenomicSignals", function(this, regions, chromosomes=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'regions':
  cl <- class(regions)[1];
  if (inherits(regions, "CopyNumberRegions")) {
    regions <- as.data.frame(regions);
    regions <- regions[,c("chromosome", "start", "stop")];
  } else if (is.matrix(regions)) {
    regions <- as.data.frame(regions);
  }
  if (!is.data.frame(regions)) {
    throw("Argument 'regions' is neither a CopyNumberRegions object, a data.frame, nor a matrix: ", cl);
  }

  # Argument 'chromosomes':
  if (is.null(chromosomes) && !is.element("chromosome", colnames(regions))) {
    # Backward compatibility only
    # This is a single-chromosome method. Assert that is the case.
    assertOneChromosome(this);
    chromosomes <- getChromosomes(this);
    chromosomes <- rep(chromosomes, times=nrow(regions));
  }

  # Argument 'regions' (again):
  if (!is.null(chromosomes)) {
    if (is.element("chromosome", names(regions))) {
      throw("Argument 'chromosomes' must not be specified when argument 'regions' already specifies chromosomes: ", hpaste(colnames(regions)));
    }
    regions <- cbind(data.frame(chromosome=chromosomes), regions);
    chromosomes <- NULL; # Not needed anymore
  }

  reqs <- c("chromosome", "start", "stop");
  missing <- reqs[!is.element(reqs, colnames(regions))];
  if (length(missing)) {
    throw("Missing fields in argument 'regions': ", hpaste(missing));
  }

  # Extract fields of interest
  regions <- regions[,reqs];

  # Assert the fields are numeric
  for (key in colnames(regions)) {
    regions[[key]] <- Arguments$getNumerics(regions[[key]], .name=key);
  }

  # Assert ordered (start,stop)
  stopifnot(all(regions[["start"]] <= regions[["stop"]]));



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For each chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chromosome <- this$chromosome;
  x <- getPositions(this);
  keep <- rep(FALSE, times=length(x));
  for (rr in seq_len(nrow(regions))) {
    region <- unlist(regions[rr,], use.names=TRUE);
    chr <- region["chromosome"];
    start <- region["start"];
    stop <- region["stop"];
    keepRR <- (chromosome == chr & start <= x & x <= stop);
    keep <- keep | keepRR;
  } # for (rr ...)

  keep <- which(keep);

  extractSubset(this, keep);
}, protected=TRUE) # extractRegions()


setMethodS3("extractChromosomes", "RawGenomicSignals", function(this, chromosomes=NULL, ...) {
  # Argument 'chromosomes':
  if (!is.null(chromosomes)) {
    chromosomes <- Arguments$getChromosomes(chromosomes);
    chromosomes <- unique(chromosomes);
    chromosomes <- sort(chromosomes);
  }

  # Nothing todo?
  if (is.null(chromosomes)) {
    return(this);
  }

  keep <- is.element(this$chromosome, chromosomes);
  keep <- which(keep);

  extractSubset(this, keep);
}, protected=TRUE) # extractChromosomes()


setMethodS3("extractChromosome", "RawGenomicSignals", function(this, chromosome, ...) {
  # Argument 'chromosome':
  chromosome <- Arguments$getChromosome(chromosome);

  extractChromosomes(this, chromosomes=chromosome, ...);
}, protected=TRUE) # extractChromosome()



setMethodS3("extractRawGenomicSignals", "default", abstract=TRUE, protected=TRUE);





setMethodS3("getPositions", "RawGenomicSignals", function(this, ...) {
  # This is a single-chromosome method. Assert that is the case.
  assertOneChromosome(this);

  x <- this$x;
  if (is.null(x)) {
    x <- seq_len(nbrOfLoci(this));
  }
  x;
})


setMethodS3("getChromosome", "RawGenomicSignals", function(this, ...) {
  # This is a single-chromosome method. Assert that is the case.
  assertOneChromosome(this);

  chr <- this$chromosome;
  if (is.null(chr)) {
    chr <- as.integer(NA);
  }
  chr <- as.integer(chr);
  chr;
})


setMethodS3("getSignals", "RawGenomicSignals", function(this, field=NULL, ...) {
  # Argument 'field':
  if (is.null(field)) {
    field <- getSignalColumnName(this);
  } else {
    field <- Arguments$getCharacter(field);
  }

  name <- field;

  # Sanity check
  if (!hasColumn(this, name, ...)) {
    throw(sprintf("Cannot get signals. No such column: %s not in (%s)", name, hpaste(getColumnNames(this, ...))));
  }

  this[[name]];
})


setMethodS3("setSignals", "RawGenomicSignals", function(this, values, ...) {
  name <- getSignalColumnName(this, ...);
  # Sanity check
  if (!hasColumn(this, name, ...)) {
    throw(sprintf("Cannot get signals. No such column: %s not in (%s)", name, hpaste(getColumnNames(this, ...))));
  }
  this[[name]] <- values;
  invisible(this);
})


setMethodS3("setWeights", "RawGenomicSignals", function(this, weights, ...) {
  # Argument 'weights':
  n <- nbrOfLoci(this);
  weights <- Arguments$getNumerics(weights, length=c(n,n), range=c(0,Inf));
  this$w <- weights;
  invisible(this);
})

setMethodS3("getWeights", "RawGenomicSignals", function(this, ...) {
  this$w;
})

setMethodS3("hasWeights", "RawGenomicSignals", function(this, ...) {
  (!is.null(getWeights(this, ...)));
})


setMethodS3("as.data.frame", "RawGenomicSignals", function(x, ..., sort=FALSE) {
  # To please R CMD check
  this <- x;

  # Sort along genome?
  if (sort) {
    this <- sort(this, ...);
  }

  NextMethod("as.data.frame");
})


setMethodS3("getDefaultLocusFields", "RawGenomicSignals", function(this, translate=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # "Default" fields
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- this;
  if (translate) {
    data <- as.data.frame(data, virtual=FALSE);
  }
  defFields <- colnames(data);

  defFields;
}, protected=TRUE) # getDefaultLocusFields()



setMethodS3("getLocusFields", "RawGenomicSignals", function(this, ...) {
  fields <- c(getDefaultLocusFields(this, ...), getVirtualColumnNames(this, ...));
  fields <- unique(fields);
  fields;
}, protected=TRUE) # getLocusFields()


setMethodS3("sort", "RawGenomicSignals", function(x, ...) {
  # To please R CMD check
  this <- x;

  # Order by (chromosome, x)
  o <- order(this$chromosome, getPositions(this));
  this[o,, drop=FALSE];
})


setMethodS3("getXY", "RawGenomicSignals", function(this,  ...) {
  # This is a single-chromosome method. Assert that's the case.
  assertOneChromosome(this);
  data <- getCXY(this, ...);
  data <- data[,c("x","y"), drop=FALSE];
  data;
}, protected=TRUE)


setMethodS3("getCXY", "RawGenomicSignals", function(this, ...) {
  data <- as.data.frame(this, ...);
  fields <- c("chromosome", "x", "y");
  data <- data[,fields,drop=FALSE];
  data;
}, protected=TRUE)






setMethodS3("kernelSmoothing", "RawGenomicSignals", function(this, xOut=NULL, ..., verbose=FALSE) {
  # This is a single-chromosome method. Assert that's the case.
  assertOneChromosome(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'xOut':
  if (!is.null(xOut)) {
    xOut <- Arguments$getNumerics(xOut);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Smoothing data set");
  x <- getPositions(this);
  y <- getSignals(this);

  if (is.null(xOut)) {
    xOut <- x;
  }

  verbose && cat(verbose, "xOut:");
  verbose && str(verbose, xOut);

  verbose && enter(verbose, "Kernel smoothing");
  verbose && cat(verbose, "Arguments:");
  args <- list(y=y, x=x, xOut=xOut, ...);
  verbose && str(verbose, args);
  ys <- kernelSmoothing(y=y, x=x, xOut=xOut, ...);
  verbose && str(verbose, ys);
  verbose && exit(verbose);


  verbose && enter(verbose, "Creating result object");
  # Allocate the correct size
  res <- newInstance(this, nrow=length(xOut));

  res$x <- xOut;
  res <- setSignals(res, ys);
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # kernelSmoothing()


setMethodS3("gaussianSmoothing", "RawGenomicSignals", function(this, sd=10e3, ...) {
  kernelSmoothing(this, kernel="gaussian", h=sd, ...);
})



setMethodS3("binnedSmoothing", "RawGenomicSignals", function(this, fields=NULL, ..., weights=getWeights(this), byCount=FALSE, verbose=FALSE) {
  # This is a single-chromosome method. Assert that's the case.
  assertOneChromosome(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  n <- nbrOfLoci(this);
  # Argument 'weights':
  if (!is.null(weights)) {
    weights <- Arguments$getNumerics(weights, length=c(n,n), range=c(0,Inf));
  }

  # Argument 'fields':
  if (is.null(fields)) {
    fields <- getSignalColumnName(this);
  } else {
    fields <- Arguments$getCharacters(fields);
    unknown <- fields[!hasColumns(this, fields)];
    if (length(unknown) > 0) {
      throw("No such field(s): ", hpaste(unknown));
    }
  }

  # Argument 'byCount':
  byCount <- Arguments$getLogical(byCount);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Smoothing data set");
  verbose && cat(verbose, "Smoothing fields: ", paste(fields, collapse=", "));
  dataY <- as.data.frame(this);
  dataY <- dataY[,fields,drop=FALSE];
  dataY <- as.matrix(dataY);
  verbose && str(verbose, dataY);

  verbose && cat(verbose, "Genomic positions:");
  x <- getPositions(this);
  verbose && str(verbose, x);
  xRange <- range(x, na.rm=TRUE);
  verbose && printf(verbose, "Range of positions: [%.0f,%.0f]\n",
                                           xRange[1], xRange[2]);


  wOut <- NULL;

  if (byCount) {
    verbose && enter(verbose, "Binned smoothing (by count)");
    # Smoothing y and x (and w).
    Y <- cbind(dataY, x=x, w=weights);
    verbose && summary(verbose, Y);
    xRank <- seq_len(nrow(Y));
    verbose && cat(verbose, "Positions (ranks):");
    verbose && str(verbose, xRank);
    verbose && cat(verbose, "Arguments:");
    args <- list(Y=Y, x=xRank, w=weights, ...);
    verbose && str(verbose, args);

    Ys <- colBinnedSmoothing(Y=Y, x=xRank, w=weights, ..., verbose=less(verbose, 10));
    verbose && str(verbose, Ys);
    verbose && summary(verbose, Ys);

    xOut <- attr(Ys, "xOut");
    verbose && str(verbose, xOut);

    # The smoothed x:s, which becomes the new target positions
    xOut <- Ys[,"x", drop=TRUE];
    verbose && str(verbose, xOut);

    # Smoothed weights
    if (!is.null(weights)) {
      wOut <- Ys[,"w", drop=TRUE];
      wOut[is.na(wOut)] <- 0;
      verbose && str(verbose, wOut);
    }

    # The smoothed y:s
    Ys <- Ys[,fields, drop=FALSE];
    verbose && str(verbose, Ys);

    # Not needed anymore
    xRank <- Y <- NULL;
    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "Binned smoothing (by position)");
    # Smoothing y (and w).
    Y <- cbind(dataY, w=weights);
    verbose && summary(verbose, Y);
    verbose && cat(verbose, "Arguments:");
    args <- list(Y=Y, x=x, w=weights, ...);
    verbose && str(verbose, args);

    Ys <- colBinnedSmoothing(Y=Y, x=x, w=weights, ..., verbose=less(verbose, 10));
    verbose && str(verbose, Ys);
    verbose && summary(verbose, Ys);

    xOut <- attr(Ys, "xOut");
    verbose && str(verbose, xOut);

    # Smoothed weights
    if (!is.null(weights)) {
      wOut <- Ys[,"w", drop=TRUE];
      wOut[is.na(wOut)] <- 0;
      verbose && str(verbose, wOut);
    }

    # The smoothed y:s
    Ys <- Ys[,fields, drop=FALSE];
    verbose && str(verbose, Ys);

    verbose && exit(verbose);
  } # if (byCount)

  verbose && enter(verbose, "Creating result object");
  # Allocate the correct size
  res <- newInstance(this, nrow=length(xOut));

  res$x <- xOut;
  res$w <- wOut;
  for (ff in fields) {
    res[[ff]] <- Ys[,ff,drop=TRUE];
  }
  verbose && print(verbose, summary(res));
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # binnedSmoothing()





setMethodS3("binnedSmoothingByField", "RawGenomicSignals", function(this, field, fields=NULL, from=xMin(this), to=xMax(this), by=NULL, length.out=NULL, byCount=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'field':
  field <- Arguments$getCharacter(field);
  if (!hasColumn(this, field)) {
    throw("Unknown field: ", field);
  }

  # Argument 'fields':
  if (is.null(fields)) {
    fields <- getSignalColumnName(this);
  } else {
    fields <- Arguments$getCharacters(fields);
    unknown <- fields[!hasColumns(this, fields)];
    if (length(unknown) > 0) {
      throw("No such field(s): ", hpaste(unknown));
    }
  }

  x <- getPositions(this);
  # Argument 'from' & 'to':
  if (is.null(from)) {
    from <- min(x, na.rm=TRUE);
  } else {
    from <- Arguments$getInteger(from);
  }
  if (is.null(to)) {
    to <- max(x, na.rm=TRUE);
  } else {
    to <- Arguments$getInteger(to, range=c(from, Inf));
  }

  # Arguments 'by' & 'length.out':
  if (is.null(by) & is.null(length.out)) {
    throw("Either argument 'by' or 'length.out' needs to be given.");
  }
  if (!is.null(by)) {
    by <- Arguments$getNumeric(by, range=c(0,to-from));
  }
  if (!is.null(length.out)) {
    length.out <- Arguments$getInteger(length.out, range=c(1,Inf));
  }

  # Argument 'byCount':
  byCount <- Arguments$getLogical(byCount);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Binning data set");
  verbose && cat(verbose, "Smoothing fields: ", paste(fields, collapse=", "));
  verbose && cat(verbose, "Stratifying by field: ", field);

  data <- as.data.frame(this);
  byValues <- data[[field]];
  setOfValues <- sort(unique(byValues), na.last=TRUE);
  verbose && cat(verbose, "Set of stratification values: ", hpaste(setOfValues));

  # Sanity check
  if (length(setOfValues) > 0.1*length(byValues)) {
    throw("Detected ridicolously large number of possible values in field '%s': (%s) [%d]", field, hpaste(setOfValues), length(setOfValues));
  }

  verbose && cat(verbose, "by: ", by);
  verbose && cat(verbose, "length.out: ", length.out);
  verbose && cat(verbose, "byCount: ", byCount);

  verbose && enter(verbose, "Find target positions");

  if (byCount) {
    verbose && enter(verbose, "By count");
    res <- sort(this);
    resOut <- binnedSmoothing(res, fields="x",
                              by=by, length.out=length.out, byCount=TRUE,
                              verbose=less(verbose, 5));
    xOut <- resOut$x;
    # Not needed anymore
    resOut <- NULL;
    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "By position");
    # Target 'x':
    if (!is.null(by)) {
      xOut <- seq(from=from, to=to, by=by);
    } else {
      xOut <- seq(from=from, to=to, length.out=length.out);
    }
    verbose && exit(verbose);
  } # if (byCount)

  verbose && cat(verbose, "xOut:");
  verbose && str(verbose, xOut);
  # Sanity check
  xOut <- Arguments$getNumerics(xOut);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate result set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocating result set");
  # Allocate results of the correct size
  res <- newInstance(this, nrow=length(xOut));

  # Target 'x' and 'y':
  res$x <- xOut;

  naValue <- as.double(NA);
  Ys <- matrix(naValue, nrow=length(xOut), ncol=length(fields));
  colnames(Ys) <- fields;
  for (ff in fields) {
    res[[ff]] <- Ys[,ff,drop=TRUE];
  }
  verbose && print(verbose, res);
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying output states (may drop short regions)");
  dataOut <- as.data.frame(res);
  # Sanity check
  stopifnot(is.element(field, colnames(dataOut)));

  byValuesOut <- dataOut[[field]];
  verbose && cat(verbose, "byValuesOut:");
  verbose && str(verbose, byValuesOut);
  setOfValuesOut <- sort(unique(byValuesOut), na.last=TRUE);
  verbose && printf(verbose, "Unique output '%s' values:\n", field);
  verbose && str(verbose, setOfValuesOut);
  verbose && exit(verbose);

  verbose && enter(verbose, "Setting up source signals");
  if (byCount) {
    # Adding ordering along genome
    gs <- sort(this);
    gs$xOrder <- seq_len(nbrOfLoci(gs));
  } else {
    gs <- this;
  }
  verbose && print(verbose, gs);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Binning (target) stratified by field
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (ss in seq_along(setOfValuesOut)) {
    byValue <- setOfValuesOut[ss];
    verbose && enter(verbose, sprintf("Value #%d (%s == '%s') of %d",
                                        ss, field, byValue, length(setOfValuesOut)));

    verbose && cat(verbose, "By value: ", byValue);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identify target loci with this byValue
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Extracting subset of (target) loci for this value");
    idxsOut <- which(is.element(byValuesOut, byValue));
    resSS <- extractSubset(res, idxsOut, verbose=less(verbose,50));
    verbose && print(verbose, resSS);
    xOutSS <- getPositions(resSS);
    verbose && str(verbose, xOutSS);
    verbose && exit(verbose);
    # Nothing to do? [Should actually never happen!]
    nbrOfBins <- nbrOfLoci(resSS);
    # Not needed anymore
    resSS <- NULL;
    if (nbrOfBins == 0) {
      verbose && cat(verbose, "No bins. Skipping byValue.");
      verbose && exit(verbose);
      next;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identifying source loci with this byValue
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Extracting subset of (source) loci with this signal value");
    byValues <- gs[[field]];
    idxsSS <- which(is.element(byValues, byValue));
    gsSS <- extractSubset(gs, idxsSS);
    verbose && print(verbose, gsSS);
    # Sanity check
    if (!all(is.element(gsSS[[field]], byValue))) {
       values <- gsSS[[field]];
       print(table(values, useNA="always"));
       print(byValue);
    }
    stopifnot(all(is.element(gsSS[[field]], byValue)));
    verbose && exit(verbose);

    # Nothing to do?
    if (nbrOfLoci(gsSS) == 0) {
      verbose && cat(verbose, "No extracted loci. Skipping byValue.");
      verbose && exit(verbose);
      next;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Extending bins to equal count sizes?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (byCount) {
      verbose && enter(verbose, "Extending bins to have equal number of loci");
      verbose && cat(verbose, "xOrder:");
      verbose && str(verbose, gsSS$xOrder);
      nbrOfLociToAdd <- (gsSS$xOrder[1] %% by) - 1;
      verbose && cat(verbose, "Number of loci to add: ", nbrOfLociToAdd);
      gsSS$xOrder <- NULL;

      if (nbrOfLociToAdd > 0) {
        verbose && enter(verbose, sprintf("Appending %d extra loci", nbrOfLociToAdd));

        # Allocate the correct size
        nbrOfLoci2 <- nbrOfLoci(gsSS) + nbrOfLociToAdd;
        gsSS2 <- newInstance(gsSS, nrow=nbrOfLoci2);

        fieldsT <- getColumnNames(gsSS, virtual=FALSE, translate=TRUE);
        verbose && cat(verbose, "Fields:");
        verbose && print(verbose, fieldsT);

        for (ff in seq_along(fieldsT)) {
          field <- fieldsT[ff];
          verbose && enter(verbose, sprintf("Field #%d ('%s') of %d", ff, field, length(fieldsT)));
          values <- gsSS[[field]];

          # Sanity check
          stopifnot(!is.null(values));

          if (is.element(field, "w")) {
            naValue <- 0;
          } else if (is.element(field, "x")) {
            naValue <- 0;
          } else {
            naValue <- NA;
          }
          storage.mode(naValue) <- storage.mode(values);
          naValues <- rep(naValue, times=nbrOfLociToAdd);
          values <- c(naValues, values);
          gsSS2[[field]] <- values;

          verbose && exit(verbose);
        } # for (ff ...)
        gsSS <- gsSS2;
        # Not needed anymore
        gsSS2 <- NULL;
        verbose && exit(verbose);
      } # if (nbrOfLociToAdd > 0)
      verbose && exit(verbose);
    } # if (byCount)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Bin loci of this byValue towards target loci (of the same byValue)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Binned smoothing of temporary object");
    if (byCount) {
      xOutSSRank <- seq(from=(1+by)/2, by=by, length.out=nbrOfBins);
      verbose && cat(verbose, "Arguments:");
      args <- list(gsSS, xOut=xOutSSRank, by=by, byCount=byCount, ...);
      verbose && str(verbose, args);
      resSS <- binnedSmoothing(gsSS, fields=fields,
                               xOut=xOutSSRank, by=by, byCount=byCount,
                               ..., verbose=less(verbose,5));
      resSS$x <- xOutSS;
    } else {
      verbose && cat(verbose, "Arguments:");
      args <- list(gsSS, xOut=xOutSS, by=by, byCount=byCount, ...);
      verbose && str(verbose, args);
      resSS <- binnedSmoothing(gsSS, fields=fields,
                               xOut=xOutSS, by=by, byCount=byCount,
                               ..., verbose=less(verbose,5));
    } # if (byCount)
    # Not needed anymore
    gsSS <- args <- NULL;
    verbose && print(verbose, resSS);
    verbose && exit(verbose);

    # Sanity check
    stopifnot(length(idxsOut) == nbrOfLoci(resSS));

    YsSS <- as.data.frame(resSS)[,fields, drop=FALSE];
    YsSS <- as.matrix(YsSS);
    Ys[idxsOut,fields] <- YsSS;
    # Not needed anymore
    resSS <- YsSS <- NULL;

    verbose && exit(verbose);
  } # for (ss ...)

  # Update 'res':
  for (ff in fields) {
    res[[ff]] <- Ys[,ff,drop=TRUE];
  }
  verbose && print(verbose, res);
  verbose && print(verbose, summary(res));

  verbose && exit(verbose);

  res;
}, protected=TRUE) # binnedSmoothingByField()




###########################################################################/**
# @set "class=RawGenomicSignals"
# @RdocMethod estimateStandardDeviation
#
# @title "Estimates the standard deviation of the raw Ys"
#
# \description{
#  @get "title" robustly or non-robustly using either a "direct" estimator
#  or a first-order difference estimator.
# }
#
# @synopsis
#
# \arguments{
#   \item{field}{A @character specifying the field to estimate.}
#   \item{method}{If \code{"diff"}, the estimate is based on the first-order
#     contigous differences of raw Ys. If \code{"direct"}, it is based
#     directly on the raw Ys.}
#   \item{estimator}{If \code{"mad"}, the robust @see "stats::mad" estimator
#     is used.  If \code{"sd"}, the @see "stats::sd" estimator is used.}
#   \item{na.rm}{If @TRUE, missing values are excluded first.}
#   \item{weights}{Locus specific weights.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a non-negative @numeric value.
# }
#
# @author
#
# \seealso{
#   @see "base::diff", @see "stats::sd", and @see "stats::mad".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("estimateStandardDeviation", "RawGenomicSignals", function(this, field=NULL, method=c("diff", "direct"), estimator=c("mad", "sd"), na.rm=TRUE, weights=getWeights(this), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'method':
  method <- match.arg(method);

  # Argument 'field':
  if (is.null(field)) {
    field <- getSignalColumnName(this);
  } else {
    field <- Arguments$getCharacter(field);
  }

  # Argument 'estimator':
  estimator <- match.arg(estimator);

  n <- nbrOfLoci(this);

  # Nothing todo?
  if (n <= 1) {
    return(as.double(NA));
  }

  # Argument 'weights':
  if (!is.null(weights)) {
    weights <- Arguments$getNumerics(weights, range=c(0,Inf), length=c(n,n));
  }

  # Get the estimator function
  if (!is.null(weights)) {
    estimator <- sprintf("weighted %s", estimator);
    estimator <- toCamelCase(estimator);
  }
  estimatorFcn <- get(estimator, mode="function");


  if (method == "diff") {
    # Sort along genome
    rgs <- sort(this);
    y <- rgs[[field]];

    # Insert NAs ("dividers") between chromosomes?
    if (nbrOfChromosomes(this) > 1) {
      chrs <- rgs$chromosome;
      dchrs <- diff(chrs);
      ats <- which(is.finite(dchrs) & dchrs > 0);
      y <- insert(y, ats=ats);  # R.utils::insert()
      if (!is.null(weights)) {
        weights <- insert(weights, ats=ats);
      }
      # Not needed anymore
      chrs <- dchrs <- ats <- NULL;
    }
    # Not needed anymore
    rgs <- NULL;

    y <- diff(y);

    # Weighted estimator?
    if (!is.null(weights)) {
      # Calculate weights per pair
      weights <- (weights[1:(n-1)]+weights[2:n])/2;
      sigma <- estimatorFcn(y, w=weights, na.rm=na.rm)/sqrt(2);
    } else {
      sigma <- estimatorFcn(y, na.rm=na.rm)/sqrt(2);
    }
  } else if (method == "direct") {
    y <- this[[field]];
    if (!is.null(weights)) {
      sigma <- estimatorFcn(y, w=weights, na.rm=na.rm);
    } else {
      sigma <- estimatorFcn(y, na.rm=na.rm);
    }
  }

  sigma;
})



setMethodS3("xSeq", "RawGenomicSignals", function(this, from=1, to=xMax(this), by=100e3, ...) {
  # This is a single-chromosome method. Assert that's the case.
  assertOneChromosome(this);

  seq(from=from, to=to, by=by);
})

setMethodS3("xRange", "RawGenomicSignals", function(this, na.rm=TRUE, ...) {
  # This is a single-chromosome method. Assert that's the case.
  assertOneChromosome(this);

  x <- getPositions(this);
  range(x, na.rm=na.rm);
})

setMethodS3("xMin", "RawGenomicSignals", function(this, ...) {
  xRange(this, ...)[1];
})

setMethodS3("xMax", "RawGenomicSignals", function(this, ...) {
  xRange(this, ...)[2];
})

setMethodS3("yRange", "RawGenomicSignals", function(this, na.rm=TRUE, ...) {
  y <- getSignals(this, ...);
  range(y, na.rm=na.rm);
})

setMethodS3("yMin", "RawGenomicSignals", function(this, ...) {
  yRange(this, ...)[1];
})

setMethodS3("yMax", "RawGenomicSignals", function(this, ...) {
  yRange(this, ...)[2];
})


setMethodS3("signalRange", "RawGenomicSignals", function(this, na.rm=TRUE, ...) {
  y <- getSignals(this, ...);
  range(y, na.rm=na.rm);
})

setMethodS3("setSigma", "RawGenomicSignals", function(this, sigma, ...) {
  sigma <- Arguments$getNumeric(sigma, range=c(0,Inf), disallow=NULL);
  this <- setBasicField(this, ".sigma", sigma);
  invisible(this);
})

setMethodS3("getSigma", "RawGenomicSignals", function(this, ..., force=FALSE) {
  sigma <- getBasicField(this, ".sigma");
  if (is.null(sigma)) {
    sigma <- estimateStandardDeviation(this, ...);
    this <- setSigma(this, sigma);
  }
  sigma;
})



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# TRICKS DURING Object -> BasicObject -> data.frame transition
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Adding clone() and clearCache() for RawGenomicSignals so that its
# methods work regardless of RawGenomicSignals extending BasicObject
# or Object (the original implementation).
setMethodS3("clone", "RawGenomicSignals", function(this, ...) {
  this;
}, protected=TRUE)

setMethodS3("clearCache", "RawGenomicSignals", function(this, ...) {
  this;
}, protected=TRUE)

setMethodS3("getBasicField", "RawGenomicSignals", function(this, key, ...) {
  attr(this, key);
}, protected=TRUE)

setMethodS3("setBasicField", "RawGenomicSignals", function(this, key, value, ...) {
  attr(this, key) <- value;
  invisible(this);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2012-03-14
# o Added argument 'field' to getSignals() and estimateStandardDeviation().
# o Added argument 'fields' to binnedSmoothing() for RawGenomicSignals.
#   Also to binnedSmoothingByField().
# o BUG FIX: binnedSmoothingByField(..., byCount=TRUE) would bin using
#   incorrect bin sets.
# 2012-03-13
# o CLEANUP: Dropped all code for backward compatibility for
#   RawGenericSignals inheriting from Object.
# 2012-03-12
# o Added binnedSmoothingByField() for RawGenomicSignals.
# o Added subset() for RawGenomicSignals.
# o Added (get|set)VirtualField().
# 2012-03-10
# o Now getDefaultLocusFields() utilizes as.data.frame().
# o Added argument 'sort' to as.data.frame() for RawGenomicSignals.
# o GENERALIZATION: Now nbrOfLoci() for RawGenomicSignals infers the
#   number of loci from the number of rows in as.data.frame(), instead
#   of assuming there exists a column names 'y'.
# 2012-03-02
# o Added newInstance() for RawGenomicSignals.
# o BUG FIX: Forgot to update kernelSmoothing() for data.frame:s.
# 2012-03-01
# o Now RawGenomicSignals inherits from data.frame().
# o Added (get|set)BasicField() for the Object -> BasicObject
#   -> data.frame transition.
# o Now all setNnn() methods returns invisible(this).
# o Preparing to support multiple-chromosome RawGenomicSignals;
#   - Added getChromosomes() and nbrOfChromosomes().
#   - Added single-chromosome assertions for methods assume that.
#   - Turn single-chromosome methods into multi-chromosome ones;
#     as.character(), (get|set)LocusFields(), sort(),
#     estimateStandardDeviation(), extractRegion(), extractRegions().
#   - Added getCXY().
#   - Added extractChromosome() and extractChromosomes().
#   - Added argument 'chromosome' to extractRegion().
#   - Added argument 'chromosomes' to extractRegions().
# o CLEANUP: Now extractRegion() and extractRegions() use extractSubset().
# 2012-02-04
# o GENERALIZATION: Now binnedSmoothing() of RawGenomicSignals default to
#   generate the same target bins as binnedSmoothing() of a numeric vector.
#   Before the bins had to be specified explicitly by the caller.
# 2011-12-15
# o ROBUSTNESS: Now binnedSmoothing(..., xOut) for RawGenomicSignals
#   guarantees to return length(xOut) loci.
# 2010-07-19
# o Now extractRegion() for RawGenomicSignals also accepts a
#   CopyNumberRegions object for argument 'regions'.
# o Added extractRegions() for RawGenomicSignals.
# 2009-10-10
# o Added setName().
# 2009-09-07
# o Added yRange(), yMin() and yMax() to RawGenomicSignals.
# 2009-07-03
# o BUG FIX: binnedSmoothing() added non existing locus field 'w'.
# 2009-06-30
# o Now binnedSmoothing() of RawGenomicSignals drops locus fields that were
#   not binned.  Ideally all locus fields (including custom ones) should be
#   binned, but we leave that for a future implementation.
# 2009-06-13
# o Now RawGenomicSignals(y=rgs) sets all locus fields in 'rgs' if it is
#   a RawGenomicSignals object.
# 2009-05-16
# o Now all methods of RawCopyNumbers() coerce numerics only if necessary,
#   i.e. it keeps integers if integers, otherwise to doubles.  This is a
#   general design of aroma.* that saves some memory.
# 2009-05-13
# o Now as.character() also reports mean distance between loci.
# o Added extractRegion() to RawGenomicSignals.
# o Now estimateStandardDeviation() takes a weighted estimate by default,
#   if weights are available.
# 2009-05-12
# o BUG FIX: extractSubset() of RawGenomicSignals did not recognize all
#   locus fields.
# o Now getSigma() returns the current standard deviation estimate, and
#   if not available, then it estimate it using estimateStandardDeviation().
# o Now binnedSmoothing() of RawGenomicSignals uses weighted estimates
#   (by default) if weights exists.
# 2009-05-10
# o Added argument 'w=NULL' to the constructor.
# o Added getWeights(), setWeights(), and hasWeights() to RawGenomicSignals.
# 2009-05-07
# o Added (get|set)(X|Y)Scale() to RawGenomicSignals.
# o Added setLocusFields().
# o Renamed getLociFields() to getLocusFields().
# o BUG FIX: lines() of RawGenomicSignals did not recognize x/yScale.
# 2009-04-06
# o BUG FIX: binnedSmoothing(..., byCount=TRUE) of RawGenomicSignals would
#   give error "[...] object "ys" not found".
# 2009-02-19
# o Renamed from RawCopyNumbers RawGenomicSignals.
# o Added argument 'byCount' to binnedSmoothing() of RawGenomicSignals.
# 2009-02-17
# o Now RawGenomicSignals() also takes another RawCopyNumbers object as
#   input.
# 2009-02-16
# o Added optional constructor argument 'name'.
# 2009-02-07
# o Added Rdoc comments and example.
# 2008-05-21
# o Added field 'chromosome' (single value).
# 2008-05-17
# o Added abstract default extractCopyNumberRegions().
# o Moved to aroma.core.
# 2008-03-31
# o Put recently added sd() and mad() into estimateStandardDeviation().
# 2008-03-10
# o Added standard deviation estimator sd() and mad() which my default
#   uses a first-order difference variance estimator.
# 2007-08-22
# o Created.  Need a generic container for holding copy number data and
#   to plot them nicely.
############################################################################
