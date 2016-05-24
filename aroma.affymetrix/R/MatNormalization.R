###########################################################################/**
# @RdocClass MatNormalization
#
# @title "The MatNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization method that corrects for systematic
#  effects in the probe intensities due to differences in the number of
#  A, C, G, and T:s and the match scores according to MAT [1].
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "AbstractProbeSequenceNormalization".}
#   \item{unitsToFit}{The units from which the normalization curve should
#     be estimated.  If @NULL, all are considered.}
#   \item{model}{A @character string specifying the model used to fit
#     the base-count effects.}
#   \item{nbrOfBins}{The number of bins to use for the variance smoothing step.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Requirements}{
#   This class requires that an aroma probe sequence file and aroma
#   match scores file is available for the chip type.
# }
#
# @author "MR"
#
# \references{
#   [1] Johnson WE, Li W, Meyer CA, Gottardo R, Carroll JS, Brown M, Liu XS.
#     \emph{Model-based analysis of tiling-arrays for ChIP-chip}, PNAS, 2006.
# }
#*/###########################################################################
setConstructorS3("MatNormalization", function(..., unitsToFit=NULL, model=c("lm"), nbrOfBins=200) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  # Argument 'nbrOfBins':
  nbrOfBins <- Arguments$getInteger(nbrOfBins, range=c(1,Inf));

  args <- list(...);
  if (is.element("numChunks", names(args))) {
    throw("Argument 'numChunks' to MatNormalization() is deprecated.  Instead, use the 'ram' option in the aroma settings, cf. http://www.aroma-project.org/settings/");
  }
  if (is.element("numBins", names(args))) {
    throw("Argument 'numBins' is deprecated.  Instead, use argument 'nbrOfBins'.");
  }

  extend(AbstractProbeSequenceNormalization(..., unitsToFit=unitsToFit), "MatNormalization",
    .model = model,
    .scaleResiduals = TRUE,
    .nbrOfBins = as.integer(nbrOfBins)
  )
})


setMethodS3("getAromaCellMatchScoreFile", "MatNormalization", function(this, ..., force=FALSE) {
  apm <- this$.apm;

  if (force || is.null(apm)) {
    dataSet <- getInputDataSet(this);
    cdf <- getCdf(dataSet);
    chipType <- getChipType(cdf, fullname=FALSE);
    nbrOfCells <- nbrOfCells(cdf);
    apm <- AromaCellMatchScoreFile$byChipType(chipType, nbrOfCells=nbrOfCells, ...);
    this$.apm <- apm;
  }

  apm;
}, protected=TRUE)



setMethodS3("getAsteriskTags", "MatNormalization", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Add model tag?
  model <- this$.model;
  tags <- c(tags, model);

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)


setMethodS3("getParameters", "MatNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  params <- c(params, list(
    model = this$.model,
    nbrOfBins = this$.nbrOfBins
  ));

  params;
}, protected=TRUE)



setMethodS3("getDesignMatrix", "MatNormalization", function(this, cells=NULL, model=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  if (is.null(cells)) {
  } else {
    # Validated below...
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving design matrix");
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);

  verbose && enter(verbose, "Locating probe-sequence annotation data");
  # Locate AromaCellSequenceFile holding probe sequences
  aps <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Locating match scores annotation data");
  # Locate AromaCellMatchScoreFile holding match scores
  apm <- getAromaCellMatchScoreFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading sequence matrix");
  sm <- readSequenceMatrix(aps, cells=cells, verbose=verbose);
  verbose && exit(verbose);
  verbose && enter(verbose, "Reading match scores");
  ms <- readColumns(apm, rows=cells, verbose=verbose);
  verbose && exit(verbose);

  verbose && enter(verbose, "Constructing design matrix");
  nT <- rowSums(sm == "T");
  G <- (sm == "G")+0;
  A <- (sm == "A")+0;
  C <- (sm == "C")+0;
  designMatrix <- cbind(nT, A, C, G, rowSums(A)^2, rowSums(C)^2, rowSums(G)^2, nT^2, log(as.integer(ms[,1])));

  # Garbage collect
  # Not needed anymore
  nT <- G <- A <- C <- ms <- sm <- NULL;
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && str(verbose, designMatrix);
  verbose && cat(verbose, "object.size(designMatrix): ",
                                             object.size(designMatrix));
  verbose && exit(verbose);

  verbose && exit(verbose);

  designMatrix;
}, private=TRUE)




setMethodS3("fitOne", "MatNormalization", function(this, df, ram=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating probe-sequence annotation data");
  # Locate AromaCellSequenceFile holding probe sequences
  aps <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Locating match scores annotation data");
  # Locate AromaCellMatchScoreFile holding match scores
  apm <- getAromaCellMatchScoreFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading 'non-missing' cells to fit");
  cellsToFit <- which( !(isMissing(aps) | isMissing(apm)) );
  nbrOfCells <- length(cellsToFit);
  verbose && cat(verbose, "Cells to fit:");
  verbose && str(verbose, cellsToFit);
  verbose && exit(verbose);

  # this code adopted from Dave Fourniers 17/08/2007 post to r-help
  # mailing list entitled "[R] Linear models over large datasets"

  # Calculate the number of cells to process per chunk. It should:
  # (1) increase with the 'ram' option.
  # (2) there is only on array, so it should be independent of everything.
  cellsPerChunk <- ceiling(ram * 1e6 + 1L);
  cellsPerChunk <- Arguments$getInteger(cellsPerChunk, range=c(1,Inf));

  ## In aroma.affymetrix v1.9.4 and before, the number of cells processed
  ## per chunk increase with the total number of cells. /HB 2011-02-15
  ## cellsPerChunk <- ceiling(nbrOfCells/numChunks) + 1;

  nbrOfChunks <- ceiling(nbrOfCells / cellsPerChunk);
  verbose && cat(verbose, "Number cells per chunk: ", cellsPerChunk);

  verbose && enter(verbose, "Reading signals to fit");
  y <- extractMatrix(df, cells=cellsToFit, verbose=less(verbose, 10));
  verbose && exit(verbose);

  verbose && enter(verbose, "Log2 transforming signals");
  y <- log2(y);
  verbose && cat(verbose, "Target log2 probe signals:");
  verbose && str(verbose, y);
  verbose && exit(verbose);

  start <- 0L;
  xtx <- xty <- 0;

  chunk <- 1L;
  while (start < nbrOfCells) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", chunk, nbrOfChunks));

    verbose && enter(verbose, "Working on indices over range");
    from <- start + 1L;
    to <- min(start+cellsPerChunk, nbrOfCells);
    indSubset <- (from:to);
    rng <- c(from, to);
    rng <- rng / nbrOfCells;
    verbose && cat(verbose, sprintf("[%g,%g]", rng[1], rng[2]));
    verbose && exit(verbose);

    verbose && enter(verbose, "Reading design matrix for this iteration");
    X <- getDesignMatrix(this, cells=cellsToFit[indSubset], verbose=verbose);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating cross products");
    xtx <- xtx + crossprod(X);
    xty <- xty + crossprod(X, y[indSubset]);
    verbose && exit(verbose);

    start <- start + cellsPerChunk;

    chunk <- chunk + 1L;
    verbose && exit(verbose);
  } # while (...)

  verbose && enter(verbose, "Solving normal equations");
  #fit <- list(xtx=xtx, xty=xty) #,beta=solve(xtx, xty))
  verbose && exit(verbose);
  fit <- list(beta=solve(xtx, xty), scaleResiduals=this$.scaleResiduals);

  fit;
}, protected=TRUE)


setMethodS3("predictOne", "MatNormalization", function(this, fit, ram=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating probe-sequence annotation data");
  # Locate AromaCellSequenceFile holding probe sequences
  aps <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Locating match scores annotation data");
  # Locate AromaCellMatchScoreFile holding match scores
  apm <- getAromaCellMatchScoreFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Allocating mu vector");
  mu <- double(nbrOfCells(aps));
  verbose && str(verbose, mu);
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading 'non-missing' cells to predict");
  cellsToPredict <- which( !(isMissing(aps) | isMissing(apm)) );
  nbrOfCells <- length(cellsToPredict);
  verbose && cat(verbose, "Cells to predict:");
  verbose && str(verbose, cellsToPredict);
  verbose && exit(verbose);

  # Calculate the number of cells to process per chunk. It should:
  # (1) increase with the 'ram' option.
  # (2) there is only on array, so it should be independent of everything.
  cellsPerChunk <- ceiling(ram * 1e6 + 1L);
  cellsPerChunk <- Arguments$getInteger(cellsPerChunk, range=c(1,Inf));

  ## In aroma.affymetrix v1.9.4 and before, the number of cells processed
  ## per chunk increase with the total number of cells. /HB 2011-02-15
  ## cellsPerChunk <- ceiling(nbrOfCells/numChunks) + 1;

  nbrOfChunks <- ceiling(nbrOfCells / cellsPerChunk);
  verbose && cat(verbose, "Number cells per chunk: ", cellsPerChunk);

  start <- 0L;

  chunk <- 1L;
  while (start < nbrOfCells) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", chunk, nbrOfChunks));

    verbose && enter(verbose, "Working on indices over range");
    from <- start + 1L;
    to <- min(start+cellsPerChunk, nbrOfCells);
    indSubset <- (from:to);
    rng <- c(from, to);
    rng <- rng / nbrOfCells;
    verbose && cat(verbose, sprintf("[%g,%g]", rng[1], rng[2]));
    verbose && exit(verbose);

    verbose && enter(verbose, "Reading design matrix for this iteration");
    X <- getDesignMatrix(this, cells=cellsToPredict[indSubset], verbose=verbose)
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating mu");
    mu[ cellsToPredict[indSubset] ] <- X %*% fit$beta;
    verbose && exit(verbose);

    start <- start + cellsPerChunk;

    chunk <- chunk + 1L;
    verbose && exit(verbose);
  } # while (...)

  # Not needed anymore
  X <- indSubset <- NULL;
  gc <- gc();

  #nbrOfBins <- this$.nbrOfBins
  #q <- quantile(mu[cellsToPredict],prob=(0:nbrOfBins)/nbrOfBins)
  #cuts<-cut(mu[cellsToPredict],breaks=q,labels=1:(length(q)-1))  # define
  #ss<-split(data.frame(resid),cuts)
  #ssvar<-sapply(ss,var)
  #v<-ssvar[as.character(cuts)]
  #for(j in 1:length(b))
  #rr<-resid/sqrt(v)

  # Return results
  mu;
}, protected=TRUE)


###########################################################################/**
# @RdocMethod process
#
# @title "Normalizes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already normalized is re-normalized,
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "MatNormalization", function(this, ..., ram=NULL, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Normalizing data set for probe-sequence effects");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already normalized");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get (and create) the output path
  outputPath <- getPath(this);

  verbose && enter(verbose, "Locating probe-sequence annotation data");
  # Locate AromaCellSequenceFile holding probe sequences
  aps <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Locating match scores annotation data");
  # Locate AromaCellMatchScoreFile holding match scores
  apm <- getAromaCellMatchScoreFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading 'non-missing' cells to fit");
  cellsToFit <- which( !(isMissing(aps) | isMissing(apm)) );
  verbose && cat(verbose, "Cells to fit:");
  verbose && str(verbose, cellsToFit);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize all arrays simultaneously
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(ds);
  verbose && enter(verbose, "Normalizing ", nbrOfArrays, " arrays");
  verbose && enter(verbose, "Path: ", outputPath);

  df <- getOneFile(ds);
  nbrOfCells <- nbrOfCells(df);

  xtx <- 0;
  xtyList <- as.list(double(nbrOfArrays));

  nbrOfCells <- length(cellsToFit);

  # Calculate the number of cells to process per chunk. It should:
  # (1) increase with the 'ram' option.
  # (2) the arrays are processed sequentially, i.e. it is constant
  #     in number of arrays.
  cellsPerChunk <- ceiling(ram * 1e6 + 1L);
  cellsPerChunk <- Arguments$getInteger(cellsPerChunk, range=c(1,Inf));

  ## In aroma.affymetrix v1.9.4 and before, the number of cells processed
  ## per chunk increase with the total number of cells. /HB 2011-02-15
  ## cellsPerChunk <- ceiling(nbrOfCells/numChunks) + 1;

  nbrOfChunks <- ceiling(nbrOfCells / cellsPerChunk);
  verbose && cat(verbose, "Number cells per chunk: ", cellsPerChunk);

  idxs <- 1:nbrOfCells;
  head <- 1:cellsPerChunk;
  count <- 1L;
  while (length(idxs) > 0) {
    verbose && enter(verbose, sprintf("Fitting chunk #%d of %d", count, nbrOfChunks));
    if (length(idxs) < cellsPerChunk) {
      head <- 1:length(idxs);
    }
    cc <- idxs[head];

    verbose && cat(verbose, "Cells: ");
    verbose && str(verbose, cellsToFit[cc]);

    verbose && enter(verbose, "Reading design matrix");
    X <- getDesignMatrix(this, cells=cellsToFit[cc], verbose=verbose);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating cross product X'X");
    xtx <- xtx + crossprod(X);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating cross product X'y for each array");
    for (ii in seq_len(nbrOfArrays)) {
      df <- ds[[ii]];
      verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                          ii, getName(df), nbrOfArrays));

      y <- extractMatrix(df, cells=cellsToFit[cc], verbose=verbose);
      y <- log2(y);
      verbose && cat(verbose, "Target log2 probe signals:");
      verbose && str(verbose, y);

      xtyList[[ii]] <- xtyList[[ii]] + crossprod(X, y);

      # Not needed anymore
      y <- NULL;
      verbose && exit(verbose);
    } # for (ii ...)
    verbose && exit(verbose);

    # Not needed anymore
    X <- NULL;

    # Next chunk
    idxs <- idxs[-head];
    count <- count + 1L;

    verbose && exit(verbose);
  }  # while (...)

  verbose && enter(verbose, "Solving for each array");
  fits <- lapply(xtyList, FUN=function(xty) {
    list(beta=solve(xtx, xty));
  });
  verbose && exit(verbose);

  # Not needed anmore
  # Not needed anymore
  xtx <- xtyList <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save model fits
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (ii in seq_len(nbrOfArrays)) {
    df <- ds[[ii]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                              ii, getName(df), nbrOfArrays));

    fullname <- getFullName(df);

    # Store model fit
    verbose && enter(verbose, "Saving model fit");
    # Store fit and parameters (in case someone are interested in looking
    # at them later; no promises of backward compatibility though).
    filename <- sprintf("%s,fit.RData", fullname);
    fitPathname <- Arguments$getWritablePathname(filename,
                                                    path=outputPath, ...);
    saveObject(fits[[ii]], file=fitPathname);
    verbose && str(verbose, fits[[ii]], level=-50);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);
  } # for (ii ...)



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create output CEL files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating (temporary) CEL files")

  pathnamesT <- listenv()
  for (ii in seq_len(nbrOfArrays)) {
    df <- ds[[ii]]
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                              ii, getName(df), nbrOfArrays))

    fullname <- getFullName(df)
    filename <- sprintf("%s.CEL", fullname)
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...)

    # Write to a temporary file (allow rename of existing one if forced)
    isFile <- isFile(pathname)
    pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose)

    # Create CEL file to store results, if missing
    verbose && enter(verbose, "Creating CEL file for results, if missing")
    createFrom(df, filename=pathnameT, path=NULL, verbose=less(verbose))
    verbose && exit(verbose)

    # Record temporary filename and reuse below
    pathnamesT[[ii]] <- pathnameT

    verbose && exit(verbose);
  } ## for (ii ...)
  ## Assert proper temporary files
  pathnamesT <- unlist(pathnamesT)
  stopifnot(length(pathnamesT) == nbrOfArrays)

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Predict
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfCells <- length(cellsToFit);

  # Calculate the number of cells to process per chunk. It should:
  # (1) increase with the 'ram' option.
  # (2) the arrays are processed sequentially, i.e. it is constant
  #     in number of arrays.
  cellsPerChunk <- ceiling(ram * 1e6 + 1L);
  cellsPerChunk <- Arguments$getInteger(cellsPerChunk, range=c(1,Inf));

  ## In aroma.affymetrix v1.9.4 and before, the number of cells processed
  ## per chunk increase with the total number of cells. /HB 2011-02-15
  ## cellsPerChunk <- ceiling(nbrOfCells/numChunks) + 1;

  nbrOfChunks <- ceiling(nbrOfCells / cellsPerChunk);
  verbose && cat(verbose, "Number cells per chunk: ", cellsPerChunk);

  start <- 0L;
  xtx <- 0;
  mu <- vector("list", length=nbrOfArrays);

  count <- 1L;
  while (start < nbrOfCells) {
    verbose && enter(verbose, sprintf("Fitting chunk #%d of %d", count, nbrOfChunks));

    verbose && enter(verbose, "Working on indices over range");
    from <- start + 1L;
    to <- min(start+cellsPerChunk, nbrOfCells);
    indSubset <- (from:to);
    rng <- c(from, to);
    rng <- rng / nbrOfCells;
    verbose && cat(verbose, sprintf("[%g,%g]", rng[1], rng[2]));
    verbose && exit(verbose);

    verbose && enter(verbose, "Set of cells");
    cellsChunk <- cellsToFit[indSubset];
    verbose && str(verbose, cellsChunk);
    verbose && exit(verbose);

    verbose && enter(verbose, "Reading design matrix for this iteration");
    X <- getDesignMatrix(this, cells=cellsChunk, verbose=verbose);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating model fits");
    xtx <- xtx + crossprod(X);

    verbose && enter(verbose, "Processing ", nbrOfArrays, " arrays")
    for (ii in seq_len(nbrOfArrays)) {
      df <- ds[[ii]]
      verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                              ii, getName(df), nbrOfArrays))

      mu <- X %*% fits[[ii]]$beta
      mu <- as.double(mu)
      verbose && str(verbose, mu)

      verbose && enter(verbose, "Updating temporary CEL file")
      pathnameT <- pathnamesT[ii]
      verbose2 <- as.logical(verbose)
      .updateCel(pathnameT, indices=cellsChunk, intensities=2^mu, verbose=verbose2)
      verbose && exit(verbose)

      verbose && exit(verbose)
    } # for (ii ...)
    verbose && exit(verbose)

    # Not needed anymore
    mu <- NULL;

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);

    start <- start + cellsPerChunk;
    #start <- nbrOfCells + 1;

    count <- count + 1L;
    verbose && exit(verbose);
  } # while (...)




  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Scale residuals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Scaling residuals")

  nbrOfBins <- this$.nbrOfBins
  verbose && cat(verbose, "Number of bins: ", nbrOfBins)
  probs <- (0:nbrOfBins) / nbrOfBins
  verbose && cat(verbose, "Quantile probabilities:")
  verbose && str(verbose, probs)
  cutLabels <- seq_len(nbrOfBins)

  res <- listenv()

  for (ii in seq_len(nbrOfArrays)) {
    df <- ds[[ii]]
    verbose && enter(verbose, "Binning predicted values, calculating and scaling residuals")
    pathnameT <- pathnamesT[ii]

    res[[ii]] %<=% {
      ## Original probe intensities
      y <- extractMatrix(df, cells=cellsToFit, verbose=verbose)
      y <- log2(y)

      ## Prediction
      mu <- .readCel(pathnameT, indices=cellsToFit, readOutliers=FALSE, readHeader=FALSE, readMasked=FALSE, verbose=less(verbose,10))$intensities
      mu <- log2(mu)

      ## Residuals
      r <- y - mu

      ## Normalize residuals
      q <- quantile(mu, probs=probs)
      cuts <- cut(mu, breaks=q, labels=cutLabels)  # define
      ss <- split(r, cuts, drop=FALSE)
      ssvar <- sapply(ss, FUN=var)
      v <- ssvar[as.character(cuts)]
      r <- r / sqrt(v)
      r <- as.double(r)

      #return(list(y=y,mu=mu,r=r))

      verbose && enter(verbose, "Updating temporary CEL file")
      verbose2 <- as.logical(verbose)
      .updateCel(pathnameT, indices=cellsToFit, intensities=2^r, verbose=verbose2)
      verbose && exit(verbose)

      # Not needed anymore
      q <- ss <- ssvar <- v <- r <- y <- NULL

      gc <- gc()
      verbose && print(verbose, gc)

      # Rename temporary file
      pathname <- popTemporaryFile(pathnameT, verbose=verbose)

      ## Create checksum file
      dfZ <- getChecksumFile(pathname)

      pathname
    } ## %<=%

    verbose && exit(verbose)
  } # for (ii ...)

  ## Resolve futures
  res <- as.list(res)
  res <- NULL

  verbose && exit(verbose)

  outputDataSet <- getOutputDataSet(this, force=TRUE)

  verbose && exit(verbose)

  invisible(outputDataSet)
})



############################################################################
# HISTORY:
# 2011-02-15 [HB]
# o Now we use updateCel(..., verbose=as.logical(verbose)) instead of
#   always verbose=TRUE.
# o HARMONIZATION: Renamed argument 'numBins' of MatNormalization to
#   'nbrOfBins'.
# o HARMONIZATION: Now MatNormalization utilized the aroma setting 'ram',
#   which replaced argument 'numChunks' which is now deprecated.
# o CLEANUP: Tidied up code.
# 2009-05-23 [HB]
# o Updated some minor format mistakes in the verbose output.
# 2008-11-28 [HB]
# o Added protected getCrossProductXTX().  Still not used.
# o Modified the first loop over chunks that calculated cross products such
#   that it is constant in number of arrays.  The processing over chunks is
#   now also done as we do it elsewhere in the package.
# o fitOne() and predictOne() are never used?!?
# o Updated the Rdocs.
# o Renamed getAromaCellMatchscoreFile() to getAromaCellMatchScoreFile().
# 2008-10-29 [MR]
# o Created from BaseCountNormalization.R
############################################################################
