###########################################################################/**
# @RdocClass LinearModelProbeSequenceNormalization
#
# @title "The LinearModelProbeSequenceNormalization class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a normalization method that corrects
#  for systematic effects in the probe intensities due to probe-sequence
#  dependent effects that can be modelled using a linear model.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "AbstractProbeSequenceNormalization".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Requirements}{
#   This class requires that an aroma probe sequence file is available
#   for the chip type.
# }
#
# \section{Memory usage}{
#  The model fitting methods of this class are bounded in memory.
#  This is done by first building up the normal equations incrementally
#  in chunks of cells.  The generation of normal equations is otherwise
#  the step that consumes the most memory.
#  When the normal equations are available, the @see "base::solve"
#  method is used to solve the equations.  Note that this algorithm is
#  still exact.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("LinearModelProbeSequenceNormalization", function(...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  extend(AbstractProbeSequenceNormalization(...), "LinearModelProbeSequenceNormalization"
  )
})



setMethodS3("getDesignMatrix", "LinearModelProbeSequenceNormalization", abstract=TRUE, protected=TRUE)



setMethodS3("getNormalEquations", "LinearModelProbeSequenceNormalization", function(this, df, cells=NULL, ram=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'df':
  df <- Arguments$getInstanceOf(df, "AffymetrixCelFile");

  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Getting annotation data files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving cell sequence annotation data file");
  acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 20));
  verbose && exit(verbose);

  # Expand 'cells'?
  if (is.null(cells)) {
    cells <- seq_len(nbrOfCells(acs));
  }

  verbose && enter(verbose, "Retrieving signal transform");
  transform <- getSignalTransform(this);
  verbose && str(verbose, transform);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identifying cells with known sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying subset of cells with known probe sequences");
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);
  n0 <- length(cells);

  isMissing <- isMissing(acs, verbose=less(verbose, 10))[cells];
  cells <- cells[!isMissing];
  # Not needed anymore
  isMissing <- NULL;
  n1 <- length(cells);
  verbose && printf(verbose, "Removed %d (%.2f%%) missing sequences out of %d\n", n0-n1, 100*(n0-n1)/n0, n0);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Loading signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading signals for these cells");

  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);
  y <- extractMatrix(df, cells=cells, drop=TRUE, verbose=less(verbose, 10));

  if (!is.null(transform)) {
    verbose && enter(verbose, "Transforming signals");
    verbose && cat(verbose, "Signals before transformation:");
    verbose && str(verbose, y);
    y <- transform(y);
    verbose && exit(verbose);
  }

  verbose && cat(verbose, "Signals to be fitted:");
  verbose && str(verbose, y);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Incrementally build up the normal equations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating number of cells per chunk");
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);
  nbrOfCells <- length(cells);
  verbose && cat(verbose, "Number of cells: ", nbrOfCells);

  verbose && cat(verbose, "RAM scale factor: ", ram);

  cellsPerChunk <- ram*1e6;
  verbose && cat(verbose, "Cells per chunk: ", cellsPerChunk);

  nbrOfChunks <- ceiling(nbrOfCells / cellsPerChunk);
  verbose && cat(verbose, "Number of chunks: ", nbrOfChunks);
  verbose && exit(verbose);

  xtx <- 0;
  xty <- 0;

  idxs <- 1:nbrOfCells;
  head <- 1:cellsPerChunk;
  count <- 1;
  while (length(idxs) > 0) {
    verbose && enter(verbose, "Processing chunk #", count, " of ", nbrOfChunks);
    if (length(idxs) < cellsPerChunk) {
      head <- 1:length(idxs);
    }
    cc <- idxs[head];

    verbose && cat(verbose, "Cells: ");
    verbose && str(verbose, cells[cc]);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Getting design matrix for subset
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Getting design matrix");
    # cache=TRUE => Cache to file.
    res <- getDesignMatrix(this, cells=cells[cc], cache=TRUE,
                                                 verbose=less(verbose, 5));
    X <- res$X;

    verbose && cat(verbose, "Design matrix:");
    verbose && str(verbose, X);

    B <- res$B;
    verbose && cat(verbose, "Basis vectors:");
    verbose && str(verbose, B);
    map <- res$map;

    factors <- res$factors;
    verbose && cat(verbose, "Factors:");
    verbose && str(verbose, factors);

    # Not needed anymore
    res <- NULL;

    gc <- gc();
    verbose && print(verbose, gc);

    yCC <- y[cc];

    # Sanity check
    stopifnot(nrow(X) == length(cells[cc]));
    stopifnot(nrow(X) == length(yCC));
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Keep only data points with finite signals
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Excluding non-finite data points");
    keep <- which(is.finite(yCC));
    yCC <- yCC[keep];
    X <- X[keep,,drop=FALSE];
    # Not needed anymore
    keep <- NULL;
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculating cross products X'X and X'y
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Calculating cross product X'X");
    xtxChunk <- crossprod(X);
    verbose && str(verbose, xtxChunk);
    xtx <- xtx + xtxChunk;
    # Not needed anymore
    xtxChunk <- NULL;
    verbose && exit(verbose);


    verbose && enter(verbose, "Calculating cross product X'y");
    xtyChunk <- crossprod(X, yCC);
    verbose && str(verbose, xtyChunk);
    xty <- xty + xtyChunk;
    # Not needed anymore
    X <- yCC <- xtyChunk <- NULL;
    verbose && exit(verbose);

    # Clean up

    # Next chunk
    idxs <- idxs[-head];
    count <- count + 1;

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # while (length(idxs) > 0);

  verbose && cat(verbose, "Normal equations:");
  verbose && str(verbose, xtx);
  verbose && str(verbose, xty);

  res <- list(xtx=xtx, xty=xty, n0=n0, n1=n1, cells=cells, map=map, B=B, factors=factors, y=y);
  # Not needed anymore
  xtx <- xty <- cells <- NULL;

  res;
}, protected=TRUE)




setMethodS3("fitOne", "LinearModelProbeSequenceNormalization", function(this, df, params=NULL, ram=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'df':
  df <- Arguments$getInstanceOf(df, "AffymetrixCelFile");

  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting normalization function for one array");
  verbose && cat(verbose, "Full name: ", getFullName(df));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Getting algorithm parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting algorithm parameters");
  if (is.null(params)) {
    params <- getParameters(this, expand=TRUE, verbose=less(verbose, 5));
    gc <- gc();
    verbose && print(verbose, gc);
  } else {
    verbose && cat(verbose, "Passed internally");
  }

  cells <- params$cellsToFit;
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);

  # Model parameters (only for display; retrieve elsewhere)
  verbose && cat(verbose, "Model: ", params$model);
  verbose && cat(verbose, "Degrees of freedom: ", params$df);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fitting model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Exact fitting of model by incrementally building the normal equations (X'X = X'y) and then solve it");


  verbose && enter(verbose, "Get normal equations X'X = X'y");
  ne <- getNormalEquations(this, df=df, cells=cells, ram=ram, verbose=verbose);
  verbose && cat(verbose, "Normal equations:");
  verbose && str(verbose, ne);

  map <- ne$map;
  B <- ne$B;
  factors <- ne$factors;

  xtx <- ne$xtx;
  xty <- ne$xty;
  # Not needed anymore
  ne <- NULL;
  verbose && exit(verbose);


  verbose && enter(verbose, "Solving normal equations");
  coefs <- solve(xtx, xty);
  coefs <- as.vector(coefs);
  verbose && cat(verbose, "Coeffients:")
  verbose && print(verbose, coefs);
  # Not needed anymore
  xtx <- xty <- NULL;
  verbose && exit(verbose);


  verbose && enter(verbose, "Restructuring results");
  params <- list();
  intercept <- TRUE;
  if (intercept) {
    params$intercept <- coefs[1];
    coefs <- coefs[-1];
  }
  df <- length(coefs)/length(factors);
  verbose && cat(verbose, "Degrees of freedom: ", df);
  idxs <- seq_len(df);
  for (kk in seq_along(factors)) {
    key <- names(factors)[kk];
    if (is.null(key)) {
      key <- sprintf("factor%02d", kk);
    }
    params[[key]] <- coefs[idxs];
    coefs <- coefs[-idxs];
  } # for (kk ...)

  fit <- list(params=params, map=map, B=B, algorithm="solve");
  class(fit) <- "ProbePositionEffects";
  verbose && exit(verbose);


  verbose && str(verbose, fit);

  verbose && exit(verbose);

  fit;
}, protected=TRUE)



setMethodS3("predictOne", "LinearModelProbeSequenceNormalization", function(this, fit, params=NULL, seqs=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Predicting model for one array");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Getting algorithm parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting algorithm parameters");
  if (is.null(params)) {
    params <- getParameters(this, expand=TRUE, verbose=less(verbose, 5));
    gc <- gc();
    verbose && print(verbose, gc);
  } else {
    verbose && cat(verbose, "Passed internally");
  }

  cells <- params$cellsToUpdate;
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);

  # Other model parameters
  model <- params$model;
  verbose && cat(verbose, "Model: ", model);
  verbose && exit(verbose);

  # Not needed anymore
  params <- NULL;

  nbrOfCells <- length(cells);


  verbose && enter(verbose, "Retrieving probe sequences");
  if (is.null(seqs)) {
    # Locate AromaCellSequenceFile holding probe sequences
    acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
    seqs <- readSequenceMatrix(acs, cells=cells, what="raw",
                                         verbose=less(verbose, 5));
    # Not needed anymore
    acs <- cells <- NULL;
    gc <- gc();
    verbose && print(verbose, gc);
  } else {
    verbose && cat(verbose, "Passed internally");
  }
  verbose && cat(verbose, "Probe-sequence matrix:");
  verbose && str(verbose, seqs);
  verbose && exit(verbose);

  verbose && enter(verbose, "Predicting mean (transformed) probe signals");
  mu <- predict(fit, seqs=seqs, verbose=less(verbose, 5));
  # Not needed anymore
  seqs <- NULL;
  verbose && str(verbose, "mu:");
  verbose && str(verbose, mu);
  verbose && summary(verbose, mu);
  verbose && exit(verbose);

  # Sanity check
  if (length(mu) != nbrOfCells) {
    throw("Internal error. Number of estimated means does not match the number of cells to be updated: ", length(mu), " != ", nbrOfCells);
  }

  # Garbage collection
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  mu;
}, protected=TRUE)


setMethodS3("getSignalTransform", "LinearModelProbeSequenceNormalization", function(this, ...) {
  NULL;
}, protected=TRUE)

############################################################################
# HISTORY:
# 2010-02-16
# o MEMORY OPTIMIZATION: Now getNormalEquations() no longer return (stray)
#   design matrix 'X' of the last processed chunk.
# 2008-12-03
# o SPEED UP: Fixed broken memoization of getNormalEquations() and fitOne().
# o SPEED UP: Now predictOne() accepts optional 'seqs'.
# o Added abstract getDesignMatrix() method.
# 2008-11-29
# o Extracted from BasePositionNormalization.R.
# o Now fitOne() takes an argument 'ram' which is passed from process().
# o The predictOne() method is looping over probe positions, which is
#   fairly memory efficient.  For this reason, we leave it as it.  We
#   can now fit a GenomeWideSNP_6 array with approx 1GB of RAM (instead
#   of 5-6GB before)!
# o Now getNormalEquations() is done in chunks. For GenomeWideSNP_6 we can
#   now generate normal equations with approx 500MB of RAM.
# o Added first step toward supporting fitting the linear model in
#   bounded memory.  This is done by setting up the normal equations and
#   using solve(xtx, xty) to estimate the parameters.  TEST: modelMethod
#   "lm.fit" and "solve" created all.equal() == TRUE output.
#   NEXT: Build up the NE incrementally.  Already without this, the memory
#   usage went down dramatically.  For a Mapping50K_Hind240 fit, the peak
#   memory usage went down from 1000MB to 380MB.  However, it is still not
#   possible to fit a GenomeWideSNP_6 on Windows Vista 32-bit.
# o Dropped the bootstrapping framework.
# 2008-07-29
# o Added support for specifying the degrees of freedom ('df') of the model.
# 2008-07-28
# o Updated to work with newer ProbeLevelTransform3.
# 2008-07-21
# o BENCHMARKING: For a GenomeWideSNP_6,Full, the BPN peaks at 5.9GB RAM.
#   This happens while fitting the model.  Prediction peaks at 3.2GB RAM.
# o Now getDesignMatrix() caches results to file.
# o Created from BaseCountNormalization.R.
############################################################################
