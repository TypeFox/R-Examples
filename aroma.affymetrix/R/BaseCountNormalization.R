###########################################################################/**
# @RdocClass BaseCountNormalization
#
# @title "The BaseCountNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization method that corrects for systematic
#  effects in the probe intensities due to differences in the number of
#  A, C, G, and T:s in the probe sequences.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "AbstractProbeSequenceNormalization".}
#   \item{model}{A @character string specifying the model used to fit
#     the base-count effects.}
#   \item{bootstrap}{If @TRUE, the model fitting is done by bootstrap in
#     order to save memory.}
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
# @author "HB"
#*/###########################################################################
setConstructorS3("BaseCountNormalization", function(..., model=c("robustSmoothSpline", "lm"), bootstrap=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  # Argument 'bootstrap':
  bootstrap <- Arguments$getLogical(bootstrap);

  if (bootstrap && model == "lm") {
    throw("Bootstrapping for the 'lm' model is not implemented.");
  }


  extend(AbstractProbeSequenceNormalization(...), "BaseCountNormalization",
    .model = model,
    .bootstrap=bootstrap,
    .chunkSize=as.integer(2.5e6),
    .maxIter=as.integer(50),
    .acc=0.005
  )
})


setMethodS3("getAsteriskTags", "BaseCountNormalization", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Add model tag?
  model <- this$.model;
  if (model != "robustSmoothSpline") {
    tags <- c(tags, model);
  }

  # Add bootstrap tag?
  if (this$.bootstrap) {
    bootstrapTag <- "B";
    tags <- c(tags, bootstrapTag);
  }

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)



setMethodS3("getParameters", "BaseCountNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  params <- c(params, list(
    model = this$.model,
    bootstrap = this$.bootstrap,
    chunkSize = this$.chunkSize,
    maxIter = this$.maxIter,
    acc = this$.acc
  ));

  params;
}, protected=TRUE)



setMethodS3("getDesignMatrix", "BaseCountNormalization", function(this, cells=NULL, model=NULL, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
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

  verbose && cat(verbose, "Model: ", model);

  # Retrieve nucleotide counts as raw values whenever
  # possible to save memory.
  if (model == "lm") {
    mode <- "integer";
    oneValue <- as.integer(1);
  } else if (model == "robustSmoothSpline") {
    mode <- "raw";
    oneValue <- as.raw(1);
  }


  verbose && enter(verbose, "Locating probe-sequence annotation data");
  # Locate AromaCellSequenceFile holding probe sequences
  acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check file cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(
    method="getDesignMatrix", class=class(this[1]),
    cells=cells,
    model=model,
    acs=list(fullname=getFullName(acs))
  );

  dirs <- c("aroma.affymetrix", getChipType(acs));
  if (!force) {
    X <- loadCache(key=key, dirs=dirs);
    if (!is.null(X)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(X);
    }
  }


  verbose && enter(verbose, "Count nucleotide bases for *all* cells");
  verbose && cat(verbose, "Chip type: ", getChipType(acs));
  designMatrix <- countBases(acs, mode=mode, verbose=less(verbose, 5));
  # Not needed anymore
  acs <- NULL;
  verbose && cat(verbose, "Nucleotide base counts:");
  verbose && str(verbose, designMatrix);
  verbose && cat(verbose, "object.size(designMatrix): ",
                                            object.size(designMatrix));
  verbose && exit(verbose);

  if (!is.null(cells)) {
    verbose && enter(verbose, "Extracing sequences of interest");
    verbose && cat(verbose, "Cells:");
    verbose && str(verbose, cells);
    nbrOfCells <- nrow(designMatrix);
    cells <- Arguments$getIndices(cells, max=nbrOfCells);
    designMatrix <- designMatrix[cells,,drop=FALSE];
    # Not needed anymore
    cells <- NULL;

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && cat(verbose, "object.size(designMatrix): ",
                                            object.size(designMatrix));
    verbose && exit(verbose);
  }

  designMatrix[,1] <- oneValue;
  verbose && cat(verbose, "Design matrix:");
  verbose && str(verbose, designMatrix);
  verbose && cat(verbose, "object.size(designMatrix): ",
                                            object.size(designMatrix));

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Cache results?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    saveCache(X, key=key, dirs=dirs);
  }

  verbose && exit(verbose);

  designMatrix;
}, private=TRUE)



setMethodS3("fitOne", "BaseCountNormalization", function(this, df, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'y' and 'X' must not contain NAs.
  fitBaseCounts <- function(y, X, subset=NULL, model=c("robustSmoothSpline", "lm"), ...) {
    # Argument 'y':

    # Argument 'X':
    if (nrow(X) != length(y)) {
      throw("The number of rows in design matrix 'X' does not match the number of observations in 'y': ", nrow(X), " != ", length(y));
    }

    if (!is.null(subset)) {
      y <- y[subset];
      X <- X[subset,,drop=FALSE];
      gc <- gc();
    }

    # Argument 'model':
    model <- match.arg(model);

    if (model == "lm") {
      requireNamespace("stats") || throw("Package not loaded: stats")
      lm.fit <- stats::lm.fit

      fitFcn <- function(X, y, ...) {
        fit <- lm.fit(x=X, y=y, ...);
        # Remove redundant parameters
        for (ff in c("residuals", "effects", "fitted.values", "qr")) {
          fit[[ff]] <- NULL;
        }
        fit;
      }
    } else if (model == "robustSmoothSpline") {
      fitFcn <- function(X, y, ...) {
        fits <- list();
        for (cc in 1:ncol(X)) {
          # Fit effect of term #cc
          if (cc == 1) {
            mu <- median(y);
            fit <- list(mu=mu);
          } else {
            # Note: 'X' may be a "raw" matrix (to save memory)
            Xcc <- as.double(X[,cc]);
            fit <- .robustSmoothSpline(x=Xcc, y=y, ...);
##            # Remove redundant parameters (although really small here)
##            for (ff in c("x", "y", "w", "yin", "lev")) {
##              fit[[ff]] <- NULL;
##            }
            mu <- predict(fit, x=Xcc)$y;
            # Not needed anymore
            Xcc <- NULL;
          }

          # Remove the effect of term #cc
          y <- y - mu;
          # Not needed anymore
          mu <- NULL;

          fits[[cc]] <- fit;
          # Not needed anymore
          fit <- NULL;
        }
        fits;
      }
    }

    fit <- fitFcn(X, y);

    fit;
  } # fitBaseCounts()


  fitSubset <- function(df, cells=NULL, ..., verbose) {
    verbose && enter(verbose, "Reading signals to fit");
    verbose && cat(verbose, "Cells:");
    verbose && str(verbose, cells);
    y <- extractMatrix(df, cells=cells, drop=TRUE, verbose=less(verbose, 10));
    verbose && exit(verbose);

    if (shift != 0) {
      verbose && enter(verbose, "Shifting signals");
      verbose && cat(verbose, "Shift: ", shift);
      y <- y + shift;
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Log2 transforming signals");
    y <- log2(y);
    verbose && cat(verbose, "Target log2 probe signals:");
    verbose && str(verbose, y);
    verbose && exit(verbose);

    verbose && enter(verbose, "Identify finite data points");
    n <- length(y);
    # Fit only finite subset
    keep <- which(is.finite(y));
    y <- y[keep];
    cells <- cells[keep];
    # Not needed anymore
    keep <- NULL;
    n2 <- length(y);
    verbose && printf(verbose, "Removed %d (%.4f%%) out of %d non-finite data points: ", n-n2, 100*(n-n2)/n, n);
    gc <- gc();
    verbose && exit(verbose);

    verbose && enter(verbose, "Fitting base-count model");
    X <- getDesignMatrix(this, cells=cells, model=model, verbose=less(verbose, 5));
    # Not needed anymore
    cells <- NULL;
    verbose && cat(verbose, "Design matrix:");
    verbose && str(verbose, X);
    gc <- gc();
    verbose && print(verbose, gc);
    fit <- fitBaseCounts(y, X=X, model=model, verbose=less(verbose, 5));
    # Not needed anymore
    y <- X <- NULL;
    verbose && print(verbose, fit);
    verbose && exit(verbose);

    fit;
  } # fitSubset()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'df':
  df <- Arguments$getInstanceOf(df, "AffymetrixCelFile");

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
  params <- getParameters(this, expand=TRUE, verbose=less(verbose, 5));
  gc <- gc();
  verbose && print(verbose, gc);

  units <- params$unitsToFit;
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);
  cells <- params$cellsToFit;
  stratifyBy <- params$typesToFit;
  verbose && cat(verbose, "stratifyBy: ", stratifyBy);
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);
  shift <- params$shift;
  verbose && cat(verbose, "Shift: ", shift);

  # Other model parameters
  model <- params$model;
  verbose && cat(verbose, "Model: ", model);
  bootstrap <- params$bootstrap;
  chunkSize <- params$chunkSize;
  maxIter <- params$maxIter;
  acc <- params$acc;
  verbose && exit(verbose);

  nbrOfCells <- length(cells);

  # Is bootstrapping necessary?
  if (chunkSize >= nbrOfCells) {
    verbose && cat(verbose, "Bootstrapping not really needed.");
    bootstrap <- FALSE;
  }

  # Bootstrap?
  if (bootstrap) {
    verbose && enter(verbose, "Fitting model using bootstrap");

    verbose && cat(verbose, "Max number of iterations: ", maxIter);
    verbose && cat(verbose, "Accuracy threshold: ", acc);

    bb <- 0;
    fit <- NULL;
    deltaMax <- Inf;
    while (deltaMax > acc && bb < maxIter) {
      bb <- bb + 1;
      verbose && enter(verbose, sprintf("Bootstrap iteration #%d", as.integer(bb)));

      verbose && enter(verbose, "Fitting model to subset of data");
      verbose && cat(verbose, "Chunk size: ", chunkSize);
      subset <- sample(1:nbrOfCells, size=chunkSize);
      cellsChunk <- cells[subset];
      # Not needed anymore
      subset <- NULL;
      cellsChunk <- sort(cellsChunk);
      verbose && cat(verbose, "Cells:");
      verbose && str(verbose, cellsChunk);
      fitB <- fitSubset(df, cells=cellsChunk, verbose=verbose);
      # Not needed anymore
      cellsChunk <- NULL;
      verbose && exit(verbose);

      verbose && enter(verbose, "Updating estimates");
      if (is.null(fit)) {
        fit <- fitB;
      } else {
        # Update parameters for each subfit
        deltaMax <- 0;
        for (kk in seq_along(fit)) {
          fitKK <- fit[[kk]];

          fitBKK <- fitB[[kk]];
          if (kk == 1) {
            # Average mean estimates
            fitKK$mu <- mean(c(fitKK$mu, fitBKK$mu), na.rm=TRUE);
          } else {
            # Average spline estimates
            # Common x:s
            x <- sort(unique(c(fitKK$x, fitBKK$x)));
            # Predict y:s based on bootstrap pool and new fit
            y <- fitKK$y;
            yB <- predict(fitBKK, x=x)$y;
            # Weight them together
            yB <- ((bb-1)*y + yB)/bb;
            delta <- mean(abs(y - yB), na.rm=TRUE);
            deltaMax <- max(c(deltaMax, delta), na.rm=TRUE);
            fitKK <- list(x=x, y=yB, delta=delta);
            # Not needed anymore
            x <- y <- yB <- delta <- NULL;
          }

          fit[[kk]] <- fitKK;

          # Not needed anymore
          fitKK <- fitBKK <- NULL;
        } # for (kk ...)
      }
      # Not needed anymore
      fitB <- NULL;
      verbose && exit(verbose);

      verbose && printf(verbose, "deltaMax < threshold: %.6f < %.6f\n", deltaMax, acc);
      verbose && printf(verbose, "iteration < maxIter: %d < %d\n", as.integer(bb), as.integer(maxIter));

      verbose && exit(verbose);
    } # while()
    converged <- (deltaMax <= acc);

    verbose && enter(verbose, "Creating final fit");
    for (kk in seq(from=2, to=length(fit))) {
      fitKK <- fit[[kk]];
      fitKK <- .robustSmoothSpline(x=fitKK$x, y=fitKK$y);
      fit[[kk]] <- fitKK;
    }
    verbose && exit(verbose);

    fit$bootstrap <- list(iter=as.integer(bb), maxIter=maxIter, converged=converged);

    verbose && exit(verbose);
  } else {
    fit <- fitSubset(df, cells=cells, verbose=verbose);
  } # if (bootstrap)

  verbose && exit(verbose);

  fit;
}, protected=TRUE)


setMethodS3("predictOne", "BaseCountNormalization", function(this, fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'y' and 'X' must not contain NAs.
  predictBaseCounts <- function(fit, X, model=c("robustSmoothSpline", "lm"), ...) {
    # Argument 'model':
    model <- match.arg(model);

    if (model == "lm") {
      predictFcn <- function(fit, X, ...) {
        coefs <- coefficients(fit);
        coefs <- as.matrix(coefs);
        yPred <- X %*% coefs;
        yPred;
      } # predictFcn()
    } else if (model == "robustSmoothSpline") {
      predictFcn <- function(fit, X, ...) {
        fits <- fit;
        yPred <- double(nrow(X));
        for (cc in 1:ncol(X)) {
          fit <- fits[[cc]];
          if (cc == 1) {
            mu <- fit$mu;
          } else {
            mu <- rep(as.double(NA), times=nrow(X));
            if (mode(X) == "raw") {
              # Note: 'X' may be a "raw" matrix (to save memory)
              idxs <- which(X[,cc] != as.raw(255));
            } else {
              idxs <- which(is.finite(X[,cc]));
            }
            x <- as.double(X[idxs,cc]);
            mu[idxs] <- predict(fit, x=x)$y;
            # Not needed anymore
            idxs <- x <- NULL;
          }
          str(mu);
          yPred <- yPred + mu;
          # Not needed anymore
          fit <- mu <- NULL;
        } # for (cc ...)
        yPred;
      } # predictFcn()
    }

    yPred <- predictFcn(fit, X);

    yPred;
  } # predictBaseCounts()


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
  params <- getParameters(this, expand=TRUE, verbose=less(verbose, 5));
  gc <- gc();
  verbose && print(verbose, gc);

  units <- params$unitsToUpdate;
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);
  cells <- params$cellsToUpdate;
  stratifyBy <- params$typesToUpdate;
  verbose && cat(verbose, "stratifyBy: ", stratifyBy);
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);
  shift <- params$shift;
  verbose && cat(verbose, "Shift: ", shift);

  # Other model parameters
  model <- params$model;
  verbose && cat(verbose, "Model: ", model);
  bootstrap <- params$bootstrap;
  chunkSize <- params$chunkSize;
  maxIter <- params$maxIter;
  acc <- params$acc;
  verbose && exit(verbose);

  nbrOfCells <- length(cells);

  verbose && enter(verbose, "Predicting mean log2 probe signals");
  X <- getDesignMatrix(this, cells=cells, model=model, verbose=less(verbose, 5));
  # Not needed anymore
  cells <- NULL;
  verbose && cat(verbose, "Design matrix:");
  verbose && str(verbose, X);

  mu <- predictBaseCounts(fit, X=X, model=model);
  # Not needed anymore
  fit <- X <- NULL;
  verbose && str(verbose, "mu:");
  verbose && str(verbose, mu);
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



############################################################################
# HISTORY:
# 2008-12-03
# o Updated getDesignMatrix() to cache results.
# 2008-07-28
# o Updated to work with newer ProbeLevelTransform3.
# 2008-07-21
# o Now BaseCountNormalization inherits from the abstract class
#   AbstractProbeSequenceNormalization, which already has process().
# 2008-07-20
# o Added an option to fit the model using bootstrap techiques.  It is not
#   obvious how to merge two smooth spline estimates for which the degrees
#   of freedom do not have to match.  Now it is done by predicting outcomes
#   at common covariates and using weighted average to combine the outcomes.
#   The weight of the pooled estimate increase by each bootstrap iteration.
# o MEMORY OPTIMIZATION: The below changes and using "raw" design matrices
#   have saved a lot of memory.  Before this, it peaked at 2.2-2.4GB of RAM
#   to normalize GWS6 data whereas now it peaks at approx 1.4GB.
# o Now fitOne() and predictOne() can handle "raw" design matrices.
# o Now using more efficient which() instead of which().
# o Added protected getDesignMatrix(), fitOne(), and predictOne().
#   This makes the process() code cleaner and we can save more memory by
#   utilizing file caching.
# o Now class inherits from ProbeLevelTransform2.
# 2008-07-19
# o Removed countBases() for this class.
# o Added 'shift' defaulting to zero.
# o Now asterisk tag is sensitive to the model and adds the model as a tag
#   for models other than robustSmoothSpline.
# 2008-07-17
# o BUG FIX: getAromaCellSequenceFile() would search using the full name of
#   the chip type, e.g. GenomeWideSNP_6,Full.
# 2008-07-16
# o Added support for fitting and updating subsets of cells and types of
#   probes according to the CDF.
# 2008-07-11
# o Updated countBases() to use new AromaCellSequenceFile class.
# 2008-06-22
# o Created from AllelicCrosstalkCalibration.R.
############################################################################
