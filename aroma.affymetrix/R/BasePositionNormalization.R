###########################################################################/**
# @RdocClass BasePositionNormalization
#
# @title "The BasePositionNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization method that corrects for systematic
#  effects in the probe intensities due to differences in positioning of
#  A, C, G, and T:s in the probe sequences.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "LinearModelProbeSequenceNormalization".}
#   \item{model}{A @character string specifying the model used to fit
#     the base-count effects.}
#   \item{df}{The degrees of freedom of the model.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB, MR"
#*/###########################################################################
setConstructorS3("BasePositionNormalization", function(..., model=c("smooth.spline"), df=5) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  # Argument 'df':
  df <- Arguments$getInteger(df, range=c(1,1e3));

  extend(LinearModelProbeSequenceNormalization(...), "BasePositionNormalization",
    .model = model,
    .df = df
  )
})


setMethodS3("getAsteriskTags", "BasePositionNormalization", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Add model tag?
  model <- this$.model;
  if (model != "smooth.spline") {
    tags <- c(tags, model);
  }

  # Add df tag?
  df <- this$.df;
  if (df != 5) {
    tags <- c(tags, sprintf("df=%d", df));
  }

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)



setMethodS3("getParameters", "BasePositionNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  params <- c(params, list(
    model = this$.model,
    df = this$.df
  ));

  params;
}, protected=TRUE)



setMethodS3("getDesignMatrix", "BasePositionNormalization", function(this, cells=NULL, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
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

  verbose && enter(verbose, "Getting algorithm parameters");
  params <- getParameters(this, expand=FALSE, verbose=less(verbose, 1));
  model <- params$model;
  df <- params$df;
  verbose && cat(verbose, "Model: ", model);
  verbose && cat(verbose, "Degrees of freedom: ", df);
  verbose && exit(verbose);

  # Locate AromaCellSequenceFile holding probe sequences
  acs <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check file cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(
    method="getDesignMatrix", class=class(this[1]),
    cells=cells,
    model=model, df=df,
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

  verbose && enter(verbose, "Reading probe sequences");
  seqs <- readSequenceMatrix(acs, cells=cells, what="raw",
                                                verbose=less(verbose, 5));
  # Not needed anymore
  acs <- NULL;
  verbose && cat(verbose, "Probe-sequence matrix:");
  verbose && str(verbose, seqs);
  verbose && exit(verbose);

  verbose && enter(verbose, "Building probe-position design matrix");
  verbose && cat(verbose, "Degrees of freedom: ", df);
  X <- getProbePositionEffectDesignMatrix(seqs, df=df,
                                               verbose=less(verbose, 5));
  # Not needed anymore
  seqs <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && cat(verbose, "Design matrix:");
  verbose && str(verbose, X);
  verbose && cat(verbose, "RAM: ", object.size(X));
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Cache results?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    saveCache(X, key=key, dirs=dirs);
  }

  verbose && exit(verbose);

  X;
}, private=TRUE)



setMethodS3("getSignalTransform", "BasePositionNormalization", function(this, ...) {
  params <- getParameters(this, expand=FALSE, ...);
  shift <- params$shift;
  # Not needed anymore
  params <- NULL;

  transform <- function(y, ...) {
    y <- y + shift;
    y <- log2(y);
    y;
  }

  transform;
}, protected=TRUE);



############################################################################
# HISTORY:
# 2008-11-29
# o Extracted the LinearModelProbeSequenceNormalization class from this
#   class.
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
