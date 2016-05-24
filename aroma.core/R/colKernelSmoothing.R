###########################################################################/**
# @set "class=matrix"
# @RdocMethod colKernelSmoothing
# @alias colKernelSmoothing
# @alias kernelSmoothing
# @alias kernelSmoothing.numeric
#
# @title "Kernel smoothing of a matrix column by column"
#
# \description{
#  @get "title".
# }
# 
# @synopsis
#
# \arguments{
#   \item{Y}{A @numeric JxI @matrix (or a @vector of length J.)}
#   \item{x}{A (optional) @numeric @vector specifying the positions of
#     the J entries. The default is to assume uniformly distributed 
#     positions.}
#   \item{w}{A optional @numeric @vector of prior weights for each of 
#     the J entries.}
#   \item{xOut}{A @numeric @vector specifying K target positions where
#      the kernel is applied.}
#   \item{kernel}{A @character string or a @function specifying the
#      kernel used.}
#   \item{h}{A single positive @numeric specifying the bandwidth of 
#      the kernel.}
#   \item{censorH}{A single positive @numeric specifying the where to 
#      truncate the kernel. If @Inf, no truncation is done.}
#   \item{na.rm}{If @TRUE, missing values are excluded, otherwise not.}
#   \item{robust}{If @TRUE, robust estimators are used, otherwise not.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @numeric KxI @matrix (or a @vector of length K).
# }
#
# %% \details{
# %%   M_i' = w*M = w*(T-R) = w*T - w*R = T_i' - R'
# %%  
# %%   Before smoothing, the reference R_i == median(T_i). 
# %%   Keep this property for R' too.
# %%  
# %%   R' = median(T_i')
# %%   T_i' = M_i' - R'
# %%  
# %%   => w*T = w*M + w*R = M' + w*R
# %% }
#
# @examples "../incl/colKernelSmoothing.Rex"
#
# @author
#
# \seealso{
#   @seemethod "colBinnedSmoothing".
# }
#
# @keyword array
# @keyword iteration
# @keyword robust
# @keyword univar 
#*/###########################################################################
setMethodS3("colKernelSmoothing", "matrix", function(Y, x=seq_len(nrow(Y)), w=NULL, xOut=x, kernel=c("gaussian", "uniform"), h, censorH=3, na.rm=TRUE, robust=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'Y'
  n <- nrow(Y);
  k <- ncol(Y);
  
  # Argument 'x'
  if (length(x) != n) {
    throw("Argument 'x' has different number of values than rows in 'Y': ", 
                                                     length(x), " != ", n);
  }

  # Argument 'w'
  if (is.null(w)) {
  } else {
    if (length(w) != n) {
      throw("Argument 'w' has different number of values that rows in 'Y': ", 
                                                       length(w), " != ", n);
    }
  }

  # Argument 'xOut'
  if (is.null(xOut)) {
    xOut <- x;
  } else {
    xOut <- Arguments$getNumerics(xOut);
  }
  nOut <- length(xOut);

  # Argument 'kernel'
  if (is.character(kernel)) {
    kernel <- match.arg(kernel);
    if (kernel == "gaussian") {
      # wKernel <- dnorm(xDiff[keep], mean=0, sd=sd);
      kernel <- function(x, h) {
        dnorm(x, mean=0, sd=h);
      }
    } else if (kernel == "uniform") {
      kernel <- function(x, h) {
        w <- rep(0, times=length(x));
        w[-h/2 <= x & x < h/2] <- 1;
        w;
      }
    }
  }
  if (is.function(kernel)) {
  } else {
    throw("Argument 'kernel' must be either a string or a function: ", 
                                                              mode(kernel));
  }

  # Arguments 'h':
  h <- Arguments$getNumeric(h, range=c(0,Inf));

  # Arguments 'censorH':
  censorH <- Arguments$getNumeric(censorH, range=c(0,Inf));

  # Arguments 'robust':
  robust <- Arguments$getLogical(robust);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Setup (precalculations)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (robust) {
    # Exclude NAs?
    naRm <- ifelse(na.rm, TRUE, NA);
  }

  censorThreshold <- censorH * h;
  isCensored <- (censorThreshold < Inf);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Smoothing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Allocate vector of smoothed signals
  naValue <- as.double(NA);
  Ys <- matrix(naValue, nrow=nOut, ncol=k);
  colnames(Ys) <- colnames(Y);

  verbose && enter(verbose, "Estimating signals at given locations");

  verbose && cat(verbose, "Output locations:");
  verbose && str(verbose, xOut);

  # At each position in 'xOut', calculate the weighted average 
  # using the kernel.
  for (kk in seq_len(nOut)) {
    if (kk %% 100 == 0)
      verbose && cat(verbose, kk);

    # Weights centered around x[kk]
    xDiff <- (x-xOut[kk]);
    if (isCensored) {
      keep <- which(abs(xDiff) <= censorThreshold);
      # Nothing to do?
      if (length(keep) == 0) {
        next;
      }
      xDiff2 <- xDiff[keep];
      Y2 <- Y[keep,,drop=FALSE];
      w2 <- w[keep];
    } else {
      xDiff2 <- xDiff;
      Y2 <- Y;
      w2 <- w;
    }

    # Kernel weights...
    wKernel <- kernel(xDiff2, h=h);

    if (!is.null(w2)) {
      # ...with prior weights?
      wKernel <- w2 * wKernel;
    }

    if (robust) {
      value <- colWeightedMedians(Y2, w=wKernel, na.rm=naRm);
    } else {
      value <- colWeightedMeans(Y2, w=wKernel, na.rm=na.rm);
    }

    # Fix: Smoothing over a window with all missing values give zeros, not NA.
    idxs <- which(value == 0);
    if (length(idxs) > 0) {
      # Are these real zeros or missing values?
      Y2 <- Y2[idxs,,drop=FALSE];
      Y2 <- !is.na(Y2);
      idxsNA <- idxs[colSums(Y2) == 0];    
      value[idxsNA] <- NA;
    }

#    verbose && str(verbose, value);
    Ys[kk,] <- value;
  } # for (kk ...)

  verbose && exit(verbose);

  Ys;
}) # colKernelSmoothing()



setMethodS3("kernelSmoothing", "numeric", function(y, ...) {
  y <- colKernelSmoothing(as.matrix(y), ...);
  dim(y) <- NULL;
  y;
})



############################################################################
# HISTORY:
# 2012-03-14
# o Now colNnnSmoothing() returns a matrix with column name as
#   in argument 'Y'.
# 2009-05-16
# o Now colKernelSmoothing() uses Arguments$getNumerics(), not 
#   getDoubles(), where possible.  This will save memory in some cases.
# 2009-02-08
# o BUG FIX: Argument 'x' had the wrong default value.
# 2009-02-05
# o Renames matrix version to colKernelSmoothing().
# o Now making use of weighted matrix averages of matrixStats.
# o Added more validation of arguments.
# o Added Rdoc comments.
# 2008-11-06
# o Added argument 'robust' to kernelSmoothing().
# o Now kernelSmoothing() returns NA for windows with all missing values.
#   Before such windows would be returned as zero, cf. sum(c()) == 0.
# 2008-11-05
# o Created kernelSmoothing() from gaussianSmoothing().
# 2008-05-21
# o Added argument 'censorSd' to sensor the kernel at a given bandwidth.
# o Added argument 'xOut' to gaussianSmoothing().
# 2007-04-08
# o Added Gaussian smoothing for columns in a matrix.
# 2007-04-02
# o Created.
############################################################################
