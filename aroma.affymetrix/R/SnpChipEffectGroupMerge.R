###########################################################################/**
# @RdocClass SnpChipEffectGroupMerge
#
# @title "The SnpChipEffectGroupMerge class"
#
# \description{
#  @classhierarchy
#
#  This class represents a method that merges SNP chip effects across groups
#  unit by unit.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#        @see "ChipEffectGroupMerge".}
#   \item{mergeStrands}{If @TRUE, group strands are merged.}
#   \item{mean}{A @character string specifying what type of averaging
#        should be applied.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# @keyword "internal"
#*/###########################################################################
setConstructorS3("SnpChipEffectGroupMerge", function(..., mergeStrands=FALSE, mean=c("arithmetic", "geometric")) {
  # Argument 'mean':
  mean <- match.arg(mean);

  extend(ChipEffectGroupMerge(...), "SnpChipEffectGroupMerge",
    mergeStrands = mergeStrands,
    .mean = mean
  )
})


setMethodS3("getMergeFunction", "SnpChipEffectGroupMerge", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mergeStrandsMatrix <- function(y, g=function(x) x, h=function(x) x, ...) {
    n <- nrow(y);
    if (n %in% c(2,4)) {
      yy <- y[1:2,,drop=FALSE];
      yy <- g(yy);
      yy <- colMeans(yy, na.rm=TRUE);
      yy <- h(yy);
      y[1,] <- yy;
      y[2,] <- 0;
      if (n == 4) {
        yy <- y[3:4,,drop=FALSE];
        yy <- log2(yy);
        yy <- colMeans(yy, na.rm=TRUE);
        yy <- 2^yy;
        y[3,] <- yy;
        y[4,] <- 0;
      }
    }
    y;
  }


  # Get the merge function
  mean <- this$.mean;
  fcn <- NULL;
  if (this$mergeStrands) {
    if (mean == "geometric") {
      fcn <- function(y) {
        mergeStrandsMatrix(y, g=log2, h=function(x) 2^x);
      }
    } else if (mean == "arithmetic") {
      fcn <- function(y) {
        mergeStrandsMatrix(y);
      }
    }
  }

  fcn;
})


setMethodS3("getAsteriskTags", "SnpChipEffectGroupMerge", function(this, ...) {
  tags <- NULL;
  if (this$mergeStrands)
    tags <- c(tags, "F+R");

  tags;
}, protected=TRUE)


setMethodS3("getParameters", "SnpChipEffectGroupMerge", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  # Get parameters of this class
  params2 <- list(
    mergeStrands = this$mergeStrands
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2007-06-11
# o BUG FIX: getMergeFunction() of SnpChipEffectGroupMerge used non-existing
#   'x...' instead of 'x'.
# 2007-02-20
# o Created.
############################################################################
