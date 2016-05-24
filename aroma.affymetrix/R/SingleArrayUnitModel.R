###########################################################################/**
# @RdocClass SingleArrayUnitModel
#
# @title "The SingleArrayUnitModel class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a unit model that fits one model per unit
#  based on signals from a single arrays.
#  The nature of a single-array unit model is that each array can be fitted
#  independently of the others.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "UnitModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("SingleArrayUnitModel", function(...) {
  extend(UnitModel(...), "SingleArrayUnitModel")
}, abstract=TRUE)




###########################################################################/**
# @RdocMethod getFitUnitGroupFunction
#
# @title "Static method to get the low-level function that fits the PLM"
#
# \description{
#  @get "title".
#  Any subclass model must provide this method, which should return
#  a @function that accepts a @numeric @vector of length K, where K
#  is the number of probes.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @function.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFitUnitGroupFunction", "SingleArrayUnitModel", abstract=TRUE, static=TRUE, private=TRUE);


setMethodS3("getFitUnitFunction", "SingleArrayUnitModel", function(this, ...) {
  # Get the fit function for a single set of intensities
  fitfcn <- getFitUnitGroupFunction(this, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the one for all blocks in a unit
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (this$probeModel == "pm-mm") {
    fitUnit <- function(unit, ...) {
      lapply(unit, FUN=function(group) {
        y <- .subset2(group, 1); # Get intensities
        y <- y[1,] - y[2,];  # PM-MM
        fitfcn(y);
      })
    }
  } else if (this$probeModel == "min1(pm-mm)") {
    fitUnit <- function(unit, ...) {
      lapply(unit, FUN=function(group) {
        y <- .subset2(group, 1); # Get intensities
        y <- y[1,] - y[2,];  # PM-MM
        y[y < 1] <- 1;       # min1(PM-MM)=min(PM-MM,1)
        fitfcn(y);
      })
    }
  } else if (this$probeModel == "pm+mm") {
    fitUnit <- function(unit, ...) {
      lapply(unit, FUN=function(group) {
        y <- .subset2(group, 1); # Get intensities
        y <- y[1,] + y[2,];  # PM+MM
        fitfcn(y);
      })
    }
  } else {
    fitUnit <- function(unit, ...) {
      lapply(unit, FUN=function(group) {
        if (length(group) > 0) {
          y <- .subset2(group, 1); # Get intensities
        } else {
          y <- NULL;
        }
        fitfcn(y);
      })
    }
  }

  fitUnit;
}, private=TRUE)





###########################################################################/**
# @RdocMethod readUnits
#
# @title "Reads data unit by unit"
#
# \description{
#  @get "title" for all or a subset of units (probeset)
#  specially structured for this PLM.
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be read. If @NULL, all units are read.}
#   \item{...}{Arguments passed to \code{getCellIndices()} of the
#     @see "AffymetrixCdfFile" class (if \code{cdf} was not specified),
#     but also to the \code{readUnits()} method of the
#     @see "AffymetrixCelFile" class.}
# }
#
# \value{
#  Returns the @list structure that \code{readUnits()} of
#  @see "AffymetrixCelFile" returns.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("readUnits", "SingleArrayUnitModel", function(this, array, units=NULL, ..., verbose=FALSE) {
  ds <- getDataSet(this);
  nbrOfArrays <- length(ds);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'array':
  array <- Arguments$getIndex(array, max=nbrOfArrays);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading probe intensities from array #", array);
  df <- ds[[array]];
  verbose && print(verbose, df);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the CDF cell indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying CDF cell indices");
  cdfUnits <- getCellIndices(this, units=units, ...);
  verbose && print(verbose, cdfUnits[1]);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the CEL intensities by units
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- getUnitIntensities(df, units=cdfUnits, ...);
  verbose && str(verbose, res[1]);

  verbose && exit(verbose);

  res;
}, private=TRUE)





############################################################################
# HISTORY:
# 2008-09-03
# o Created from MultiArrayUnitModel.R.
############################################################################
