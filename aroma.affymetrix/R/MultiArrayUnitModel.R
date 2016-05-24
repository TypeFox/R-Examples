###########################################################################/**
# @RdocClass MultiArrayUnitModel
#
# @title "The MultiArrayUnitModel class"
#
# \description{
#  @classhierarchy
#
#  This abstract class represents a unit model that fits one model per unit
#  based on signals for all arrays in the data set.
#  The nature of a multi-array unit model is that all arrays must be
#  available at the time of the fit and the estimated parameters will
#  depend on the data from all arrays.
#  Thus, if the signals in one array changes the model has to be refitted.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "UnitModel".}
#   \item{listOfPriors}{A @list of priors to be used when fitting
#    the model.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("MultiArrayUnitModel", function(..., listOfPriors=NULL) {
  this <- extend(UnitModel(...), "MultiArrayUnitModel");
  if (!is.null(listOfPriors)) {
    this <- setListOfPriors(this, listOfPriors);
  }
  this;
}, abstract=TRUE)


setMethodS3("validate", "MultiArrayUnitModel", function(this, ...) {
  ds <- getDataSet(this);
  if (is.null(ds))
    return(invisible(TRUE));

  if (length(ds) < 2) {
    priors <- getListOfPriors(this);
    hasPriors <- (!is.null(priors));
    if (!hasPriors) {
      throw("This ", class(this)[1], " requires at least 2 arrays: ",
                                                            length(ds));
    }
  }

  invisible(TRUE);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getFitUnitGroupFunction
#
# @title "Static method to get the low-level function that fits the PLM"
#
# \description{
#  @get "title".
#  Any subclass model must provide this method, which should return
#  a @function that accepts an IxK @matrix.
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
setMethodS3("getFitUnitGroupFunction", "MultiArrayUnitModel", abstract=TRUE, static=TRUE, private=TRUE);



setMethodS3("getFitUnitFunction", "MultiArrayUnitModel", function(this, ...) {
  # Get the fit function for a single set of intensities
  fitfcn <- getFitUnitGroupFunction(this, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the one for all blocks in a unit
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (this$probeModel == "pm-mm") {
    fitUnit <- function(unit, ...) {
      lapply(unit, FUN=function(group) {
        y <- .subset2(group, 1); # Get intensities
        y <- y[1,,] - y[2,,];  # PM-MM
        fitfcn(y, ...);
      })
    }
  } else if (this$probeModel == "min1(pm-mm)") {
    fitUnit <- function(unit, ...) {
      lapply(unit, FUN=function(group) {
        y <- .subset2(group, 1); # Get intensities
        y <- y[1,,] - y[2,,];  # PM-MM
        y[y < 1] <- 1;       # min1(PM-MM)=min(PM-MM,1)
        fitfcn(y, ...);
      })
    }
  } else if (this$probeModel == "pm+mm") {
    fitUnit <- function(unit, ...) {
      lapply(unit, FUN=function(group) {
        y <- .subset2(group, 1); # Get intensities
        y <- y[1,,] + y[2,,];  # PM+MM
        fitfcn(y, ...);
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
        fitfcn(y, ...);
      })
    }
  }

  fitUnit;
}, private=TRUE) # getFitUnitFunction()





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
#     @see "AffymetrixCelSet" class.}
# }
#
# \value{
#  Returns the @list structure that \code{readUnits()} of
#  @see "AffymetrixCelSet" returns.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("readUnits", "MultiArrayUnitModel", function(this, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  ds <- getDataSet(this);
  verbose && enter(verbose, "Reading probe intensities from ", length(ds), " arrays");

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
  res <- getUnitIntensities(ds, units=cdfUnits, ...);
  verbose && str(verbose, res[1]);

  verbose && exit(verbose);


  res;
}, private=TRUE)




setMethodS3("getListOfPriors", "MultiArrayUnitModel", function(this, ...) {
  this$.listOfPriors;
})

setMethodS3("setListOfPriors", "MultiArrayUnitModel", function(this, sets, ...) {
  # Arguments 'sets':
  if (!is.list(sets)) {
    throw("Argument 'sets' is not a list: ", mode(sets)[1]);
  }

  this$.listOfPriors <- sets;

  invisible(this);
})


setMethodS3("readPriorsByUnits", "MultiArrayUnitModel", function(this, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  priorsList <- getListOfPriors(this);
  nbrOfPriors <- length(priorsList);

  # Nothing to do?
  if (nbrOfPriors == 0)
    return(NULL);


  verbose && enter(verbose, "Reading prior parameters by unit");

  verbose && enter(verbose, "Reading prior parameters");
  verbose && cat(verbose, "Number of parameter sets: ", nbrOfPriors);
  res <- list();
  for (kk in seq_len(nbrOfPriors)) {
    priors <- priorsList[[kk]];
    verbose && enter(verbose, sprintf("Prior set #%d ('%s') of %d", kk, getName(priors), length(priorsList)));
    values <- readUnits(priors, units=units, ..., verbose=less(verbose, 1));
    verbose && str(verbose, values[1]);
    res[[kk]] <- values;
    # Not needed anymore
    values <- priors <- NULL;
    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  if (length(priorsList) > 0) {
    unit <- vector("list", nbrOfPriors);
    names(unit) <- names(priorsList);
    unit0 <- unit;

    verbose && enter(verbose, "Merging prior parameters unit by unit");
    unitNames <- names(res[[1]]);
    res2 <- vector("list", length(unitNames));
    names(res2) <- unitNames;
    for (uu in seq_along(res2)) {
      unit <- unit0;
      for (kk in seq_len(nbrOfPriors)) {
        unit[[kk]] <- res[[kk]][[uu]];
      } # for (kk ...)
      res2[[uu]] <- unit;
    } # for (uu ...)
    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  res2;
}, private=TRUE)





############################################################################
# HISTORY:
# 2012-01-14
# o Now the fit function returned by getFitUnitFunction() for
#   MultiArrayUnitModel passes '...'.
# o Now validate() for MultiArrayUnitModel accepts single-array data sets,
#   iff priors are set.
# o Added argument 'listOfPriors' to MultiArrayUnitModel().
# 2008-12-08
# o Added readPriorsByUnits().
# 2008-09-03
# o Added getFitUnitGroupFunction() model, which is a better name than
#   getFitFunction().
# o Created from ProbeLevelModel.R.
############################################################################
