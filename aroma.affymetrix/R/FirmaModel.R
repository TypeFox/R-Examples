###########################################################################/**
# @RdocClass FirmaModel
#
# @title "The FirmaModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the FIRMA (Finding Isoforms using RMA) alternative
#  splicing model.
#
# }
#
# @synopsis
#
# \arguments{
#   \item{rmaPlm}{An @RmaPlm object.}
#   \item{summaryMethod}{A @character specifying what summarization method should be used.}
#   \item{operateOn}{A @character specifying what statistic to operate on.}
#   \item{...}{Arguments passed to constructor of @see "UnitModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "KS, HB"
#*/###########################################################################
setConstructorS3("FirmaModel", function(rmaPlm=NULL, summaryMethod=c("median", "upperQuartile", "max"), operateOn=c("residuals", "weights"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'rmaPlm':
  if (!is.null(rmaPlm)) {
    rmaPlm <- Arguments$getInstanceOf(rmaPlm, "ProbeLevelModel");
    rmaPlm <- Arguments$getInstanceOf(rmaPlm, "ExonRmaPlm");

    # Assert that the RmaPlm has 'mergeGroups=TRUE'.
    params <- getParameters(rmaPlm);
    if (!params$mergeGroups) {
      throw("Cannot setup FirmaModel. The probe-level model must be for transcripts (mergeGroups=TRUE), not exons.");
    }
  }

  # Argument 'summaryMethod':
  summaryMethod <- match.arg(summaryMethod);

  # Argument 'operateOn':
  operateOn <- match.arg(operateOn);

  extend(UnitModel(...), "FirmaModel",
     .plm = rmaPlm,
     summaryMethod = summaryMethod,
     operateOn = operateOn,
     "cached:.fs" = NULL,
     "cached:.paFile" = NULL,
     "cached:.chipFiles" = NULL,
     "cached:.lastPlotData" = NULL
   );
})


setMethodS3("getAsteriskTags", "FirmaModel", function(this, collapse=NULL, ...) {
  # Returns 'U' (but allow for future extensions)
  tags <- NextMethod("getAsteriskTags", collapse=NULL);
  tags[1] <- "FIRMA";

  # Append parameter tags

  # Create one tag for (summaryMethod, operateOn)
  # EP on 2007-12-08:
  #  if operateOn==weights
  #         and summaryMethod=upperQuartile, then add to FIRMA tag=uqwt
  #         and summaryMethod=median, then add to FIRMA tag=medwt
  #         and summaryMethod=max, then add to FIRMA tag=maxwt
  #  if operateOn==residuals (all new)
  #         and summaryMethod=mean, then add to FIRMA tag=meanres
  #         and summaryMethod=median, then add to FIRMA tag=medres
  codes <- c("upperQuartile"="uq", "mean"="mean", "median"="med", "max"="max");
  smCode <- codes[this$summaryMethod];
  codes <- c("residuals"="res", "weights"="wt");
  ooCode <- codes[this$operateOn];
  smooTag <- paste(smCode, ooCode, sep="");
  tags <- c(tags, smooTag);

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)


setMethodS3("getTags", "FirmaModel", function(this, collapse=NULL, ...) {
  # "Pass down" tags from the "input" data set, which "happens to be"
  # the chip types, i.e. getChipEffectSet(plm).  This is why we can't
  # call NextMethod() here, because that is using getDataSet(plm).
  # /HB 2007-12-13
  plm <- getPlm(this);
  ces <- getChipEffectSet(plm);
  inputTags <- getTags(ces);

  # Get class specific tags
  tags <- this$.tags;

  # In case this$.tags is not already split
  tags <- Arguments$getTags(tags, collapse=NULL);

  # Expand asterisk tags
  if (any(tags == "*")) {
    tags[tags == "*"] <- getAsteriskTags(this, collapse=",");
  }

  # Combine input tags and local tags
  tags <- c(inputTags, tags);

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    if (length(tags) > 0)
      tags <- unlist(strsplit(tags, split=","));
  }

  # No tags?
  if (length(tags) == 0)
    tags <- NULL;

  tags;
})


setMethodS3("getPlm", "FirmaModel", function(this, ...) {
  this$.plm;
})

setMethodS3("getDataSet", "FirmaModel", function(this, ...) {
  getDataSet(this$.plm);
})

setMethodS3("getCdf", "FirmaModel", function(this, ...) {
  getCdf(this$.plm);
})

setMethodS3("getName", "FirmaModel", function(this, ...) {
  getName(this$.plm, ...);
})



setMethodS3("as.character", "FirmaModel", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  ds <- getDataSet(this);
  s <- c(s, sprintf("Data set: %s", getName(ds)));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(ds))));
  s <- c(s, sprintf("Input tags: %s", getTags(this$.plm, collapse=",")));
  s <- c(s, sprintf("Output tags: %s", getTags(this, collapse=",")));
  s <- c(s, sprintf("Parameters: %s", getParametersAsString(this)));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)


setMethodS3("calculateWeights", "FirmaModel", function(this, ...) {
  calculateWeights(this$.plm, ...);
}, protected=TRUE)



setMethodS3("getFileSetClass", "FirmaModel", function(static, ...) {
  FirmaSet;
}, static=TRUE, private=TRUE)


setMethodS3("getRootPath", "FirmaModel", function(this, ...) {
  "firmaData";
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getFirmaSet
# @aliasmethod getFirmaScores
#
# @title "Gets the set of FIRMA results for this model"
#
# \description{
#  @get "title".
#  There is one chip-effect file per array.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "FirmaSet" object.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFirmaSet", "FirmaModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  fs <- this$.fs;
  if (!is.null(fs))
    return(fs);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Let the parameter object know about the CDF structure, because we
  # might use a modified version of the one in the CEL header.
  ds <- getDataSet(this);
  if (length(ds) == 0)
    throw("Cannot create FIRMA results file. The CEL set is empty.");

  verbose && enter(verbose, "Getting FIRMA results set from data set");
  # Inherit the (monocell) CDF
  cdf <- getCdf(ds);
  cdfMono <- getMonocellCdf(cdf);

  # Gets the Class object
  clazz <- getFileSetClass(this);
  fs <- clazz$fromDataSet(dataSet=ds, path=getPath(this), cdf=cdfMono,
         pattern=",FIRMAscores[.](c|C)(e|E)(l|L)$", verbose=less(verbose));
  verbose && exit(verbose);

  # Store in cache
  this$.fs <- fs;

  fs;
})

setMethodS3("getFirmaScores", "FirmaModel", function(this, ...) {
  getFirmaSet(this, ...);
}, protected=TRUE)


###########################################################################/**
# @RdocMethod calculateResidualSet
# @aliasmethod calculateResiduals
#
# @title "Gets the set of residuals corresponding to the PLM"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "ResidualSet" object.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("calculateResidualSet", "FirmaModel", function(this, ...) {
  calculateResidualSet(this$.plm, ...)
}, protected=TRUE)




###########################################################################/**
# @RdocMethod getFitUnitGroupFunction
#
# @title "Static method to get the low-level function that fits the PLM"
#
# \description{
#  @get "title".
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
setMethodS3("getFitUnitGroupFunction", "FirmaModel", function(this, ...) {
  if(this$operateOn == "weights") {
    if (this$summaryMethod == "upperQuartile") {
      fitfcn <- function(y, ...) {
        ## J <- length(y)
        list(
          intensities=2^(1 - colQuantiles(y, probs=0.75)),
          stdvs=1, pixels=1
        )
      }
    }
    else if (this$summaryMethod == "median") {
      fitfcn <- function(y, ...) {
        ## J <- length(y)
        ## 1 - median(y)
        list(
          intensities=2^(1-colMedians(y)),
          stdvs=1, pixels=1
        )
      }
    }
    else if (this$summaryMethod == "max") {
      fitfcn <- function(y, ...) {
        ## J <- length(y)
        list(
          intensities=2^(1-colMaxs(y)),
          stdvs=1, pixels=1
        )
      }
    }
    else {
      fitfcn <- function(y, ...) {
        ## J <- length(y)
        list(intensities=1, stdvs=1, pixels=1)
      }
    }
  } else {
    if (this$summaryMethod == "median") {
      fitfcn <- function(y, ...) {
        ## J <- length(y)
        list(intensities=2^colMedians(y), stdvs=1, pixels=1)
      }
    }
    else if (this$summaryMethod == "mean") {
      fitfcn <- function(y, ...) {
        ## J <- length(y)
        list(intensities=2^colMeans(y), stdvs=1, pixels=1)
      }
    }
    else {
      fitfcn <- function(y, ...) {
        ## J <- length(y)
        list(intensities=1, stdvs=1, pixels=1)
      }
    }
  }
  fitfcn;
}, protected=TRUE)


setMethodS3("getFitUnitFunction", "FirmaModel", function(this, ...) {

  fitfcn <- getFitUnitGroupFunction(this, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Function to fit all groups (exons) within a unit
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(identical(this$operateOn, "residuals")) {
    fitUnit <- function(unit, ...) {
      u.mad <- mad(log2(unlist(unit, use.names=FALSE)), center=0);
      lapply(unit, FUN=function(group) {
        if (length(group) > 0) {
          y <- .subset2(group, 1); # Get intensities
        } else {
          y <- NULL;
        }
        y <- log2(y);              # convert to log scale
        fitfcn(y/u.mad);
      })
    } # fitUnit()
  } else {
    fitUnit <- function(unit, ...) {
      lapply(unit, FUN=function(group) {
        if (length(group) > 0) {
          y <- .subset2(group, 1); # Get intensities
        } else {
          y <- NULL;
        }
        y <- log2(y);              # convert to log scale
        fitfcn(y);
      })
    } # fitUnit()
  }

  fitUnit;
}, private=TRUE)



###########################################################################/**
# @RdocMethod findUnitsTodo
#
# @title "Identifies non-fitted units"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer @vector of unit indices.
# }
#
# \seealso{
#   Internally this methods calls the same method for the
#   @see "ChipEffectSet" class.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("findUnitsTodo", "FirmaModel", function(this, ...) {
  fs <- getFirmaScores(this);
  findUnitsTodo(fs, ...);
}, private=TRUE)



###########################################################################/**
# @RdocMethod fit
#
# @title "Estimates the model parameters"
#
# \description{
#  @get "title" for all or a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be fitted.
#     If @NULL, all units are considered.
#     If \code{remaining}, only non-fitted units are considered.
#   }
#   \item{...}{Arguments passed to \code{readUnits()}.}
#   \item{ram}{A @double indicating if more or less units should
#     be loaded into memory at the same time.}
#   \item{force}{If @TRUE, already fitted units are re-fitted, and
#     cached data is re-read.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @integer @vector of indices of the units fitted,
#  or @NULL if no units was (had to be) fitted.
# }
#
# \details{
#  All estimates are stored to file.
#
#  The non-array specific parameter estimates together with standard deviation
#  estimates and convergence information are stored in one file.
#
#  The parameter estimates specific to each array, typically "chip effects",
#  are stored in array specific files.
#
#   Data set specific estimates [L = number of probes]:
#    phi [L doubles] (probe affinities), sd(phi) [L doubles],
#    isOutlier(phi) [L logicals]
#
#   Algorithm-specific results:
#    iter [1 integer], convergence1 [1 logical], convergence2 [1 logical]
#    dTheta [1 double]
#    sd(eps) - [1 double] estimated standard deviation of the error term
#
#   Array-specific estimates [K = nbr of arrays]:
#    theta [K doubles] (chip effects), sd(theta) [K doubles],
#    isOutlier(theta) [K logicals]
#
#   For each array and each unit group, we store:
#     1 theta, 1 sd(theta), 1 isOutlier(theta), i.e. (float, float, bit)
#   => For each array and each unit (with \eqn{G_j} groups), we store:
#     \eqn{G_j} theta, \eqn{G_j} sd(theta), \eqn{G_j} isOutlier(theta),
#   i.e. \eqn{G_j}*(float, float, bit).
#   For optimal access we store all thetas first, then all sd(theta) and the
#   all isOutlier(theta).
#   To keep track of the number of groups in each unit, we have to have a
#   (unit, ngroups) map.  This can be obtained from getUnitNames() for the
#   AffymetrixCdfFile class.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("fit", "FirmaModel", function(this, units="remaining", ..., ram=NULL, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the some basic information about this model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getDataSet(this$.plm);
  cdf <- getCdf(ds);
  if (this$operateOn == "weights") {
    ws <- calculateWeights(this, verbose = verbose)
  } else {
    ws <- calculateResidualSet(this, verbose = verbose)
  }
  nbrOfArrays <- length(ds);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  doRemaining <- FALSE;
  if (is.null(units)) {
  } else if (is.numeric(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  } else if (identical(units, "remaining")) {
    doRemaining <- TRUE;
  } else {
    throw("Unknown mode of argument 'units': ", mode(units));
  }

  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting model of class ", class(this)[1]);

  verbose && print(verbose, this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify units to be fitted
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    units <- 1:nbrOfUnits;
  } else if (doRemaining) {
    verbose && enter(verbose, "Identifying non-estimated units")
    units <- findUnitsTodo(this, verbose=less(verbose));
    nbrOfUnits <- length(units);
    verbose && exit(verbose);
  } else {
    # Fit only unique units
    units <- unique(units);
    nbrOfUnits <- length(units);
  }
  verbose && printf(verbose, "Getting FIRMA results for %d units.\n", nbrOfUnits);

  # Identify which of the requested units have *not* already been estimated
  # for all arrays
  if (!doRemaining) {
    if (force) {
      verbose && printf(verbose, "All of these are forced to be fitted.\n");
    } else {
      units <- findUnitsTodo(this, units=units, verbose=less(verbose));
      nbrOfUnits <- length(units);
      verbose && printf(verbose, "Out of these, %d units need to be fitted.\n", nbrOfUnits);
    }
  }

  # Nothing more to do?
  if (nbrOfUnits == 0)
    return(NULL);

  fs <- getFirmaScores(this, verbose=less(verbose));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit the model in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the model-fit function
  fitUnit <- getFitUnitFunction(this);

  # Number of units to load into memory and fit at the same time
  bytesPerChunk <- 100e6;       # 100Mb
  bytesPerUnitAndArray <- 400;  # Just a rough number; good enough?
  bytesPerUnit <- nbrOfArrays * bytesPerUnitAndArray;
  unitsPerChunk <- ram * bytesPerChunk / bytesPerUnit;
  unitsPerChunk <- as.integer(max(unitsPerChunk,1));

  idxs <- 1:nbrOfUnits;
  head <- 1:unitsPerChunk;
  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && cat(verbose, "Number units per chunk: ", unitsPerChunk);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Garbage collect
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  clearCache(ds);
  clearCache(ws);
  clearCache(fs);
  gc <- gc();
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the model-fit function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitUnit <- getFitUnitFunction(this);

  startTime <- processTime();
  timers <- list(total=0, read=0, fit=0, writeFs=0, gc=0);

  count <- 1;
  while (length(idxs) > 0) {
    tTotal <- processTime();

    verbose && enter(verbose, "Fitting chunk #", count, " of ", nbrOfChunks);
    if (length(idxs) < unitsPerChunk) {
      head <- 1:length(idxs);
    }
    uu <- idxs[head];

    verbose && cat(verbose, "Units: ");
    verbose && str(verbose, units[uu]);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get the CEL intensities by units
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tRead <- processTime();
    y <- readUnits(ws, units=units[uu], ..., force=force, cache=FALSE, verbose=less(verbose));
    timers$read <- timers$read + (processTime() - tRead);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit the model
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Calculating FIRMA scores");
    tFit <- processTime();
    fit <- lapply(y, FUN=fitUnit);
    #fit <- lapply(fit, FUN=normalizeFits);
    timers$fit <- timers$fit + (processTime() - tFit);
    y <- NULL; # Not needed anymore (to minimize memory usage)
    verbose && str(verbose, fit[1]);
    verbose && exit(verbose);

    # Garbage collection
    tGc <- processTime();
    gc <- gc();
    timers$gc <- timers$gc + (processTime() - tGc);
    verbose && print(verbose, gc);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Store FIRMA scores
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing FIRMA results");
    tWriteFs <- processTime();
    updateUnits(fs, units=units[uu], data=fit, verbose=less(verbose));
    timers$writeFs <- timers$writeFs + (processTime() - tWriteFs);
    verbose && exit(verbose);

    fit <- NULL; # Not needed anymore

    # Next chunk
    idxs <- idxs[-head];
    count <- count + 1;

    # Garbage collection
    tGc <- processTime();
    gc <- gc();
    timers$gc <- timers$gc + (processTime() - tGc);
    verbose && print(verbose, gc);

    timers$total <- timers$total + (processTime() - tTotal);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # ETA
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (verbose) {
      # Clarifies itself once in a while (in case running two in parallel).
      verbose && print(verbose, this);

      # Fraction left
      fLeft <- length(idxs) / nbrOfUnits;
      # Time this far
      lapTime <- processTime() - startTime;
      t <- Sys.time() - lapTime[3];
      printf(verbose, "Started: %s\n", format(t, "%Y%m%d %H:%M:%S"));
      # Estimated time left
      fDone <- 1-fLeft;
      timeLeft <- fLeft/fDone * lapTime;
      t <- timeLeft[3];
      printf(verbose, "Estimated time left: %.1fmin\n", t/60);
      # Estimate time to arrivale
      eta <- Sys.time() + t;
      printf(verbose, "ETA: %s\n", format(eta, "%Y%m%d %H:%M:%S"));
    }
    verbose && exit(verbose);

  } # while (length(idxs) > 0) loop over chunks

  totalTime <- processTime() - startTime;
  if (verbose) {
    nunits <- length(units);
    t <- totalTime[3];
    printf(verbose, "Total time for all units across all %d arrays: %.2fs == %.2fmin\n", nbrOfArrays, t, t/60);
    t <- totalTime[3]/nunits
    printf(verbose, "Total time per unit across all %d arrays: %.2fs/unit\n", nbrOfArrays, t);
    t <- totalTime[3]/nunits/nbrOfArrays;
    printf(verbose, "Total time per unit and array: %.3gms/unit & array\n", 1000*t);
    t <- nbrOfUnits(cdf)*totalTime[3]/nunits/nbrOfArrays;
    printf(verbose, "Total time for one array (%d units): %.2fmin = %.2fh\n", nbrOfUnits(cdf), t/60, t/3600);
    t <- nbrOfUnits(cdf)*totalTime[3]/nunits;
    printf(verbose, "Total time for complete data set: %.2fmin = %.2fh\n", t/60, t/3600);
    # Get distribution of what is spend where
    t <- lapply(timers, FUN=function(timer) unname(timer[3]));
    t <- unlist(t);
    t <- 100 * t / t["total"];
    printf(verbose, "Fraction of time spent on different tasks: Fitting: %.1f%%, Reading: %.1f%%, Writing: %.1f%%, Explicit garbage collection: %.1f%%\n", t["fit"], t["read"], t["writeFs"], t["gc"]);
  }

  ## Create checksum files
  fsZ <- getChecksumFileSet(fs)

  invisible(units);
})


############################################################################
# HISTORY:
# 2011-11-09
# o ROBUSTIFICATION: Now FirmaModel(plm) asserts that 'plm' is an
#   ExonRmaPlm object (before any ProbeLevelModel would work), and
#   that the PLM is setup for transcripts (mergeGroups=TRUE), not exons.
# o Now the fit functions utilizes the 'matrixStats' package.
# o CLEANUP: Dropped unused expressions from the fit functions.
# 2008-05-31
# o Removed an obsolete debug print() statement.
# 2008-04-09 [HB]
# o Added calculateResidualSet() and made calculateResiduals() call it.
# 2008-02-28 [HB]
# o Added getFirmaSet() (to replace getFirmaScores()?).
# 2007-12-13 [HB]
# o Updated getAsteriskTags() and getTags().
# 2007-12-10 [HB]
# o Now getFirmaScores() of FirmaModel infers the monocell CDF from
#   the CDF of the input data set and uses that when retrieving the
#   chip-effect CEL set.  In other words, if the CDF is overridden for
#   the input data set, it will also be overridden (with the corresponding
#   monocell CDF) in the chip-effect set.  Before the monocell CDF was
#   always inferred from the CEL header, if the CEL file existed.
# 2007-12-07 [HB]
# o UPDATE: Renamed the "root" directory of FirmaModel to firmaData/
#   (formely modelFirmaModel/).
# o BUG FIX: Tags from the input data set of FirmaModel were lost.
# o Added getAsteriskTag().
# 2007-07-01 [HB]
# o Added 'cache=FALSE' to findUnitsTodo() in fit() so that the results are
#   not stored in memory for every file.
# o Now clearCache() clears the private field '.fs'.
# o Added more code to fit() to better clean up the memory.
# 2007-02-09
# o Created (based largely on ProbeLevelModel.R).
############################################################################
