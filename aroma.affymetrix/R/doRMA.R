###########################################################################/**
# @RdocDefault doRMA
# @alias doRMA.AffymetrixCelSet
#
# @title "Robust Multichip Analysis (RMA)"
#
# \description{
#  @get "title" based on [1].
#  The algorithm is processed in bounded memory, meaning virtually
#  any number of arrays can be analyzed on also very limited computer
#  systems.
#  The method replicates the results of @see "affyPLM::fitPLM"
#  (package \pkg{affyPLM}) with great precision.
# }
#
# \usage{
#   @usage doRMA,AffymetrixCelSet
#   @usage doRMA,default
# }
#
# \arguments{
#  \item{csR, dataSet}{An @see "AffymetrixCelSet" (or the name of an @see "AffymetrixCelSet").}
#  \item{arrays}{A @integer @vector specifying the subset of arrays
#   to process.  If @NULL, all arrays are considered.}
#  \item{flavor}{A character string specifying what model fitting algorithm to be used, cf. @see "RmaPlm".}
#  \item{uniquePlm}{If @TRUE, the log-additive probe-summarization model
#   is done on probeset with \emph{unique} sets of probes.
#   If @FALSE, the summarization is done on "as-is" probesets as
#   specified by the CDF.}
#  \item{drop}{If @TRUE, the summaries are returned, otherwise
#   a named @list of all intermediate and final results.}
#  \item{verbose}{See @see "Verbose".}
#  \item{...}{Additional arguments used to set up @see "AffymetrixCelSet" (when argument \code{dataSet} is specified).}
# }
#
# \value{
#   Returns a named @list, iff \code{drop == FALSE}, otherwise
#   only @see "ChipEffectSet" object.
# }
#
# \references{
#  [1] Irizarry et al.
#      \emph{Summaries of Affymetrix GeneChip probe level data}.
#      NAR, 2003, 31, e15.\cr
# }
#
# @author "HB"
#*/###########################################################################
setMethodS3("doRMA", "AffymetrixCelSet", function(csR, arrays=NULL, flavor=c("affyPLM", "oligo"), uniquePlm=FALSE, drop=TRUE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csR':
  className <- "AffymetrixCelSet";
  if (!inherits(csR, className)) {
    throw(sprintf("Argument 'csR' is not a %s: %s", className, class(csR)[1]));
  }

  # Argument 'arrays':
  if (!is.null(arrays)) {
    throw("Not supported. Argument 'arrays' should be NULL.");
    arrays <- Arguments$getIndices(arrays, max=length(csR));
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'uniquePlm':
  uniquePlm <- Arguments$getLogical(uniquePlm);

  # Argument 'drop':
  drop <- Arguments$getLogical(drop);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "RMA");
  verbose && cat(verbose, "Arguments:");
  arraysTag <- seqToHumanReadable(arrays);
  verbose && cat(verbose, "arrays:");
  verbose && str(verbose, arraysTag);
  verbose && cat(verbose, "Fit PLM on unique probe sets: ", uniquePlm);

  # Backward compatibility
  ram <- list(...)$ram;
  if (!is.null(ram)) {
    ram <- Arguments$getDouble(ram, range=c(0,Inf));
    verbose && cat(verbose, "ram: ", ram);
    warning("Argument 'ram' of doRMA() is deprecated. Instead use setOption(aromaSettings, \"memory/ram\", ram).");
    oram <- setOption(aromaSettings, "memory/ram", ram);
    on.exit({
      setOption(aromaSettings, "memory/ram", oram);
    });
  }

  verbose && cat(verbose, "Data set");
  verbose && print(verbose, csR);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(arrays)) {
    verbose && enter(verbose, "RMA/Extracting subset of arrays");
    csR <- extract(csR, arrays, onDuplicates="error");
    verbose && cat(verbose, "Data subset");
    verbose && print(verbose, csR);
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check if the final results are already available?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (drop) {
    verbose && enter(verbose, "Checking whether final results are available or not");

    # The name, tags and chip type and array names of the results to look for
    dataSet <- getFullName(csR);
    cdf <- getCdf(csR);
    chipType <- getChipType(cdf, fullname=FALSE);

    # The fullnames of all arrays that should exist
    fullnames <- getFullNames(csR);

    # RMA tags
    tags <- c("RBC");
    tags <- c(tags, "QN");
    tags <- c(tags, "RMA");
    if (flavor != "affyPLM") tags <- c(tags, flavor);

    # Try to load the final TCN data set
    ces <- tryCatch({
      cesT <- ChipEffectSet$byName(dataSet, tags=tags, chipType=chipType);
      extract(cesT, fullnames, onMissing="error", onDuplicates="error");
    }, error=function(ex) { NULL });

    # Done?
    if (!is.null(ces)) {
      verbose && cat(verbose, "Already done.");
      verbose && print(verbose, ces);
      verbose && exit(verbose);
      verbose && exit(verbose);
      return(ces);
    } # if (!is.null(ces))

    verbose && exit(verbose);
  } # if (drop)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # RMA
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # List of objects to be returned
  res <- list();
  if (!drop) {
    res <- c(res, list(csR=csR));
  }

  verbose && enter(verbose, "RMA/Background correction (normal & exponential mixture model)");
  bc <- RmaBackgroundCorrection(csR);
  verbose && print(verbose, bc);
  csB <- process(bc, verbose=verbose);
  verbose && print(verbose, csB);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(bc=bc, csB=csB));
  }

  # Clean up
  # Not needed anymore
  csR <- bc <- NULL;
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "RMA/Rank-based quantile normalization (PM-only)");
  qn <- QuantileNormalization(csB, typesToUpdate="pm");
  verbose && print(verbose, qn);
  csN <- process(qn, verbose=verbose);
  verbose && print(verbose, csN);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(qn=qn, csN=csN));
  }
  # Clean up
  # Not needed anymore
  csB <- qn <- NULL;
  gc <- gc();
  verbose && print(verbose, gc);


  verbose && enter(verbose, "RMA/Probe summarization (log-additive model)");
  verbose && cat(verbose, "Fit PLM on unique probe sets: ", uniquePlm);

  if (uniquePlm) {
    verbose && enter(verbose, "Probe-summarization using a \"unique\" CDF requested");

    verbose && enter(verbose, "Getting \"unique\" CDF (with non-unique probes dropped)")
    cdf <- getCdf(csN);
    verbose && cat(verbose, "CDF:");
    verbose && print(verbose, cdf);
    cdfU <- getUniqueCdf(cdf, verbose=less(verbose, 5));
    verbose && cat(verbose, "CDF with non-unique probes dropped:");
    verbose && print(verbose, cdfU);
    verbose && exit(verbose)

    if (equals(cdfU, cdf)) {
      verbose && cat(verbose, "The \"unique\" CDF equals the original CDF: Skipping.");
    } else {
      csNU <- convertToUnique(csN, verbose=verbose);
      verbose && print(verbose, csNU);
      csN <- csNU;
    }
    verbose && exit(verbose);
  }

  plm <- RmaPlm(csN, flavor=flavor);
  verbose && print(verbose, plm);
  if (length(findUnitsTodo(plm)) > 0) {
    units <- fit(plm, verbose=verbose);
    verbose && str(verbose, units);
    # Not needed anymore
    units <- NULL;
  }
  verbose && print(verbose, gc);
  ces <- getChipEffectSet(plm);
  verbose && print(verbose, ces);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(ces=ces, plm=plm));
  }

  # Clean up
  # Not needed anymore
  plm <- csN <- NULL;
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  # Return only the final output data set?
  if (drop) {
    res <- ces;
  }

  res;
}) # doRMA()


setMethodS3("doRMA", "default", function(dataSet, ..., verbose=FALSE) {
  .require <- require
  .require("aroma.affymetrix") || throw("Package not loaded: aroma.affymetrix")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "RMA");

  verbose && enter(verbose, "RMA/Setting up CEL set");
  csR <- AffymetrixCelSet$byName(dataSet, ..., verbose=less(verbose, 50),
                                                  .onUnknownArgs="ignore");
  verbose && print(verbose, csR);
  verbose && exit(verbose);

  res <- doRMA(csR, ..., verbose=verbose);

  # Clean up
  # Not needed anymore
  csR <- NULL;
  gc <- gc();

  verbose && exit(verbose);

  res;
})


############################################################################
# HISTORY:
# 2013-07-03
# o Added argument 'flavor' or doRMA().
# 2013-05-02
# o Removed argument 'ram' in favor of aroma option 'memory/ram'.
# 2013-04-29
# o SPEEDUP: Now doRMA() returns much quicker, iff already done.
# o DOCUMENTATION: Added help("doRMA").
# 2011-04-04
# o Added argument 'drop' to doRMA().  If FALSE, all intermediate data
#   sets and models are returned in a named list, otherwise only the
#   final data set.
# o Added argument 'uniquePlm' to doRMA().
# 2010-06-16
# o Created from doCRMAv1.R.
############################################################################
