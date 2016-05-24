###########################################################################/**
# @RdocDefault doFIRMA
# @alias doFIRMA.AffymetrixCelSet
#
# @title "Finding Isoforms using Robust Multichip Analysis (FIRMA)"
#
# \description{
#  @get "title" based on [1].
# }
#
# \usage{
#   @usage doFIRMA,AffymetrixCelSet
#   @usage doFIRMA,default
# }
#
# \arguments{
#  \item{csR, dataSet}{An @see "AffymetrixCelSet" (or the name of an @see "AffymetrixCelSet").}
#  \item{...}{Additional arguments passed to @see "FirmaModel",
#             and to set up @see "AffymetrixCelSet" (when
#             argument \code{dataSet} is specified).}
#  \item{flavor}{A @character string specifying the flavor of FIRMA to use.}
#  \item{drop}{If @TRUE, the FIRMA scores are returned, otherwise
#   a named @list of all intermediate and final results.}
#  \item{verbose}{See @see "Verbose".}
# }
#
# \value{
#   Returns a named @list, iff \code{drop == FALSE}, otherwise
#   only @see "FirmaSet" object (containing the FIRMA scores).
# }
#
# \section{Using a custom exon-by-transcript CDF}{
#   It is strongly recommended to use a custom CDF, e.g. "core",
#   "extended" or "full" [1].  To use a custom CDF, set it before
#   calling this method, i.e. \code{setCdf(csR, cdf)}.
#   Do not set the standard "non-supported" Affymetrix CDF
#   (see also Section 'Flavors').
# }
#
# \section{Flavors}{
#   If \code{flavor == "v1b"} (default), then the standard
#   "non-supported" Affymetrix CDF is used for background correction
#   and the quantile normalization steps, and the custom CDF
#   is used for the probe summarization and everything that follows.
#   The advantage of this flavor is that those two first preprocessing
#   steps will remain the same if one later changes to a different
#   custom CDF.
#
#   If \code{flavor == "v1a"}, then the custom CDF is used throughout
#   all steps of FIRMA, which means that if one changes the custom CDF
#   all steps will be redone.
# }
#
# \references{
#  [1] E. Purdom, K. Simpson, M. Robinson, J. Conboy, A. Lapuk & T.P. Speed,
#      \emph{FIRMA: a method for detection of alternative splicing from
#            exon array data}, Bioinformatics, 2008.\cr
# }
#
# @author "HB"
#*/###########################################################################
setMethodS3("doFIRMA", "AffymetrixCelSet", function(csR, ..., flavor=c("v1b", "v1a"), drop=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csR':
  className <- "AffymetrixCelSet";
  if (!inherits(csR, className)) {
    throw(sprintf("Argument 'csR' is not a %s: %s", className, class(csR)[1]));
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'drop':
  drop <- Arguments$getLogical(drop);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);




  verbose && enter(verbose, "FIRMA");
  verbose && cat(verbose, "Flavor: ", flavor);
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, list(...));

  # Backward compatibility
  ram <- list(...)$ram;
  if (!is.null(ram)) {
    ram <- Arguments$getDouble(ram, range=c(0,Inf));
    verbose && cat(verbose, "ram: ", ram);
    warning("Argument 'ram' of doFIRMA() is deprecated. Instead use setOption(aromaSettings, \"memory/ram\", ram).");
    oram <- setOption(aromaSettings, "memory/ram", ram);
    on.exit({
      setOption(aromaSettings, "memory/ram", oram);
    });
  }

  # List of objects to be returned
  res <- list();
  if (!drop) {
    res <- c(res, list(csR=csR));
  }

  verbose && cat(verbose, "Data set:");
  verbose && print(verbose, csR);

  # Tag for the custom CDF
  customCdfTags <- getTags(cdf);
  customCdfTag <- customCdfTags[1];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (0) Use the standard CDF for RMA preprocessing(?)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdfTag <- customCdfTag;
  if (is.element(flavor, c("v1b"))) {
    verbose && cat(verbose, "Detected flavor: ", flavor);
    verbose && enter(verbose, "FIRMA/Retrieving the standard CDF");

    cdf <- getCdf(csR);
    chipTypeS <- getChipType(cdf, fullname=FALSE);
    cdfS <- AffymetrixCdfFile$byChipType(chipTypeS);

    if (!equals(cdfS, cdf)) {
      verbose && print(verbose, cdf);
      setCdf(csR, cdfS);
      on.exit({
        # Make sure to undo, e.g. if interrupted.
        setCdf(csR, cdf);
      });
      cdfTag <- NULL;
    } else {
      verbose && cat(verbose, "Same as the custom CDF. Skipping.");
    }

    # Not needed anymore
    # Not needed anymore
    chipTypeS <- NULL;

    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (1) Background correction
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "FIRMA/RMA/Background correction");
  verbose && cat(verbose, "Using CDF: ", getFullName(getCdf(csR)));

  bc <- RmaBackgroundCorrection(csR, tags=c("*", cdfTag));
  verbose && print(verbose, bc);
  csB <- process(bc, verbose=verbose);
  verbose && print(verbose, csB);

  if (!drop) {
    res <- c(res, list(bc=bc, csB=csB));
  }

  # Not needed anymore
  bc <- NULL;
  gc <- gc();
  verbose && print(verbose, gc);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (2) Quantile normalization
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "FIRMA/RMA/Rank-based quantile normalization (PM-only)");
  verbose && cat(verbose, "Using CDF: ", getFullName(getCdf(csB)));

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

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (3) Quantile normalization
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "FIRMA/RMA/Probe summarization of transcripts (log-additive model)");

  if (is.element(flavor, c("v1b"))) {
    if (!equals(cdfS, cdf)) {
      # Always use the custom CDF in what follows
      setCdf(csN, cdf);
    }
    # Not needed anymore
    # Not needed anymore
    cdfS <- NULL;
  }
  verbose && cat(verbose, "Using CDF: ", getFullName(getCdf(csN)));

  # PLM tags
  cdfTag <- customCdfTag;
  plmTags <- c("*", cdfTag);
  verbose && cat(verbose, "Tag added: ", cdfTag);

  plm <- ExonRmaPlm(csN, mergeGroups=TRUE, tags=plmTags);
  verbose && print(verbose, plm);

  units <- fit(plm, verbose=verbose);
  verbose && cat(verbose, "Units fitted:");
  verbose && str(verbose, units);

  if (!drop) {
    res <- c(res, list(plm=plm));
  }

  # Clean up
  # Not needed anymore
  cdfTag <- plmTags <- units <- NULL;
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (4) Alternative splicing analysis
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "FIRMA/Alternative Splicing Analysis");

  # Setup FIRMA model
  firma <- FirmaModel(plm, ..., .onUnknownArgs="warning");
  verbose && print(verbose, firma);

  units <- fit(firma, verbose=verbose);
  verbose && cat(verbose, "Units processed:");
  verbose && str(verbose, units);

  fs <- getFirmaScores(firma, verbose=less(verbose, 10));
  verbose && str(verbose, "FIRMA scores:");
  verbose && print(verbose, fs);

  if (!drop) {
    res <- c(res, list(firma=firma, fs=fs));
  }

  # Return only the final output data set?
  if (drop) {
    res <- fs;
  }

  # Clean up
  # Not needed anymore
  firma <- fs <- NULL;
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # doFIRMA()



setMethodS3("doFIRMA", "default", function(dataSet, ..., verbose=FALSE) {
  .require <- require
  .require("aroma.affymetrix") || throw("Package not loaded: aroma.affymetrix")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "FIRMA");

  verbose && enter(verbose, "FIRMA/Setting up CEL set");
  csR <- AffymetrixCelSet$byName(dataSet, ..., verbose=less(verbose, 50),
                                                  .onUnknownArgs="ignore");
  verbose && print(verbose, csR);
  verbose && exit(verbose);

  res <- doFIRMA(csR, ..., verbose=verbose);

  # Clean up
  # Not needed anymore
  csR <- NULL;
  gc <- gc();

  verbose && exit(verbose);

  res;
}) # doFIRMA()


############################################################################
# HISTORY:
# 2013-06-02
# o BUG FIX: doFIRMA() would give <simpleError in UseMethod("setCdf"):
#   no applicable method for 'setCdf' applied to an object of class
#   "NULL"> - a bug introduced in v2.9.3.
# 2013-05-02
# o Removed argument 'ram' in favor of aroma option 'memory/ram'.
# 2011-11-10
# o Added Rdoc comments.
# o Created from doGCRMA.R.
############################################################################
