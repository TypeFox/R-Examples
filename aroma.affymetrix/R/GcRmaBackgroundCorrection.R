###########################################################################/**
# @RdocClass GcRmaBackgroundCorrection
#
# @title "The GcRmaBackgroundCorrection class"
#
# \description{
#  @classhierarchy
#
#  This class represents the GCRMA background adjustment function.
#
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#      @see "ProbeLevelTransform".}
#   \item{indicesNegativeControl}{Locations of any negative control
#       probes (e.g., the anti-genomic controls on the human exon array).
#       If @NULL and \code{type == "affinities"}, then all non-PM probes
#       are used as the negative controls.}
#   \item{affinities}{A @numeric @vector of probe affinities, usually as
#       calculated by \code{computeAffinities()} of the
#       @see "AffymetrixCdfFile" class.}
#   \item{type}{Type (flavor) of background correction, which can
#       be either \code{"fullmodel"} (uses MMs; requires that the chip type
#       has PM/MM pairs) or \code{"affinities"} (uses probe sequence only).}
#   \item{gsbAdjust}{If @TRUE, adjustment for specific binding is done,
#       otherwise not.}
#   \item{opticalAdjust}{If @TRUE, correction for optical effect is done
#       first, utilizing @see "OpticalBackgroundCorrection".}
#   \item{gsbParameters}{Additional argument passed to the internal
#       \code{bgAdjustGcrma()} method.}
#   \item{seed}{An (optional) @integer specifying a temporary random seed
#     to be used during processing.  The random seed is set to its original
#     state when done.  If @NULL, it is not set.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \references{
#  [1] Z. Wu, R. Irizarry, R. Gentleman, F.M. Murillo & F. Spencer.
#      \emph{A Model Based Background Adjustment for Oligonucleotide
#      Expression Arrays}, JASA, 2004.\cr
# }
#
# @author "KS, HB"
#*/###########################################################################
setConstructorS3("GcRmaBackgroundCorrection", function(..., indicesNegativeControl=NULL, affinities=NULL, type=c("fullmodel", "affinities"), opticalAdjust=TRUE, gsbAdjust=TRUE, gsbParameters=NULL, seed=NULL) {

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'indicesNegativeControl':
  if (!is.null(indicesNegativeControl)) {
    indicesNegativeControl <- Arguments$getIndices(indicesNegativeControl);
  }

  # Argument 'affinities':
  if (!is.null(affinities)) {
    affinities <- Arguments$getNumerics(affinities);
  }

  # Argument 'type':
  type <- match.arg(type);

  # Argument 'opticalAdjust':
  opticalAdjust <- Arguments$getLogical(opticalAdjust);

  # Argument 'gsbAdjust':
  gsbAdjust <- Arguments$getLogical(gsbAdjust);

  # Argument 'gsbParameters':

  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed);
  }


  extend(BackgroundCorrection(..., typesToUpdate="pm"), "GcRmaBackgroundCorrection",
    .indicesNegativeControl=indicesNegativeControl,
    .affinities=affinities,
    .type=type,
    .opticalAdjust=opticalAdjust,
    .gsbAdjust=gsbAdjust,
    .gsbParameters=gsbParameters,
    .seed=seed
  );
})


setMethodS3("getParameters", "GcRmaBackgroundCorrection", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters");

  # Get parameters of this class
  params2 <- list(
    indicesNegativeControl = this$.indicesNegativeControl,
    affinities = this$.affinities,
    type = this$.type,
    opticalAdjust = this$.opticalAdjust,
    gsbAdjust = this$.gsbAdjust,
    gsbParameters = this$.gsbParameters,
    seed = this$.seed
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, protected=TRUE)


setMethodS3("calculateAffinities", "GcRmaBackgroundCorrection", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Computing probe affinities (independent of data)");

  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);

  # Alternative #1: Using ACS annotation file
  affinities <- NULL;
  tryCatch({
    affinities <- computeAffinitiesByACS(cdf, ..., verbose=less(verbose));
  }, error = function(ex) {});

  if (is.null(affinities)) {
    # Alternative #2: Using Affymetrix probe-tab files (deprecated)
    affinities <- computeAffinities(cdf, ..., verbose=less(verbose));
  }

  verbose && printf(verbose, "RAM: %.2fMB\n", object.size(affinities)/1024^2);

  verbose && exit(verbose);

  affinities;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod process
#
# @title "Performs background correction"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the output data set.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "GcRmaBackgroundCorrection", function(this, ..., force=FALSE, verbose=FALSE) {
  requireNamespace("gcrma") || throw("Package not loaded: gcrma")
  bg.adjust.fullmodel <- gcrma::bg.adjust.fullmodel
  bg.adjust.affinities <- gcrma::bg.adjust.affinities

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Background correcting data set")

  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already background corrected")
    outputDataSet <- getOutputDataSet(this)
    verbose && exit(verbose)
    return(outputDataSet)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this)

  # Get the output path
  outputPath <- getPath(this)

  # Get algorithm parameters
  params <- getParameters(this)
  opticalAdjust <- params$opticalAdjust
  gsbAdjust <- params$gsbAdjust
  type <- params$type
  indicesNegativeControl <- params$indicesNegativeControl
  seed <- params$seed

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Optical background correction?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (opticalAdjust) {
    obg <- OpticalBackgroundCorrection(ds)
    dsOBG <- process(obg, ..., verbose=verbose)
    ds <- dsOBG
    # Not needed anymore
    obg <- dsOBG <- NULL
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get/calculate affinities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  affinities <- params$affinities
  if (is.null(affinities)) {
    affinities <- calculateAffinities(this, verbose=less(verbose))
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # estimate specific binding (GSB, in gcrma terminology)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  parametersGsb <- NULL
  if (gsbAdjust) {
    verbose && enter(verbose, "Estimating specific binding parameters (data dependent)")
    parametersGsb <- calculateParametersGsb(ds, affinities=affinities, seed=seed, ..., verbose=verbose)
    verbose && cat(verbose, "parametersGsb:")
    verbose && print(verbose, parametersGsb)
    verbose && exit(verbose)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # NSB correction for each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- length(ds)
  verbose && enter(verbose, sprintf("Adjusting %d arrays", nbrOfArrays))

  # Background correction parameters
  fast <- TRUE
  if (fast) {
    k <- 6
    stretch <- 1.15
  } else {
    k <- 0.5
    stretch <- 1
  }
  rho <- 0.7

  cdf <- getCdf(ds)
  chipType <- getChipType(cdf)
  pmCells <- mmCells <- apm <- amm <- anc <- NULL

  dataFiles <- listenv()
  for (ii in seq_along(ds)) {
    df <- ds[[ii]]
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", ii, getName(df), nbrOfArrays))

    filename <- basename(getPathname(df))
    filename <- gsub("[.]cel$", ".CEL", filename)  # Only output upper case!
    pathname <- Arguments$getWritablePathname(filename, path=outputPath,
                                              mustNotExist=FALSE)

    # Already processed?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Already processed. Skipping.")
      # Assert valid file
      dfOut <- newInstance(df, pathname)
      setCdf(dfOut, cdf)
      dataFiles[[ii]] <- pathname
      verbose && exit(verbose)
      next
    }


    if (is.null(mmCells)) {
      verbose && enter(verbose, "Identifying PM and MM cells")
      key <- list(method="bgAdjustGcrma", class=class(df)[1], chipType=chipType, source="gcrma")
      dirs <- c("aroma.affymetrix", chipType)
      indices <- loadCache(key=key, dirs=dirs)
      if (is.null(indices)) {
        indices <- getCellIndices(cdf, useNames=FALSE, unlist=TRUE)
        saveCache(indices, key=key, dirs=dirs)
      }

      # Identify PM & MM cell indices
      # Ordered according to CEL file [whilst isPm() is ordered as the CDF]
      pmCells <- indices[isPm(cdf)]
      mmCells <- indices[!isPm(cdf)]
      verbose && cat(verbose, "Number of PM probes: ", length(pmCells))
      verbose && cat(verbose, "Number of MM probes: ", length(mmCells))
      indices <- NULL ## Not needed anymore
      verbose && exit(verbose)
    }


    verbose && enter(verbose, "Retrieving PM and MM affinities")
    apm <- affinities[pmCells]
    amm <- affinities[mmCells]
    verbose && cat(verbose, "PM affinities:")
    verbose && str(verbose, apm)
    verbose && cat(verbose, "MM affinities:")
    verbose && str(verbose, amm)
    if (!is.null(indicesNegativeControl)) {
      anc <- affinities[indicesNegativeControl]
      verbose && cat(verbose, "Negative-control affinities:")
      verbose && str(verbose, anc)
    }
    verbose && exit(verbose)


    dataFiles[[ii]] %<=% {
      verbose && enter(verbose, "Retrieving PM and MM signals")
      # PM & MM signals
      pm <- getData(df, indices=pmCells)$intensities
      verbose && cat(verbose, "PM signals:")
      verbose && str(verbose, pm)
      mm <- getData(df, indices=mmCells)$intensities
      verbose && cat(verbose, "MM signals:")
      verbose && str(verbose, mm)
      verbose && exit(verbose)

      if (!is.null(indicesNegativeControl)) {
        verbose && enter(verbose, "Retrieving negative-control signals")
        ncs <- getData(df, indices=indicesNegativeControl)$intensities
        verbose && cat(verbose, "Negative-control signals:")
        verbose && str(verbose, ncs)
        stopifnot(length(ncs) == length(anc))
        verbose && exit(verbose)
      } else {
        ncs <- NULL
      }


      # adjust background - use original GCRMA functions to avoid errors from
      # re-coding
      if (type == "fullmodel") {
        verbose && enter(verbose, "Full GCRMA model background adjustment")

        verbose && cat(verbose, "Number of PMs: ", length(pm))
        verbose && cat(verbose, "Number of MMs: ", length(mm))
        pm <- bg.adjust.fullmodel(pms=pm, mms=mm, ncs=ncs, apm=apm, amm=amm, anc=anc, index.affinities=seq_len(length(pm)), k=k, rho=rho, fast=fast)

        verbose && exit(verbose)
      } else if (type == "affinities") {
        verbose && enter(verbose, "Affinity-based GCRMA model background adjustment")

        if (is.null(indicesNegativeControl)) {
          verbose && cat(verbose, "Using mismatch probes (MMs) as negative controls")
          ncs <- mm
          anc <- amm
          nomm <- FALSE  # However...
          # AD HOC: If there are equal number of PMs and MMs, then we assume they
          # are matched, and we will allow gcrma::bg.adjust.affinities() to subset
          # also MMs using 'index.affinities', otherwise not.  See its code:
          #  if (!nomm)
          #    parameters <- bg.parameters.ns(ncs[index.affinities], anc, apm)
          #  else
          #    parameters <- bg.parameters.ns(ncs, anc, apm)
          # However, since we use 'index.affinities' to use all PMs, this special
          # case would equal using all MMs as well, and in that case we can equally
          # well use nomm=TRUE (see its code).
          # /HB 2010-10-02
          nomm <- TRUE
        } else {
          verbose && cat(verbose, "Using a specified set of negative controls")
          nomm <- TRUE
        }

        verbose && enter(verbose, "Dropping perfect-match probes (PMs) with missing signals or missing affinities")

        n0 <- length(pm)
        keepA <- (!is.na(pm))
        n <- sum(keepA)
        verbose && printf(verbose, "Number of finite PMs: %d out of %d (%.1f%%)\n", n, n0, 100*n/n0)

        keepB <- (!is.na(apm))
        n <- sum(keepB)
        verbose && printf(verbose, "Number of finite PM affinities: %d out of %d (%.1f%%)\n", n, n0, 100*n/n0)

        keep <- which(keepA & keepB)
        n <- length(keep)
        verbose && printf(verbose, "Number of PMs kept: %d out of %d (%.1f%%)\n", n, n0, 100*n/n0)

        pm <- pm[keep]
        pmCells <- pmCells[keep]
        apm <- apm[keep]
        # Not needed anymore
        keep <- NULL
        verbose && exit(verbose)


        verbose && enter(verbose, "Dropping negative controls with missing signals or missing affinities")
        n0 <- length(ncs)
        keepA <- (!is.na(ncs))
        n <- sum(keepA)
        verbose && printf(verbose, "Number of finite negative controls: %d out of %d (%.1f%%)\n", n, n0, 100*n/n0)

        keepB <- (!is.na(anc))
        n <- sum(keepB)
        verbose && printf(verbose, "Number of finite negative-control affinities: %d out of %d (%.1f%%)\n", n, n0, 100*n/n0)

        keep <- which(keepA & keepB)
        n <- length(keep)
        verbose && printf(verbose, "Number of negative controls kept: %d out of %d (%.1f%%)\n", n, n0, 100*n/n0)

        ncs <- ncs[keep]
        anc <- anc[keep]
        # Not needed anymore
        keep <- NULL
        verbose && exit(verbose)


        verbose && cat(verbose, "Number of PMs: ", length(pm))
        verbose && cat(verbose, "Number of negative controls: ", length(ncs))

        minNbrOfNegControls <- 1L
        if (length(ncs) < minNbrOfNegControls) {
          throw(sprintf("Cannot perform GCRMA background (type=\"affinities\") correction: The number (%d) of negative control is too small.", length(ncs)))
        }

        pm <- bg.adjust.affinities(pms=pm, ncs=ncs, apm=apm, anc=anc, index.affinities=seq_len(length(pm)), k=k, fast=fast, nomm=nomm)

        verbose && exit(verbose)
      } # if (type == ...)

      mm <- NULL  ## Not needed anymore

      # if specific binding correction requested, carry it out
      if (gsbAdjust && !is.null(parametersGsb)) {
        verbose && enter(verbose, "Adjusting for specific binding")
             #> GSB.adj
             #function (Yin, subset, aff, fit1, k = k)
             #{
             #    y0 = Yin[subset]
             #    y0.adj = k + 2^(-fit1[2] * (aff - mean(aff))) * (y0 - k)
             #    Yin[subset] = y0.adj
             #    Yin
             #}
        #pm <- 2^(log2(pm) - parametersGsb[2]*apm + mean(parametersGsb[2]*apm))  # this is what it used to be
        pm <- k + 2^(-parametersGsb[2] * (apm - mean(apm, na.rm=TRUE))) * (pm - k)
        verbose && exit(verbose)
      }

      apm <- amm <- NULL ## Not needed anymore

      # don't understand this, but it was in original bg.adjust.gcrma(), so
      # we will keep it. /KS
      if (stretch != 1) {
        verbose && enter(verbose, "Stretching")
        verbose && cat(verbose, "Stretch factor: ", stretch)
        mu <- mean(log(pm), na.rm=TRUE)
        pm <- exp(mu + stretch * (log(pm) - mu))
        verbose && exit(verbose)
      }


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Store the adjusted PM signals
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Writing adjusted probe signals")

      # Write to a temporary file (allow rename of existing one if forced)
      isFile <- (force && isFile(pathname))
      pathnameT <- pushTemporaryFile(pathname, isFile=isFile)

      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing")
      createFrom(df, filename=pathnameT, path=NULL, verbose=less(verbose))
      verbose && exit(verbose)

      verbose && enter(verbose, "Writing adjusted intensities")
      verbose && cat(verbose, "Number of cells (PMs only): ", length(pmCells))
      .updateCel(pathnameT, indices=pmCells, intensities=pm)
      verbose && exit(verbose)

      # Rename temporary file
      popTemporaryFile(pathnameT, verbose=verbose)

      ## Create checksum file
      dfZ <- getChecksumFile(pathname)

      verbose && exit(verbose)


      # Assert correctness
      dfOut <- newInstance(df, pathname)
      setCdf(dfOut, cdf)

      pathname
    } ## %<=%


    verbose && exit(verbose)
  } # for (ii ...)

  verbose && exit(verbose)

  ## Resolve futures
  dataFiles <- as.list(dataFiles)
  dataFiles <- NULL

  ## Garbage collect
  gc <- gc()
  verbose && print(verbose, gc)

  verbose && exit(verbose)

  # Gets the output data set
  outputDataSet <- getOutputDataSet(this)

  outputDataSet
})



############################################################################
# HISTORY:
# 2012-11-20
# o Now process() calculates affinities.
# o Added calculateAffinities() to GcRmaBackgroundCorrection.
# 2010-10-01
# o Now GcRmaBackgroundCorrection tries to calculate probe affinites based
#   on ACS annotation files and then as a backup/backward compatibility
#   it uses Affymetrix probe-tab files.
# 2010-09-26
# o Added explicit descriptions to the arguments list of the Rdocs.
# o ROBUSTNESS: Added more validation of the arguments passed to
#   the GcRmaBackgroundCorrection constructor.
# 2007-08-24
# o BUG FIX: Forgot to pass argument '.deprecated=FALSE' to bgAdjustGcrma()
#   because the latter is deprecated at the user-API level.
# 2007-03-21
# o Created.
############################################################################
