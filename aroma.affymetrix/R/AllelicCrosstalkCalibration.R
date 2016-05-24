###########################################################################/**
# @RdocClass AllelicCrosstalkCalibration
#
# @title "The AllelicCrosstalkCalibration class"
#
# \description{
#  @classhierarchy
#
#  This class represents a calibration function that transforms the
#  probe-level signals such that the signals from the two alleles are
#  orthogonal.
#  The method fits and calibrates PM signals only.  MM signals will not
#  affect the model fitting and are unaffected.
# }
#
# @synopsis
#
# \arguments{
#   \item{dataSet}{An @see "AffymetrixCelSet".}
#   \item{...}{Arguments passed to the constructor of
#     @see "ProbeLevelTransform".}
#   \item{model}{A @character string for quickly specifying default
#     parameter settings.}
#   \item{rescaleBy}{A @character string specifying what sets of cells
#     should be rescaled towards a target average, if any.
#     Default is to rescale all cells together.
#     If \code{"none"}, no rescaling is done.}
#   \item{targetAvg}{The signal(s) that either the average of the sum
#     (if one target value) or the average of each of the alleles
#     (if two target values) should have after calibration.
#     Only used if \code{rescaleBy != "none"}.}
#   \item{subsetToAvg}{The indices of the cells (taken as the intersect of
#     existing indices) used to calculate average in order to rescale to
#     the target average. If @NULL, all probes are considered.}
#   \item{mergeShifts}{If @TRUE, the shift of the probe sequence
#     relative to the SNP position is ignored, otherwise not.}
#   \item{B}{An @integer specifying by how many nucleotides the allelic
#     groups should be stratified by. If zero, all SNPs are put in one
#     group.}
#   \item{flavor}{A @character string specifying what algorithm is used
#     to fit the crosstalk calibration.}
#   \item{alpha, q, Q, lambda}{Model fitting parameters.}
#   \item{pairBy}{A @character string specifying how allele probe pairs
#     are identified.}
# }
#
# \section{What probe signals are updated?}{
#   Calibration for crosstalk between allele signals applies by definition
#   only SNP units.
#   Furthermore, it is only SNP units with two or four unit groups that
#   are calibrated.  For instance, in at least on custom SNP CDFs we
#   know of, there is a small number of SNP units that have six groups.
#   \emph{Currently these units are not calibrated (at all).}
#   It is only PM probes that will be calibrated.
#   Note that, non-calibrated signals will be saved in the output files.
# }
#
# \section{What probe signals are used to fit model?}{
#   All PM probe pairs are used to fit the crosstalk model.
#   In the second step where signals are rescaled to a target average,
#   it is possible to specify the set of cells that should be included
#   when estimating the target average.
# }
#
# \section{Important about rescaling towards target average}{
#   Rescaling each allele-pair group (e.g. AC, AG, AT, CG, CT, GC)
#   towards a target average (\code{rescaleBy="groups"})
#   \emph{must not} be used for multi-enzyme chip types,
#   e.g. GenomeWideSNP\_6.
#   If still done, due to confounded effects of non-perfect enzyme
#   mixtures etc, there will be a significant bias between raw CNs
#   for SNPs and CN probes.
#   Instead, for such chip types \emph{all probe signales} should be
#   rescale together towards the target average (\code{rescaleBy="all"}).
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("AllelicCrosstalkCalibration", function(dataSet=NULL, ..., model=c("asis", "auto", "CRMA", "CRMAv2"), rescaleBy=c("auto", "groups", "all", "none"), targetAvg=c(2200, 2200), subsetToAvg="-XY", mergeShifts=TRUE, B=1, flavor=c("sfit", "expectile"), alpha=c(0.1, 0.075, 0.05, 0.03, 0.01), q=2, Q=98, lambda=2, pairBy=c("CDF", "sequence")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  extraTags <- NULL;

  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "AffymetrixCelSet");

    # Argument 'model':
    model <- match.arg(model);
    if (model == "auto") {
      chipType <- getChipType(cdf);
      if (regexpr("^Mapping[0-9]+K_", chipType) != -1) {
        model <- "CRMA";
      } else if (regexpr("^GenomeWideSNP_", chipType) != -1) {
        model <- "CRMAv2";
      } else {
        model <- "CRMAv2";
      }
    }

    if (model == "CRMA") {
      alpha <- c(0.1, 0.075, 0.05, 0.03, 0.01);
      q <- 2;
      Q <- 98;
      lambda <- 2;
      mergeShifts <- TRUE;
      B <- 1;
      rescaleBy <- "auto";
      pairBy <- "CDF";
    } else if (model == "CRMAv2") {
      alpha <- c(0.1, 0.075, 0.05, 0.03, 0.01, 0.0025, 1e-3, 1e-4);
      q <- 2;
      Q <- 98;
      lambda <- 2;
      mergeShifts <- TRUE;
      B <- 1;
      rescaleBy <- "auto";
      pairBy <- "sequence";
    }

    cdf <- getCdf(dataSet);

    # Argument 'rescaleBy':
    rescaleBy <- match.arg(rescaleBy);
    if (rescaleBy == "auto") {
      # First, hardwired...
      # Extract the CDF
      chipType <- getChipType(cdf);
      if (regexpr("^Mapping[0-9]+K_", chipType) != -1) {
        rescaleBy <- "groups";
      } else if (regexpr("^GenomeWideSNP_", chipType) != -1) {
        rescaleBy <- "all";
      } else if (regexpr("^Cyto(genetics|ScanHD)_Array$", chipType) != -1) {
        rescaleBy <- "all";
      } else if (regexpr("^MOUSEDIVm520650$", chipType) != -1) {
        rescaleBy <- "all";
      } else {
        # Heuristics so that we can work with "future/unknown" chip types.
        types <- getUnitTypes(cdf);
        # 5 == Copy Number
        hasCns <- is.element(5, types);
        # Not needed anymore
        types <- NULL;

        if (hasCns) {
          rescaleBy <- "all";
        } else {
          rescaleBy <- "groups";
        }
      }
    }

    if (rescaleBy == "all") {
      targetAvg <- targetAvg[1];
      extraTags <- c(extraTags, rescaleBy="ra");  # rescale all together
    } else if (rescaleBy == "groups") {
      # For backward compatibility, no extra tag (at the moment) /2007-12-01
      extraTags <- c(extraTags);
    } else if (rescaleBy == "none") {
      extraTags <- c(extraTags, rescaleBy="rn");  # rescale nothing
    }

    # Argument 'targetAvg':
    if (!is.null(targetAvg)) {
      targetAvg <- Arguments$getDoubles(targetAvg, range=c(0, Inf));
      if (!length(targetAvg) %in% 1:2) {
        throw("Argument 'targetAvg' must be of length one or two: ", length(targetAvg));
      }
    }

    # Argument 'subsetToAvg':
    if (is.null(subsetToAvg)) {
    } else if (is.character(subsetToAvg)) {
      if (subsetToAvg %in% c("-X", "-Y", "-XY")) {
      } else {
        throw("Unknown value of argument 'subsetToAvg': ", subsetToAvg);
      }
      extraTags <- c(extraTags, subsetToAvg=subsetToAvg);
    } else {
      subsetToAvg <- Arguments$getIndices(subsetToAvg, max=nbrOfCells(cdf));
      subsetToAvg <- unique(subsetToAvg);
      subsetToAvg <- sort(subsetToAvg);
    }
  }

  # Argument 'mergeShifts':
  mergeShifts <- Arguments$getLogical(mergeShifts);

  # Argument 'B':
  B <- Arguments$getInteger(B, range=c(0,3));
  if (B > 1) {
    if (B %% 2 != 1) {
      throw("Argument 'B' must be zero or and odd integer: ", B);
    }
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'pairBy':
  pairBy <- match.arg(pairBy);

  # Argument 'alpha':
  alpha <- Arguments$getDoubles(alpha, range=c(0,1));

  if (flavor == "sfit") {
    algorithmParameters <- list(alpha=alpha, q=q, Q=Q);
  } else if (flavor == "expectile") {
    alpha <- alpha[length(alpha)];
    algorithmParameters <- list(alpha=alpha, lambda=lambda);
  }


  extend(ProbeLevelTransform(dataSet=dataSet, ...), "AllelicCrosstalkCalibration",
    "cached:.setsOfProbes" = NULL,
    "cached:.subsetToAvgExpanded" = NULL,
    .rescaleBy = rescaleBy,
    .targetAvg = targetAvg,
    .subsetToAvg = subsetToAvg,
    .mergeShifts = mergeShifts,
    .B = B,
    .flavor = flavor,
    .algorithmParameters = algorithmParameters,
    .pairBy = pairBy,
    .extraTags = extraTags
  )
})


setMethodS3("getAsteriskTags", "AllelicCrosstalkCalibration", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # shift tag?
  if (!this$.mergeShifts) {
    tags <- c(tags, "byShift");
  }
  # flavor tag
  flavor <- this$.flavor;
  if (flavor != "sfit") {
    tags <- c(tags, flavor);
  }

  # B tag?
  B <- as.integer(this$.B);
  if (B != 1) {
    tags <- c(tags, sprintf("B=%d", B));
  }

  # Extra tags?
  tags <- c(tags, this$.extraTags);

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)



setMethodS3("getSubsetToAvg", "AllelicCrosstalkCalibration", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  subsetToAvg <- this$.subsetToAvg;

  # Expand?
  if (is.character(subsetToAvg)) {
    if (subsetToAvg %in% c("-X", "-Y", "-XY")) {
      verbose && enter(verbose, "Identify subset of units from genome information");
      verbose && cat(verbose, "subsetToAvg: ", subsetToAvg);

      # Look up in cache
      subset <- this$.subsetToAvgExpanded;
      if (is.null(subset)) {
        dataSet <- getInputDataSet(this);
        cdf <- getCdf(dataSet);

        # Get the genome information (throws an exception if missing)
        gi <- getGenomeInformation(cdf);
        verbose && print(verbose, gi);

        # Identify units to be excluded
        if (subsetToAvg == "-X") {
          subset <- getUnitsOnChromosomes(gi, 23, .checkArgs=FALSE);
        } else if (subsetToAvg == "-Y") {
          subset <- getUnitsOnChromosomes(gi, 24, .checkArgs=FALSE);
        } else if (subsetToAvg == "-XY") {
          subset <- getUnitsOnChromosomes(gi, 23:24, .checkArgs=FALSE);
        }

        verbose && cat(verbose, "Units to exclude: ");
        verbose && str(verbose, subset);

        # Identify the cell indices for these units
        subset <- getCellIndices(cdf, units=subset,
                                 useNames=FALSE, unlist=TRUE);
        verbose && cat(verbose, "Cells to exclude: ");
        verbose && str(verbose, subset);

        # The cells to keep
        subset <- setdiff(1:nbrOfCells(cdf), subset);

        verbose && cat(verbose, "Cells to include: ");
        verbose && str(verbose, subset);

        # Store
        this$.subsetToAvgExpanded <- subset;
      }

      subsetToAvg <- subset;
      # Not needed anymore
      subset <- NULL;

      verbose && exit(verbose);
    }
  }

  subsetToAvg;
}, protected=TRUE);



setMethodS3("getParameters", "AllelicCrosstalkCalibration", function(this, expand=TRUE, ...) {
  # Get parameters from super class
  params <- NextMethod("getParameters", expand=expand);

  params <- c(params, list(
    rescaleBy = this$.rescaleBy,
    targetAvg = this$.targetAvg,
    subsetToAvg = this$.subsetToAvg,
    mergeShifts = this$.mergeShifts,
    B = this$.B,
    flavor = this$.flavor,
    algorithmParameters = this$.algorithmParameters
  ));

  # Expand?
  if (expand) {
    params$subsetToAvg <- getSubsetToAvg(this);
  }

  params;
}, protected=TRUE)




setMethodS3("rescale", "AllelicCrosstalkCalibration", function(this, yAll, params, ...) {
  if (params$rescaleBy == "all") {
    yAll <- rescaleByAll(this, yAll=yAll, params=params, ...);
  } else if (params$rescaleBy == "groups") {
    yAll <- rescaleByGroups(this, yAll=yAll, params=params, ...);
  } else if (params$rescaleBy == "none") {
    # Do nothing.
  }

  yAll;
}, protected=TRUE)


setMethodS3("rescaleByAll", "AllelicCrosstalkCalibration", function(this, yAll, params, ..., verbose=FALSE) {
  targetAvg <- params$targetAvg;
  nt <- length(targetAvg);
  if (nt != 1) {
    throw("In order rescale towards a global target average ('rescaleBy' == \"all\"), argument 'targetAvg' must be a scalar: ", paste(targetAvg, collapse=","));
  }

  if (verbose) {
    enter(verbose, "Rescaling toward target average");
    printf(verbose, "Rescale by: %s\n", params$rescaleBy);
    cat(verbose, "Target average: ", targetAvg);
    if (!is.null(params$subsetToAvg)) {
      cat(verbose, "Using subset of cells for estimate of target average:");
      src <- attr(params$subsetToAvg, "src");
      if (is.null(src))
         src <- params$subsetToAvg;
      str(verbose, src);
    }
    cat(verbose, "yAll: ");
    str(verbose, yAll);
  }

  # Total number of values
  n0 <- length(yAll);

  # Average of *all* values
  yAvg0 <- median(yAll, na.rm=TRUE);

  if (!is.null(params$subsetToAvg)) {
    y <- yAll[params$subsetToAvg];
    n <- length(y);
    if (n == 0) {
      throw("Cannot rescale to target average. There are no cells to average over.");
    }
    yAvg <- median(y, na.rm=TRUE);
    verbose && printf(verbose, "yAvg (using %d/%.1f%% summed pairs): %.2f of %.2f (%.1f%%)\n", n, 100*n/n0, yAvg, yAvg0, 100*yAvg/yAvg0);
    # Not needed anymore
    y <- n <- NULL;
  } else {
    yAvg <- yAvg0;
    verbose && printf(verbose, "yAvg (100%%): %.2f\n", yAvg);
  }

  if (!is.finite(yAvg)) {
    throw("Cannot rescale to target average. Signal average is non-finite: ", yAvg);
  }

  b <- targetAvg/yAvg;
  verbose && printf(verbose, "Scale factor: %.2f\n", b);
  yAll <- b*yAll;
  fit <- list(all=list(b=b));
  verbose && exit(verbose);

  attr(yAll, "fit") <- fit;

  verbose && cat(verbose, "Rescaling parameter estimates:");
  verbose && str(verbose, fit);
  verbose && exit(verbose);

  yAll;
}, protected=TRUE)


setMethodS3("rescaleByGroups", "AllelicCrosstalkCalibration", function(this, yAll, params, setsOfProbes, ..., verbose=FALSE) {
  basepairs <- names(setsOfProbes$snps);
  nbrOfPairs <- length(basepairs);

  # Infer method?
  targetAvg <- params$targetAvg;
  nt <- length(targetAvg);
  if (nt == 1) {
    method <- "sum";
  } else if (nt == 2) {
    method <- "allele";
  } else {
    throw("In order rescale towards target average, argument 'targetAvg' must be either a scalar or a vector of length two: ", paste(targetAvg, collapse=","));
  }

  if (verbose) {
    enter(verbose, "Rescaling toward target average");
    printf(verbose, "Rescale by: %s\n", params$rescaleBy);
    cat(verbose, "Target average(s): ", paste(targetAvg, collapse=", "));
    printf(verbose, "Method: %s\n", method);
    if (!is.null(params$subsetToAvg)) {
      cat(verbose, "Using subset of cells for estimate of target average:");
      src <- attr(params$subsetToAvg, "src");
      if (is.null(src))
         src <- params$subsetToAvg;
      str(verbose, src);
    }
  }


  snps <- vector("list", nbrOfPairs);
  names(snps) <- basepairs;
  subset <- list(snps=snps, nonSNPs=NULL);
  # Not needed anymore
  snps <- NULL;
  fit <- list(
    dimY = dim(yAll),
    params = params,
    method = method,
    subset = subset
  );
  # Not needed anymore
  subset <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Rescale based on y = yA+yB
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (method == "sum") {
    for (kk in seq_len(nbrOfPairs)) {
      name <- basepairs[kk];
      verbose && enter(verbose, sprintf("Allele-pair group #%d ('%s') of %d", kk, name, nbrOfPairs));

      # Get data pairs
      idxAB <- setsOfProbes$snps[[name]];
      idxAB <- matrix(idxAB, ncol=2, byrow=FALSE);

      # Sum y=yA+yB
      y <- yAll[idxAB[,1]]+yAll[idxAB[,2]];
      n <- length(y);
      n0 <- n;

      # Calculate current average
      yAvg <- median(y, na.rm=TRUE);
      yAvg0 <- yAvg;
      # Not needed anymore
      y <- NULL;

      if (!is.null(params$subsetToAvg)) {
        keep <- matrix((idxAB %in% params$subsetToAvg), ncol=2, byrow=FALSE);
        keep <- (keep[,1] & keep[,2]);
        idxAB <- idxAB[keep,,drop=FALSE];
        # Not needed anymore
        keep <- NULL;

        # Sum y=yA+yB
        y <- yAll[idxAB[,1]]+yAll[idxAB[,2]];
        n <- length(y);

        if (n == 0) {
          throw("Cannot rescale to target average. After taking the intersect of the subset of cells to be used, there are no cells left.");
        }

        yAvg <- median(y, na.rm=TRUE);
        # Not needed anymore
        y <- NULL;

        verbose && printf(verbose, "yAvg (using %d/%.1f%% summed pairs): %.2f of %.2f (%.1f%%)\n", n, 100*n/n0, yAvg, yAvg0, 100*yAvg/yAvg0);
      } else {
        verbose && printf(verbose, "yAvg (100%%): %.2f\n", yAvg);
      }

      if (!is.finite(yAvg))
        throw("Cannot rescale to target average. Signal average is non-finite: ", yAvg);

      # Rescale
      b <- targetAvg/yAvg;
      verbose && printf(verbose, "scale factor: %.2f\n", b);

      idxAB <- setsOfProbes$snps[[name]];
      yAll[idxAB] <- b*yAll[idxAB];

      # Not needed anymore
      idx <- NULL;

      fit$subset$snps[[name]] <- list(b=b);

      verbose && exit(verbose);
    } # for (kk in ...)
  } # if (method == "sum")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Rescale based on (yA,yB)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (method == "allele") {
    for (kk in seq_len(nbrOfPairs)) {
      name <- basepairs[kk];
      verbose && enter(verbose, sprintf("Allele-pair group #%d ('%s') of %d", kk, name, nbrOfPairs));

      # Get data pairs
      idxAB <- setsOfProbes$snps[[name]];
      idxAB <- matrix(idxAB, ncol=2, byrow=FALSE);

      # Default scale factors
      b <- c(1,1);

      # For each allele
      for (cc in 1:2) {
        idx <- idxAB[,cc];
        y <- yAll[idx];
        n <- length(y);
        n0 <- n;

        # Calculate current average
        yAvg <- median(y, na.rm=TRUE);
        yAvg0 <- yAvg;
        # Not needed anymore
        y <- NULL;

        if (!is.null(params$subsetToAvg)) {
          idx <- idxAB[,cc];
          idx <- intersect(idx, params$subsetToAvg);
          y <- yAll[idx];
          n <- length(y);

          if (n == 0) {
            throw("Cannot rescale to target average. After taking the intersect of the subset of cells to be used, there are no cells left.");
          }

          yAvg <- median(y, na.rm=TRUE);
          # Not needed anymore
          y <- NULL;

          verbose && printf(verbose, "yAvg (using %d/%.1f%% pairs): %.2f of %.2f (%.1f%%)\n", n, 100*n/n0, yAvg, yAvg0, 100*yAvg/yAvg0);
        } else {
          verbose && printf(verbose, "yAvg (100%%): %.2f\n", yAvg);
        }

        if (!is.finite(yAvg))
          throw("Cannot rescale to target average. Signal average is non-finite: ", yAvg);

        # Rescale
        b[cc] <- targetAvg[cc]/yAvg;
        verbose && printf(verbose, "scale factor: %.2f\n", b[cc]);

        idx <- idxAB[,cc];
        yAll[idx] <- b[cc]*yAll[idx];

        # Not needed anymore
        idx <- NULL;
      } # for (cc ...)

      fit$subset$snps[[name]] <- list(b=b);
      verbose && exit(verbose);
    } # for (kk in ...)
  } # if (method == "allelic")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Rescaling non-SNP cells
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Non-SNP cells");
  idx <- setsOfProbes$nonSNPs;

  if (!is.null(params$subsetToAvg)) {
    idx <- intersect(idx, params$subsetToAvg);
  }

  y <- yAll[idx];
  n <- length(y);
  if (n > 0) {
    yAvg <- median(y, na.rm=TRUE);
    # Not needed anymore
    y <- NULL;
    if (!is.finite(yAvg))
      throw("Cannot rescale to target average. Signal average is non-finite: ", yAvg);

    # Rescale (to half of the allele target averages)
    #  AA : P(2A,0B)=0.25
    #  AB : P(1A,1B)=0.50
    #  BB : P(0A,2B)=0.25
    # => P(A=0)=0.25, P(A=1)=0.5, P(A=2)=0.25, ...
    # => E[A] = 0*0.25 + 1*0.5 + 2*0.25 = 1
    #    E[B] = 0*0.25 + 1*0.5 + 2*0.25 = 1
    # => E[cA] = 1c, E[dB] = 1d
    # => E[cA+cB] = c+d
    # => E[A+B] = E[A]+E[B] = 2
    # => E[c(A+B)] = 2c
    if (method == "allele") {
      #  c = targetAvg[1], d = targetAvg[2]
      #  => (c+d) = sum(targetAvg)
      targetAvg <- sum(params$targetAvg);
    } else if (method == "sum") {
      #  2c = targetAvg
      targetAvg <- params$targetAvg;
    }
    b <- targetAvg/yAvg;
    verbose && printf(verbose, "scale factor: %.2f\n", b);
    yAll[idx] <- b*yAll[idx];
    fit$subset$nonSNPs <- list(b=b);
  } # if (n > 0)
  # Not needed anymore
  idx <- NULL;
  verbose && exit(verbose);

  attr(yAll, "fit") <- fit;

  verbose && cat(verbose, "Rescaling parameter estimates:");
  verbose && str(verbose, fit);
  verbose && exit(verbose);

  yAll;
}, protected=TRUE)




###########################################################################/**
# @RdocMethod process
#
# @title "Calibrates the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already calibrated is re-calibrated,
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "AllelicCrosstalkCalibration", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calibrating data set for allelic cross talk");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already calibrated");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get algorithm parameters
  params <- getParameters(this);

  # Get (and create) the output path
  outputPath <- getPath(this);

  # Model flavor and parameters
  flavor <- params$flavor;
  algorithmParameters <- params$algorithmParameters;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For hooks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hookName <- "process.AllelicCrosstalkCalibration";


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Precalculate some model fit parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Compressing model parameter to a short format");
  paramsShort <- params;
  paramsShort$subsetToAvg <- NULL;
#  paramsShort$subsetToAvgIntervals <- seqToIntervals(params$subsetToAvg);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibrate each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(ds);
  nbrOfArrays <- length(ds);
  verbose && enter(verbose, "Calibrating ", nbrOfArrays, " arrays");
  verbose && cat(verbose, "Path: ", outputPath);

  res <- listenv()

  for (kk in seq_len(nbrOfArrays)) {
    df <- ds[[kk]]
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                              kk, getName(df), nbrOfArrays))

    fullname <- getFullName(df)
    filename <- sprintf("%s.CEL", fullname)
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...)

    # Already calibrated?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Calibrated data file already exists: ", pathname)
      dfC <- newInstance(df, pathname)
      setCdf(dfC, cdf)
      res[[kk]] <- pathname
      verbose && exit(verbose)
      next
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculations used by all samples
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    setsOfProbes <- getSetsOfProbes(this, verbose=less(verbose, 1))
    verbose && cat(verbose, "setsOfProbes:")
    verbose && str(verbose, setsOfProbes)


    res[[kk]] %<=% {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Reading data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Reading all probe intensities")
      yAll <- getData(df, fields="intensities", ...)$intensities
      verbose && exit(verbose)

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Calibrating
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ## callHooks(sprintf("%s.onBegin", hookName), df=df, setsOfProbes=setsOfProbes, ...);

      modelFit <- list(
        paramsShort=paramsShort
      )

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Fitting each allelic basepair
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      basepairs <- names(setsOfProbes$snps)
      nbrOfPairs <- length(basepairs)
      fits <- vector("list", length=nbrOfPairs)
      names(fits) <- basepairs
      verbose && enter(verbose, "Fitting calibration model")

      for (pp in seq_len(nbrOfPairs)) {
        name <- basepairs[pp]
        verbose && enter(verbose, sprintf("Allele probe-pair group #%d ('%s') of %d", pp, name, nbrOfPairs))
        basepair <- unlist(strsplit(name, split=""))
        idx <- setsOfProbes$snps[[name]]

        verbose && enter(verbose, "Fitting")
        y <- matrix(yAll[idx], ncol=2, byrow=FALSE)
        verboseL <- (verbose && isVisible(verbose, -50))
        # Not needed anymore
        idx <- NULL

        verbose && cat(verbose, "Model/algorithm flavor: ", flavor)
        if (flavor == "sfit") {
          alpha <- algorithmParameters$alpha
          q <- algorithmParameters$q
          Q <- algorithmParameters$Q
          verbose & cat(verbose, "Model parameters:")
          verbose & str(verbose, list(alpha=alpha, q=q, Q=Q))
          verbose & cat(verbose, "Number of data points: ", nrow(y))
          if (nrow(y) > 10) {
            fit <- fitGenotypeCone(y, flavor=flavor, alpha=alpha, q=q, Q=Q,
                                                         verbose=verboseL)
          } else {
            fit <- NULL
            verbose & cat(verbose, "Cannot fit model: too few data points. Skipping this group: ", name)
          }
        } else if (flavor == "expectile") {
          alpha <- algorithmParameters$alpha
          lambda <- algorithmParameters$lambda
          verbose & cat(verbose, "Model parameters:")
          verbose & str(verbose, list(alpha=alpha, lambda=lambda))
          verbose & cat(verbose, "Number of data points: ", nrow(y))
          fit <- fitGenotypeCone(y, flavor=flavor, alpha=alpha,
                                        lambda=lambda, verbose=verboseL)
        }
        verbose && print(verbose, fit, level=-5)
        fits[[name]] <- fit
        verbose && exit(verbose)

        callHooks(sprintf("%s.onFitOne", hookName), df=df, y=y, fit=fit, ...)

        # Not needed anymore
        y <- fit <- NULL ## Not needed anymore
        gc <- gc()

        verbose && exit(verbose)
      } # for (pp in seq_len(nbrOfPairs))
      verbose && exit(verbose)


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Estimate offset for non-SNP PM cells
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ns <- sapply(fits, FUN=function(fit) fit$dimData[1])

      # Alt (1) Weighted average of all offset estimates
      w <- ns / sum(ns)
      origins <- sapply(fits, FUN=function(fit) fit$origin)
      verbose && cat(verbose, "Estimated origins:")
      verbose && print(verbose, origins)
      origins <- colMeans(origins, na.rm=TRUE)
      offset <- sum(w*origins, na.rm=TRUE)
      verbose && printf(verbose, "Weighted average offset: %.2f\n", offset)
      # Not needed anymore
      origins <- w <- NULL

      fit <- list(
        offset = offset,
        ns = ns,
        dimData = length(setsOfProbes$nonSNPs)
      )
      fits[["nonSNPs"]] <- fit
      # Not needed anymore
      fit <- ns <- offset <- NULL

      # Store allelic crosstalk model fits
      modelFit$accFits <- fits


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Harmonizing parameter estimates?
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Here we can harmonize the estimates, e.g. make all offset estimates
      # the same.  This might be useful if we want to correct other probes
      # not included above such as CN probes on SNP 6.0. /HB 2007-09-05

      callHooks(sprintf("%s.onFit", hookName), df=df, y=y, basepair=basepair, ...)

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Backtransforming (calibrating)
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Backtransforming (calibrating) data")
      for (pp in seq_len(nbrOfPairs)) {
        name <- basepairs[pp]
        verbose && enter(verbose, sprintf("Allele basepair #%d ('%s') of %d", pp, name, nbrOfPairs))

        idx <- setsOfProbes$snps[[name]]
        y <- matrix(yAll[idx], ncol=2, byrow=FALSE)
        fit <- fits[[name]]
        if (!is.null(fit)) {
          yC <- backtransformGenotypeCone(y, fit=fit)
        } else {
          verbose && cat(verbose, "Cannot do backtransformation because there were to few data points in group to fit anything: ", name)
          yC <- y
        }
        yAll[idx] <- yC

        ## callHooks(sprintf("%s.onUpdated", hookName), df=df, y=y, basepair=basepair, fit=fits[[name]], yC=yC,...)
        # Not needed anymore
        idx <- y <- yC <- NULL
        gc <- gc()

        verbose && exit(verbose)
      } # for (pp in seq_len(nbrOfPairs))
      verbose && exit(verbose)


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Correcting offset for all non-SNP cells
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      cells <- setsOfProbes$nonSNPs
      if (length(cells) > 0) {
        verbose && enter(verbose, "Correcting offset for all non-SNP cells")
        verbose && cat(verbose, "Cells:")
        verbose && str(verbose, cells)
        offset <- fits[["nonSNPs"]]$offset
        verbose && cat(verbose, "Offset: ", offset)
        yAll[cells] <- yAll[cells] - offset
        verbose && exit(verbose)
      }
      # Not needed anymore
      cells <- fits <- NULL

      # Garbage collect
      gc <- gc()


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Rescaling toward target average?
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!is.null(params$targetAvg)) {
        yAll <- rescale(this, yAll=yAll, params=params,
                           setsOfProbes=setsOfProbes, verbose=less(verbose))
        fit <- attr(yAll, "fit")
        fit$params <- NULL
        fit$paramsShort <- paramsShort
        modelFit$rescaleFit <- fit
        # Not needed anymore
        fit <- NULL

        # Garbage collect
        gc <- gc()
        verbose && print(verbose, gc)
      }



      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Store model fit
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Store fit and parameters (in case someone are interested in looking
      # at them later; no promises of backward compatibility though).
      filename <- sprintf("%s,fit.RData", fullname)
      fitPathname <- Arguments$getWritablePathname(filename, path=outputPath, ...)

      saveObject(modelFit, file=fitPathname)
      verbose && str(verbose, modelFit, level=-50)
      # Not needed anymore
      modelFit <- NULL

        # Garbage collect
      gc <- gc()
      verbose && print(verbose, gc)


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Storing data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Storing calibrated data")

      # Write to a temporary file (allow rename of existing one if forced)
      isFile <- (force && isFile(pathname))
      pathnameT <- pushTemporaryFile(pathname, isFile=isFile, verbose=verbose)

      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing")
      createFrom(df, filename=pathnameT, path=NULL, verbose=less(verbose))
      verbose && exit(verbose)

      # Write calibrated data to file
      verbose2 <- -as.integer(verbose)-2
      .updateCel(pathnameT, intensities=yAll, verbose=verbose2)

      # Not needed anymore
      yAll <- verbose2 <- NULL

      # Rename temporary file
      popTemporaryFile(pathnameT, verbose=verbose)

      gc <- gc()
      verbose && print(verbose, gc)
      verbose && exit(verbose)

      # Assert validity of the calibrated data file
      dfC <- newInstance(df, pathname)
      setCdf(dfC, cdf)

      ## Create checksum file
      dfCZ <- getChecksumFile(dfC)

      ## callHooks(sprintf("%s.onExit", hookName), df=df, dfC=dfC, ...);

      dfC
    } ## %<=%

    verbose && exit(verbose);
  } # for (kk in seq_len(nbrOfArrays))
  verbose && exit(verbose);

  ## Not needed anymore
  ds <- setsOfProbes <- NULL

  ## Resolve futures
  res <- as.list(res)
  res <- NULL

  ## Garbage collect
#  clearCache(this);
  gc <- gc();
  verbose && print(verbose, gc);


  outputDataSet <- getOutputDataSet(this, force=TRUE);

  verbose && exit(verbose);

  invisible(outputDataSet);
})




setMethodS3("plotBasepair", "AllelicCrosstalkCalibration", function(this, array, basepairs=NULL, what=c("before", "after"), ..., plotFcn=NULL, xlim=c(-500,65535), ylim=xlim, linesFcn=NULL, lwd=4, lcol="red", scale=1, force=FALSE, verbose=FALSE) {
  linesAllelicCrosstalk <- function(a, B, max=1e5, ...) {
    bA <- B[1,2]/B[1,1];
    bB <- B[2,1]/B[2,2];
    lines(x=a[1]+c(0,1)*max, y=a[2]+c(0,1)*max*bA, ...);
    lines(x=a[1]+c(0,1)*max*bB, y=a[2]+c(0,1)*max, ...);
  } # linesAllelicCrosstalk()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'what':
  what <- match.arg(what);

  if (what == "before") {
    cs <- getInputDataSet(this);
  } else if (what == "after") {
    cs <- getOutputDataSet(this);
  }
  if (is.null(cs)) {
    throw("Requested data set does not exist: ", what);
  }

  # Argument 'array':
  array <- Arguments$getIndex(array, max=length(cs));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the cell indices for each possible allele basepair.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying PM cell indices for each possible allele basepair");
  setsOfProbes <- getSetsOfProbes(this, verbose=less(verbose, 1));
  verbose && exit(verbose);

  # Argument 'basepair':
  knownBasepairs <- names(setsOfProbes$snps);
  if (is.null(basepairs)) {
    basepairs <- knownBasepairs;
  } else {
    if (!all(basepairs %in% knownBasepairs)) {
      throw("Argument 'basepairs' refers one or several unknown basepairs: ", paste(basepairs, collapse=", "));
    }
  }

  # Argument 'plotFcn':
  if (is.null(plotFcn)) {
    plotFcn <- function(..., pch=NA, transformation=function(x) x^0.33) {
      smoothScatter(..., pch=pch, transformation=transformation);
    }
  } else if (!is.function(plotFcn)) {
    throw("Argument 'plotFcn' is not a function: ", mode(plotFcn));
  }


  # Argument 'linesFcn':
  if (is.null(linesFcn)) {
    linesFcn <- linesAllelicCrosstalk;
  } else if (!is.function(linesFcn)) {
    throw("Argument 'linesFcn' is not a function: ", mode(linesFcn));
  }

  # Get the data file
  cf <- cs[[array]];
  verbose && print(verbose, cf);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load model parameter estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (what == "before") {
    fullname <- getFullName(cf);
    filename <- sprintf("%s,fit.RData", fullname);
    pathname <- Arguments$getWritablePathname(filename, path=getPath(this), ...);
    if (!isFile(pathname)) {
      throw("File containing parameter estimates not found: ", pathname);
    }

    modelFit <- loadObject(pathname);
    missing <- basepairs[!(basepairs %in% names(modelFit$accFits))];
    if (length(missing) > 0) {
      throw("Loaded model fit does not contain estimates for some basepairs: ",
                                                paste(missing, collapse=", "));
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading all probe intensities");
  yAll <- getData(cf, fields="intensities", ...)$intensities;
  yAll <- scale * yAll;
  verbose && summary(verbose, yAll);
  verbose && exit(verbose);

  hookName <- "process.AllelicCrosstalkCalibration";

  nbrOfPairs <- length(basepairs);
  verbose && enter(verbose, "Plotting (PMA,PMB) for ", nbrOfPairs, " basepair(s) ", what, " calibration");
  for (kk in seq_along(basepairs)) {
    name <- basepairs[kk];
    verbose && enter(verbose, "Plotting (PMA,PMB) for basepair ", name, " in array ", array, " ('", getName(cf), "')");

    alleles <- strsplit(name, split="")[[1]];

    # Extracting data
    idx <- setsOfProbes$snps[[name]];
    y <- matrix(yAll[idx], ncol=2, byrow=FALSE);
    colnames(y) <- alleles;

    # Exclude data points outside plot area (faster)
    y[y < xlim[1]-0.1*diff(xlim) | y > xlim[2] + 0.1*diff(xlim)] <- NA;
    keep <- (!is.na(y[,1]) & !is.na(y[,2]));
    y <- y[keep,,drop=FALSE];
    # Not needed anymore
    keep <- NULL;

    plotFcn(y, xlim=xlim, ylim=ylim, ...);

    if (!is.null(linesFcn)) {
      if (what == "before") {
        fit <- modelFit$accFits[[name]];
        verbose && cat(verbose, "Parameter estimates:", level=-20);
        verbose && str(verbose, fit, level=-20);
      } else if (what == "after") {
        fit <- list(origin=c(0,0), W=matrix(c(1,0,0,1), ncol=2, byrow=FALSE));
      }

      linesFcn(a=fit$origin, B=fit$W, lwd=lwd, col=lcol, ...);
    }

    # Not needed anymore
    y <- idx <- fit <- NULL; # Not needed anymore
    gc <- gc();

    verbose && exit(verbose);
  } # for (kk in ...)
  verbose && exit(verbose);

}, protected=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NOT USED
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getDataPairs", "AllelicCrosstalkCalibration", function(this, array, cs=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (is.null(cs)) {
    cs <- getInputDataSet(this);
  }

  verbose && enter(verbose, "Identifying cell indices for each possible allele basepair");
  setsOfProbes <- getSetsOfProbes(this, verbose=less(verbose, 1));
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading all probe intensities");
  cf <- cs[[array]];
  yAll <- getData(cf, fields="intensities", ...)$intensities;
  verbose && exit(verbose);

  nbrOfPairs <- length(setsOfProbes$snps);
  res <- vector("list", nbrOfPairs);
  names(res) <- names(setsOfProbes$snps);

  verbose && enter(verbose, "Extracting data pairs");
  for (kk in seq_len(nbrOfPairs)) {
    name <- names(setsOfProbes$snps)[kk];
    basepair <- unlist(strsplit(name, split=""));
    idx <- setsOfProbes$snps[[name]];
    y <- matrix(yAll[idx], ncol=2, byrow=FALSE);
    colnames(y) <- c("A", "B");
    res[[kk]] <- y;
    # Not needed anymore
    y <- NULL;
  }
  verbose && exit(verbose);

  res;
}, protected=TRUE);



############################################################################
# HISTORY:
# 2011-11-19
# o SPEEDUP: Constructor for AllelicCrosstalkCalibration no longer calls
#   the very slow hasUnitTypes() for AffymetrixCdfFile, but instead only
#   the much faster getUnitTypes().
# 2011-03-14
# o SPEEDUP: Added MOUSEDIVm520650 to the set of predefined chip types
#   in AllelicCrosstalkCalibration.
# 2010-05-10
# o Now the constructor AllelicCrosstalkCalibration() is set to recognize
#   the Cytogenetics_Array chip type.  This avoids having to scan the CDF
#   for unit types and check for SNPs, which is slow and not really wanted
#   for a constructor function.  The allele pairs are inferred from the
#   CDF (by default).
# 2010-02-15
# o MEMORY OPTIMIZATION: Now process() of AllelicCrosstalkCalibration
#   clears the in-memory cache when finished.
# 2009-09-04
# o Now smoothScatter() is loaded via aroma.core.
# 2008-12-17
# o process() had a nested for loop with the same iteration variable 'kk'
#   was used for in the inner and the outer loop.  I was surprised to find
#   that this worked and is valid in R, cf. help(for).  However, for
#   readability I've changed it to use to 'kk' and 'pp'.
# 2008-12-10
# o BUG FIX: Now process() avoids sets of pairs with too few probe pairs.
#   This could happen because of the new getSetsOfProbes() working off
#   the probe sequence files.
# 2008-11-28
# o Added argument 'model' for quick specification of default parameter
#   settings according to the CRMA or CRMA v2 model.
# o Added argument 'pairBy' to specify how the allele probe pairs are
#   identified.  Historically we inferred this from the CDF, but it is
#   safer and more generic to do this from the probe sequences, which then
#   requires an ACS cell-sequence annotation file.
# 2008-08-31
# o BUG FIX: The allele pairs identified was not correct for GWS arrays.
# o Updated AllelicCrosstalkCalibration to support flavor 'expectile' too.
# 2008-08-30
# o Added argument 'mergeShifts=TRUE' and 'B=1'.  Currently B=0 and B=1
#   is supported.
# 2008-08-29
# o Added protected getSetsOfProbes().  By overriding this method, other
#   sets of probes can be used.  This function might later also recognize
#   user specified function generating the sets.
# 2008-08-04
# o Added support to fit the genotype "cone" using the 'expectile' package
#   instead of the 'sfit' package. This is controlled by the 'flavor'
#   argument of the constructor.
# 2008-07-14
# o Now explicitly using matrix(..., byrow=FALSE).
# 2008-05-30
# o BUG FIX: The constructor of AllelicCrosstalkCalibration used
#   non-defined variable 'verbose'.
# 2008-02-21
# o Now SNPs and CN probes are infered from getUnitTypes(cdf) and no longer
#   from the unit names.
# 2008-02-14
# o Now 'verbose' is passed as a logical argument to fitGenotypeCone().
# 2007-12-01
# o MEMORY OPTIMIZATION: Added clearCache() to AllelicCrosstalkCalibration.
# o BUG FIX: The AllelicCrosstalkCalibration introduced in previous version
#   was broken for 10K (maybe 100K and 500K as well).
# o Now 'subsetToAvg' of AllelicCrosstalkCalibration accepts '-XY' (and
#   '-X' and '-Y') for automatic look up of all units and exclude those
#   that are on ChrX and ChrY.  Note, '-XY' will work on all chip types,
#   also older ones for which there are no ChrY units.
# o The new constructor argument 'rescaleBy' now sets a "subtag", e.g.
#   'ACC,ra' where the 'ra' indicates that 'rescaleBy=all' was used.
# 2007-11-27
# o Added specific getAsteriskTag() for AllelicCrosstalkCalibration.
# o Added the option to rescale towards a target average of *all* probes.
# o Now rescale() throws an error if length(params$targetAvg) is not 1 or 2.
#   This should never happend, but was added just in case.
# 2007-09-14
# o BUG FIX: Now the target average of non-SNP probes is half of the target
#   average of alleles.
# 2007-09-09
# o Added alpha version of plotBasepair() to AllelicCrosstalkCalibration.
# 2007-09-08
# o Now AllelicCrosstalkCalibration corrects also non-SNP PM cells by
#   substracting a global offset and rescaling towards target average.
#   The global offset is calculated as the weighted average of all
#   allelic offsets.  This is the simplest way to incorporate a calibration
#   for non-SNP cells.  A more advanced version would be to stratify by
#   middle nucleotide, and use the corresponding estimates from the SNP
#   cells, but that would require information about the middle nucleotide
#   for all cells, which got some overhead.  Hopefully the simpler version
#   is good enough.
# 2007-09-05
# o Now the rescaling can be done either on (yA,yB) separately or on
#   y=yA+yB.  If targetAvg has two values the former, otherwise the latter.
# o Now AllelicCrosstalkCalibration recognizes argument 'subsetToAvg'.
# o Now process() stores the crosstalk settings and estimated parameters
#   to file. May be useful if one wants to go back and look at the details.
#   One day we might get around to store this information in the CEL file
#   headers.
# o Now process() first fits the crosstalk model for all basepairs, then
#   backtransform the signals, then optional rescale signals to target
#   average, then saves the calibrated signals.
# o SPEED UP: Now getAlleleProbePairs() is only called if data needs to be
#   calibrated, i.e. if already calibrated it is not loaded.
# o CLEAN UP: Now the code of calibrateAllelicCrosstalk() is included here.
# 2007-03-29
# o Now 'targetAvg' defaults to 2200 so that allele A and allele B signals
#   are rescaled to be one the same scale.  If so,  it does not make sense
#   to do background correction afterwards.
# o Added getParameters().
# o Added support for arguments 'targetAvg', 'alpha', 'q', and 'Q'.
# 2006-12-08
# o Now this class inherits from the ProbePreprocessing class.
# o Now this pre-processor output results to probeData/.
# o Renamed from AllelicCrosstalkCalibrator.
# 2006-11-18
# o Removed version and subversion tags, and related functions.
#   Now getTags() returns the tags of the input data set plus any tags
#   of this instance.
# 2006-11-02
# o Created from QuantileNormalizer.R.
############################################################################
