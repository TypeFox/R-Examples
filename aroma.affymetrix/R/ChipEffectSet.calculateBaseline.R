###########################################################################/**
# @set "class=ChipEffectSet"
# @RdocMethod calculateBaseline
#
# @title "Estimates the baseline signal chromosome by chromosome"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{chromosomes}{An @integer @vector specifying for which chromsosomes
#     the baseline should be estimated.
#     If @NULL, all chromosomes are considered.}
#   \item{ploidy}{An @integer specifying the ploidy that the baseline
#     should have.}
#   \item{defaultPloidy}{An @integer specifying the default ploidy of
#     chromosomes that have not explicitly been allocated one.}
#   \item{all}{If @TRUE, signals are averaged also for cells that are not
#     on the genome.}
#   \item{force}{If @TRUE, the CEL file that stores the is recreated.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{...}{Not used.}
# }
#
# @author "HB"
#
# \seealso{
#   @see "getAverageFile".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("calculateBaseline", "ChipEffectSet", function(this, chromosomes=NULL, ploidy=2, defaultPloidy=NA, all=FALSE, force=FALSE, verbose=FALSE, ...) {
  cdf <- getCdf(this);
  gi <- getGenomeInformation(cdf);
  allChromosomes <- getChromosomes(gi);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chromosomes':
  if (is.null(chromosomes)) {
    chromosomes <- allChromosomes;
  } else {
    chromosomes <- Arguments$getChromosomes(chromosomes,
                                                range=range(allChromosomes));
  }

  # Argument 'ploidy':
  ploidy <- Arguments$getInteger(ploidy, range=c(1,8));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Estimating the baseline signals for each chromosome");


  verbose && enter(verbose, "Getting CEL file to store baseline signals");
  csBaseline <- getBaseline(this, force=force, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && enter(verbose, "Checking for non-estimated cells");
  ds <- getData(csBaseline, fields="intensities", verbose=less(verbose))$intensities;
  ncells <- length(ds);
  todo <- which(isZero(ds));
  ntodo <- length(todo);
  # Not needed anymore
  ds <- NULL;

  verbose && printf(verbose, "Found %d (%.1f%%) non-estimated cells.\n",
                                                  ntodo, 100*ntodo/ncells);
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  n <- length(this);
  for (chromosome in chromosomes) {
    verbose && enter(verbose, "Chromosome ", chromosome);

    verbose && enter(verbose, "Identifying units on chromosome");
    units <- getUnitsOnChromosome(gi, chromosomes=chromosome);
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);
    verbose && exit(verbose);

    verbose && enter(verbose, "Identifying cells for these units");
    cells <- getCellIndices(this, units=units);
    # Not needed anymore
    units <- NULL;
    cells <- unlist(cells, use.names=FALSE);
    cells <- sort(cells);
    ncells <- length(cells);
    verbose && cat(verbose, "Cells:");
    verbose && str(verbose, cells);
    verbose && exit(verbose);

    if (!force) {
      verbose && enter(verbose, "Checking for non-estimated loci");
      cells <- intersect(cells, todo);
      nkeep <- length(cells);

      verbose && printf(verbose, "Found %d (%.1f%%) non-estimated loci.\n",
                                              nkeep, 100*nkeep/ncells);
      verbose && exit(verbose);
      if (nkeep == 0) {
        verbose && cat(verbose, "Baseline averages already exist for all loci on this chromosome.");
        verbose && exit(verbose);
        # Not needed anymore
        cells <- NULL;
        next;
      }
    }

    verbose && enter(verbose, "Identifying samples that have the baseline ploidy and those that have not");
    ploidies <- sapply(this, FUN=getPloidy, chromosome=chromosome,
                                              defaultValue=defaultPloidy);
    isBaseline <- (ploidies == ploidy);
    nB <- sum(isBaseline, na.rm=TRUE);
    # Number of samples with non-baseline ploidies.
    nM <- n - nB;
    verbose && printf(verbose, "Number of samples with ploidy %d: %d\n",
                                                              ploidy, nB);
    verbose && printf(verbose, "Number of other samples: %d\n", nM);

    # Assert that there are samples with the baseline ploidy
    if (nB == 0) {
      throw("Cannot estimate baseline signals. No samples with ploidy ", ploidy, " available.");
    }
    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Baseline samples
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Processing samples with baseline ploidy");

    verbose && enter(verbose, "Extracting subset of samples");
    # Baseline samples
    csB <- extract(this, which( isBaseline), onDuplicates="error");
    verbose && printf(verbose, "Baseline samples (with ploidy %d):\n", ploidy);
    verbose && print(verbose, csB);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating average");
    csBavg <- getAverageFile(csB, indices=cells, force=force, verbose=less(verbose));
    verbose && exit(verbose);

    verbose && enter(verbose, "Reading the average signals");
    muBs <- getData(csBavg, indices=cells, fields="intensities", verbose=less(verbose))$intensities;
    # Not needed anymore
    csBavg <- NULL;
    verbose && str(verbose, muBs);
#    verbose && cat(verbose, "Summary of log2(mu2s)");
#    verbose && print(verbose, summary(log2(muBs)));
    verbose && exit(verbose);

    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Other samples?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (nM > 0) {
      verbose && enter(verbose, "Processing samples with non-baseline ploidies");
      verbose && enter(verbose, "Extracting subset of samples");
      csM <- extract(this, which(!isBaseline), onDuplicates="error");
      verbose && cat(verbose, "All other samples:");
      verbose && print(verbose, csM);
      verbose && exit(verbose);

      verbose && enter(verbose, "Calculating average");
      csMavg <- getAverageFile(csM, indices=cells, force=force, verbose=less(verbose));
      verbose && exit(verbose);

      verbose && enter(verbose, "Reading the average signals");
      muMs <- getData(csMavg, indices=cells, fields="intensities", verbose=less(verbose))$intensities;
      # Not needed anymore
      csMavg <- NULL;
      verbose && str(verbose, muMs);
      verbose && exit(verbose);

      verbose && exit(verbose);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Estimating the shift
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Note, all of this is on the intensity and not the log scale.

      verbose && enter(verbose, "Estimating the baseline bias correction");

      # 1) Get the differences the two groups for each locus.
      cs <- (muBs / muMs);
      verbose && str(verbose, cs);
      verbose && cat(verbose, "Summary of log2(cs*)");
      verbose && print(verbose, summary(log2(cs)));

      # 2) Get the average difference across all loci.
      c <- median(cs, na.rm=TRUE);
      # Not needed anymore
      cs <- NULL;
      verbose && printf(verbose, "Bias correction: log2(c*)=%.3f\n", log2(c));
      verbose && exit(verbose);


      # 3) Weighted average of the two groups
      verbose && enter(verbose, "Estimating the weighted average of the two groups at each locus");

      # The estimate of the baseline according to the non-baseline samples
      muBs2 <- muMs * c;
      # Not needed anymore
      muMs <- NULL;
      verbose && cat(verbose, "Summary of log2(muBs*)");
      verbose && print(verbose, summary(log2(muBs2)));

      # The weights for the two groups
      wB <- nB/n;
      wM <- 1-wB;

      ds <- wB*muBs + wM*muBs2;
      # Not needed anymore
      muBs <- muBs2 <- NULL;
      verbose && exit(verbose);
    } else {
      ds <- muBs;
      verbose && cat(verbose, "All samples have baseline ploidy and are used to estimate the baseline signals.");
    }

    verbose && cat(verbose, "Summary of baseline signals log2(ds)");
    verbose && print(verbose, summary(log2(ds)));

    verbose && enter(verbose, "Storing baseline signals");
    ds <- cbind(intensities=ds, cell=cells);
    muBs <- updateDataFlat(csBaseline, data=ds, verbose=less(verbose));
    # Not needed anymore
    ds <- NULL;
    verbose && exit(verbose);

    # Mark cells as done
    todo <- setdiff(todo, cells);

    # Not needed anymore
    cells <- NULL;

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # for (chromosome ...)


  if (all) {
    verbose && enter(verbose, "Calculate the average signal for all cells not on a chromosome");
    cells <- todo;
    # Not needed anymore
    todo <- NULL;
    ncells <- length(cells);
    verbose && cat(verbose, "Number of remaining cells: ", length(cells));
    if (ncells > 0) {
      verbose && enter(verbose, "Checking for non-estimated cells");
      keep <- todo[cells];
      nkeep <- length(keep);
      cells <- cells[keep];
      # Not needed anymore
      keep <- NULL;

      verbose && printf(verbose, "Found %d (%.1f%%) non-estimated loci.\n",
                                                  nkeep, 100*nkeep/ncells);
      verbose && exit(verbose);

      if (nkeep > 0) {
        csRavg <- getAverageFile(this, indices=cells, force=force, verbose=less(verbose));

        verbose && enter(verbose, "Reading the average signals");
        ds <- getData(csRavg, indices=cells, fields="intensities", verbose=less(verbose))$intensities;
        # Not needed anymore
        csRavg <- NULL;
        verbose && str(verbose, ds);
        verbose && exit(verbose);

        verbose && enter(verbose, "Storing baseline signals");
        ds <- cbind(intensities=ds, cell=cells);
        muBs <- updateDataFlat(csBaseline, data=ds, verbose=less(verbose));
        # Not needed anymore
        ds <- NULL;
        verbose && exit(verbose);
      } # if (nkeep > 0)
    } # if (ncells > 0)
  } # if (all)

  verbose && exit(verbose);

  csBaseline;
}, protected=TRUE) # calculateBaseline()




############################################################################
# HISTORY:
# 2008-04-17
# o Now calling getUnitsOnChromosome() with argument 'chromosomes'.
# 2007-06-11
# o Removed never used 'muRs'.
# 2007-03-30
# o BUG FIX: Sometimes cell indices would become NAs.
# 2007-03-24
# o Now the average of non-genomic cells are also calculated if 'all=TRUE'.
#   Thus, if all samples have baseline ploidy, calculateBaseline() and
#   getAverage() should give the same results.
# 2007-03-22
# o TO DO: Estimate standard errors just like getAverage() does.
# o Unless 'force=TRUE', only cells for which the average has not already
#   been estimate are estimated.  This way calculateBaseline() is faster
#   if called multiple times.
# o Added getBaseline().
# o First working version of calculateBaseline(). Method now creates a CEL
#   files to store the estimates.
# 2007-03-16
# o Created.  See ploidy4.ps paper.
############################################################################
