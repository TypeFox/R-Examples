setMethodS3("onFitAddGenotypeCalls", "default", function(fit, ...) {
}, protected=FALSE)

setMethodS3("onFitAddGenotypeCalls", "GladModel", function(gladFit, callList, arrayName, resScale=1e6, ylim=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  gladFit <- Arguments$getInstanceOf(gladFit, "profileCGH");

  # Nothing to do?
  if (is.null(callList)) {
    return();
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  pv <- gladFit$profileValues;

  # Identify chip types included
  chipTypes <- sort(unique(pv$chipType));
  if (length(chipTypes) == 0) {
    warning("Could not add genotype call tracks.  No chip type information available.");
    return();
  }

  callCols <- c("-"="lightgray", AA="red", AB="blue", BB="red", NC="orange");

  # Extract the chromosome from the GLAD fit object
  pv <- gladFit$profileValues;

  verbose && enter(verbose, "Adding genotype calls");

  # Figure out where to put the genotype track?
  if (is.null(ylim)) {
    ylim <- par("usr")[3:4];
    ylim <- ylim + c(+1,-1)*0.04*diff(ylim);
    ylim <- ylim + c(+1,-1)*0.04*diff(ylim);
  }

  # Get (chip type, unit) information
  chipType <- pv$chipType;
  unit <- pv$unit;

  # Add call tracks for each chip type available
  for (chipType in chipTypes) {
    # Got genotype calls for this chip type?
    callSet <- callList[[chipType]];
    if (is.null(callSet))
      next;

    # Got chip-effect estimates for this chip type?
    idxs <- which(pv$chipType == chipType);
    if (length(idxs) == 0)
      next;

    # Got genotype calls for this array?
    if (!arrayName %in% getNames(callSet))
      next;

    # Get subset of genotype calls for this array & chromosome.
    units <- pv$unit[idxs];
    call <- callSet[units, arrayName];
    call <- as.character(call);

    # Get the positions of these calls
    x <- pv$PosBase[idxs];
    x <- x/resScale;

    # Plot the homozygote/heterozygote genotype tracks
    y <- rep(ylim[1], length(callCols));
    names(y) <- names(callCols);
    y["AB"] <- y["AB"] + 0.02*diff(ylim);
    y <- y[call];
    points(x,y, pch=".", cex=2, col=callCols[call]);

    # Not needed anymore
    idxs <- call <- callSet <- units <- x <- y <- NULL;
  } # for (chipType ...)

  verbose && exit(verbose);
}, protected=TRUE) # onFitAddGenotypeCalls()




##############################################################################
# HISTORY:
# 2012-11-29
# o Added missing 'verbose' argument to onFitAddGenotypeCalls() for GladModel.
# 2007-09-04
# o Made plot() of GladModel more generic.  I'm trying to move this method
#   up to CopyNumberSegmentationModel, which means several of the plot
#   functions implemented in GLAD has to be generalized and rewritten.
# 2007-08-20
# o Initial tests show that the updated GladModel gives identical results.
# o Now GladModel inherits from CopyNumberSegmentationModel.
# 2007-05-10
# o BUG FIX: getRegions() and getLog2Ratios() would give an error if a subset
#   of the chromosomes where queried.
# o Added more verbose output to getRegions().
# 2007-04-12
# o Now plot() of the GladModel writes the chip type annotation in black and
#   not in gray as before.
# 2007-03-24
# o BUG FIX: getPath() created the root path before trying to expand
#   Windows shortcuts.
# 2007-03-19
# o Now asterisk tags are handles dynamically, and not by the constructor.
# o Added getAsteriskTags().
# o Now the constructor expects a ChipEffectSetTuple instead of a list of
#   ChipEffectSet:s.
# o Updated code to internally make use of the ChipEffectSetTuple class.
# 2007-03-15
# o Updated the GladModel to only work with chromosome indices (integers).
# o Now the GladModel infers the set of possible chromosomes from the
#   GenomeInformation file.
# 2007-03-12
# o BUG FIX: getFullNames() of GladModel would give 'Error in getName(ceList[[
#   1]]) : no applicable method for "getName"' if there was not hybridization
#   for the *first* chip type in a set of multiple chip types.
# o BUG FIX: fit() of GladModel would give 'Error in FUN(X[[2]], ...) : no
#   applicable method for "getTags"' if there were not data files for all
#   chip types.
# 2007-02-20
# o Added getFullNames(), which for each tuple (across chip types) returns the
#   sample name of the tuple, together with all *common* tags across all
#   chip types.  Tags existing in only some of the chip types are ignored.
# 2007-02-16
# o Now the default version of the human genome is 'hg17' and not 'hg18'.
#   The reason for this is that the dChip annotation files are 'hg17'. We
#   still have to figure out a way to do version control for this in the
#   package.  Maybe it won't be a problem as soon as we start using the
#   annotation packages of Bioconductor.  On the to do list...
# o Added arguments 'organism' and 'db' to getRegions().
# 2007-02-15
# o Now getChipTypes() sorts the chip types in lexicographic order before
#   merging.  This guarantees the same result regardsless of order of the
#   input list.
# o Added getReportPath().
# o Path is now back to <rootPath>/<data set>,<tags>/<chipType>/.
# o Reports are written to reports/<data set>/<tags>/<chipType>/glad/.
# 2007-02-06
# o Updated the path to <rootPath>/<dataSetName>/<tags>/<chipType>/<set>/.
# 2007-01-25
# o Added so that plot() generates fixed sized horizontal margins (50px),
#   so that we can infer the relative genomic location from the horisontal
#   pixel location.
# 2007-01-21
# o Added a better error message when plot() fails to locate hgChromosomes.txt.
# 2007-01-17
# o Now argument 'arrays' can be either a vector of indices or array names,
#   or NULL.
# o Added indexOfArrays().
# 2007-01-16
# o Now NULL values for arguments 'arrays' and 'chromosomes' of fit() defaults
#   to all arrays and all chromosomes, respectively.
# o BUG FIX: writeRegions() would give an error if no regions was found.
# 2007-01-15
# o Now fit(..., force=TRUE) also calls getReferenceFiles(..., force=force).
# o Added some more Rdoc comments.
# 2007-01-10
# o Now plot() of GladModel is search for 'hgChromosomes.txt' in both
#   annotations/ and the package installation directory.
# 2007-01-07
# o Renamed MultiGladModel to GladModel fully replacing the older class.
# 2006-12-20
# o Now the class accepts any ChipEffectSet, not only CnChipEffectSet objects.
#   CnChipEffectSet objects are still validated specially, if used.
# 2006-12-17
# o BUG FIX: The new fitDone() in plot() choked on chr 23 (should be 'X').
# 2006-12-15
# o This class should be considered temporary, because we might design a
#   ChipEffectSet class that can contain multiple chip types, but treated as
#   if it contained one chip type, so it can be passed to the current
#   GladModel class.  However, such a class design will require multiple
#   inheritance etc, which will take time to develope.
# o Created from GladModel.R with history as below:
# 2006-11-29
# o Added chip type annotation to plot() and option to plot to screen.
# 2006-11-27
# o Added argument 'flat' to getRegions().
# 2006-11-24
# o Now the fit() function of GladModel stores the profileCGH object as a
#   binary XDR file in the default path, see getPath().
# 2006-11-23
# o Added writeWig().
# 2006-11-22
# o Added writeRegions().
# o Added fit(), plot(), and getRegions().
# o Re-created from the CnAnalyzer class from 2006-10-31.
##############################################################################
