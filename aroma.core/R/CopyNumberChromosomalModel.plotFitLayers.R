setMethodS3("plotFitLayers", "CopyNumberChromosomalModel", function(this, FUN, path, xlim=NULL, ..., pixelsPerMb=3, zooms=2^(0:6), height=400, xmargin=c(50,50), imageFormat="current", transparent=FALSE, skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'FUN':
  if (!is.function(FUN)) {
    throw("Arguments 'FUN' is not a function: ", class(FUN)[1]);
  }

  # Argument 'pixelsPerMb':
  pixelsPerMb <- Arguments$getDouble(pixelsPerMb, range=c(0.001,9999));

  # Argument 'zooms':
  zooms <- Arguments$getIntegers(zooms, range=c(1,9999));
  zooms <- unique(zooms);

  # Argument 'height':
  height <- Arguments$getInteger(height, range=c(1,4096));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Get genome annotation data (chromosome lengths etc)
  genome <- getGenomeData(this);

  # In units of 10^unit bases (default is Mb)
  unit <- 6;

  # Default 'ylim'
  ylim <- c(-1,+1)*3;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the PNG device
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(imageFormat)) {
    imageFormat <- "current";
  }

  resScale <- 1;
  if (identical(imageFormat, "current")) {
    plotDev <- NULL;
    zooms <- zooms[1];
  } else if (identical(imageFormat, "screen")) {
    screenDev <- function(pathname, width, height, ..., xpinch=50, ypinch=xpinch) {
      # Dimensions are in pixels. Rescale to inches
      width <- width/xpinch;
      height <- height/ypinch;
      dev.new(width=width, height=height, ...);
    }

    # When plotting to the screen, use only the first zoom
    zooms <- zooms[1];
    plotDev <- screenDev;
  } else if (identical(imageFormat, "png")) {
    pngDev <- findPngDevice(transparent=TRUE);
    plotDev <- pngDev;
    if (identical(pngDev, png2))
      resScale <- 2;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Define the plot function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hookName <- "onFit.CopyNumberSegmentationModel";
  on.exit({
    setHook(hookName, NULL, action="replace");
  })

  setHook(hookName, function(fit, chromosome, fullname) {
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    tryCatch({
      # Extract the array name from the full name
      arrayFullName <- gsub("^(.*),chr[0-9][0-9].*$", "\\1", fullname);
      arrayName <- gsub("^([^,]*).*$", "\\1", arrayFullName);

      # Infer the length (in bases) of the chromosome
      nbrOfBases <- genome$nbrOfBases[chromosome];
      widthMb <- nbrOfBases / 10^unit;

      # Argument 'xlim' missing?
      if (is.null(xlim)) {
        xlim <- c(0, widthMb);
      }

      verbose && enter(verbose, sprintf("Plotting %s for chromosome %02d [%.2fMB]", arrayName, chromosome, widthMb));

      for (zz in seq_along(zooms)) {
        zoom <- zooms[zz];

        # Create the pathname to the file
        imgName <- sprintf("%s,chr%02d,x%04d.%s",
                          arrayFullName, chromosome, zoom, imageFormat);
        pathname <- filePath(path, imgName);

        # pngDev() (that is bitmap()) does not accept spaces in pathnames
        pathname <- gsub(" ", "_", pathname);
        if (!imageFormat %in% c("screen", "current")) {
          if (skip && isFile(pathname)) {
            next;
          }
        }

        # Calculate width in pixels from MBs
        width <- round(zoom * widthMb * pixelsPerMb + sum(xmargin));

        # Plot to PNG file
        verbose && printf(verbose, "Pathname: %s\n", pathname);
        verbose && printf(verbose, "Dimensions: %dx%d\n", width, height);

        args <- list(cns=this, fit=fit, chromosome=chromosome, xlim=xlim, ylim=ylim, unit=unit, width=width, height=height, zoom=zoom, pixelsPerMb=pixelsPerMb, nbrOfBases=nbrOfBases, ..., verbose=less(verbose,1));

        if (!is.null(plotDev))
          plotDev(pathname, width=width, height=height);

        if (transparent) {
          par(bg=NA, xaxs="r");
        } else {
          par(xaxs="r");
        }

        tryCatch({
          do.call(FUN, args=args);
        }, error = function(ex) {
          print(ex);
        }, finally = {
          if (!imageFormat %in% c("screen", "current"))
            dev.off();
        });
      } # for (zz in ...)
    }, error = function(ex) {
      cat("ERROR caught in ", hookName, "():\n", sep="");
      print(ex);
    }, finally = {
      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose);
    }) # tryCatch()
  }, action="replace")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Start fitting and plotting
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit(this, ..., .retResults=FALSE, verbose=verbose);

  invisible();
}, protected=TRUE) # plotFitLayers()



##############################################################################
# HISTORY:
# 2013-07-20
# o CLEANUP: Replaces an x11() with a dev.new().
# 2010-12-02
# o CLEANUP: Dropped any usage getChromosomeLength().
# 2007-10-09
# o Added plotCytobandLayers().
# 2007-09-15
# o Now the cytoband is only drawn for some genomes, which currently is
#   hardwired to the "Human" genome.
# 2007-09-04
# o Finally, now plot() works pretty much the same for GladModel as for
#   the new CbsModel.
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
