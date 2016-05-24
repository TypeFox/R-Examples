setMethodS3("plot", "CopyNumberSegmentationModel", function(x, xlim=NULL, ..., pixelsPerMb=3, zooms=2^(0:6), pixelsPerTick=2.5, height=400, xmargin=c(50,50), imageFormat="current", skip=TRUE, path=NULL, callList=NULL, verbose=FALSE) {
  # To please R CMD check.
  this <- x;

  # If 'RColorBrewer' is missing, the below tryCatch() will plot the framwork
  # but not the data points and "nicely" catch the error and report it in
  # the verbose output.  This will cause PNGs with no data points.  Thus,
  # here we assert that RColorBrewer is available.
  # See thread '[aroma.affymetrix] No plot CNV analysis (Myriam)' on
  # Jun 22-July 1, 2009. /HB 2009-07-01
  pkg <- "RColorBrewer";
  require(pkg, character.only=TRUE) || throw("Package not loaded: ", pkg);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  drawXAxisRuler <- function(xrange=NULL, ticksBy=1) {
    xlim <- par("usr")[1:2];
    if (is.null(xrange))
      xrange <- xlim;

    for (kk in 1:3) {
      at <- seq(from=xrange[1], to=xrange[2], by=ticksBy*c(1,5,10)[kk]);
      keep <- (at >= xlim[1] & at <= xlim[2]);
      at <- at[keep];
      tcl <- c(0.2,0.4,0.6)[kk];
      lwd <- c(1,1,2)[kk];
      for (ss in c(1,3))
        axis(side=ss, at=at, tcl=tcl, lwd=lwd, labels=FALSE);
    }
    cxy <- par("cxy");
    text(x=at, y=par("usr")[3]-0.5*cxy[2], labels=at, srt=90,
                                       adj=1, cex=1, xpd=TRUE);
  } # drawXAxisRuler()



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pixelsPerMb':
  pixelsPerMb <- Arguments$getDouble(pixelsPerMb, range=c(0.001,9999));

  # Argument 'zooms':
  zooms <- Arguments$getIntegers(zooms, range=c(1,9999));
  zooms <- unique(zooms);

  # Argument 'pixelsPerTick':
  pixelsPerTick <- Arguments$getDouble(pixelsPerTick, range=c(1,256));

  # Argument 'height':
  height <- Arguments$getInteger(height, range=c(1,4096));

  # Argument 'callList':
  chipTypes <- getChipTypes(this);
  if (length(callList) > 0) {
    if (!is.list(callList))
      callList <- list(callList);

    if (length(callList) != length(chipTypes)) {
      throw("Number of elements in argument 'callList' does not match the number of chip types: ", length(callList), " != ", length(chipTypes));
    }

    if (is.null(names(callList)))
      names(callList) <- chipTypes;

    for (chipType in chipTypes) {
      callSet <- callList[[chipType]];
      if (!is.null(callSet)) {
        callSet <- Arguments$getInstanceOf(callSet, "GenotypeCallSet");

        if (getChipType(callSet) != chipType) {
          throw("Argument 'callList' contains a GenotypeCallSet for a different chip type than the corresponding copy-number set: ", getChipType(callSet), " != ", chipType);
        }
      }
    }
  }

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
  # Output path
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The report path
  if (is.null(path)) {
    path <- getReportPath(this);
  }
  path <- Arguments$getWritablePath(path);

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
    pngDev <- findPngDevice(transparent=FALSE);
    plotDev <- pngDev;
    if (identical(pngDev, png2))
      resScale <- 2;
  }


  # Get chip type (used to annotate the plot)
  chipType <- getChipType(this);

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

      # Figure out what chromosome is fitted
#      chromosome <- unique(fit$profileValues$Chromosome); # Should only be one!
##      chromosome[chromosome == "23"] <- "X";  # TODO

      # Infer the length (in bases) of the chromosome
##      chromosomeIdx <- match(chromosome, getChromosomes(this));
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

        # Calculate MBs per ticks
        ticksBy <- 10^ceiling(log10(pixelsPerTick / (zoom * pixelsPerMb)));

        # Calculate width in pixels from MBs
        width <- round(zoom * widthMb * pixelsPerMb + sum(xmargin));

        # Plot to PNG file
        verbose && printf(verbose, "Pathname: %s\n", pathname);
        verbose && printf(verbose, "Dimensions: %dx%d\n", width, height);
        verbose && printf(verbose, "Ticks by: %f\n", ticksBy);

        if (!is.null(plotDev))
          plotDev(pathname, width=width, height=height);
        tryCatch({
          verbose && enter(verbose, "Plotting graph");
          opar <- par(xaxs="r");
          suppressWarnings({
            # Create empty plot
            verbose && enter(verbose, "Creating empty plot");
            newPlot(this, xlim=xlim, ylim=ylim, flavor="ce", unit=unit, ...);
            verbose && exit(verbose);

#            plotProfile2(fit, unit=unit, ylim=c(-1,+1)*3, xmargin=xmargin, resScale=resScale, flavor="ce", Bkp=FALSE, Smoothing=NULL, cnLevels=NULL, plotband=FALSE, ...);

            verbose && enter(verbose, "Adding axes and rulers");
            # Add ruler
            drawXAxisRuler(xrange=c(0,nbrOfBases)/10^unit, ticksBy=ticksBy);
            verbose && exit(verbose);

            # Add cytoband to graph (optional; class specific)
            if (!identical(this$.plotCytoband, FALSE)) {
              verbose && enter(verbose, "Adding cytoband");
              drawCytoband(this, chromosome=chromosome, unit=unit);
              verbose && exit(verbose);
            }

            # Add CN=1,2,3 lines to graph
            verbose && enter(verbose, "Adding CN grid lines");
            cnLevels <- c(1/2,1,3/2);
            for (level in cnLevels) {
              abline(h=log2(level), col="blue", lty=2);
            }
            verbose && exit(verbose);

            # Extract raw CN estimates
            rawCns <- extractRawCopyNumbers(fit);
            verbose && print(verbose, rawCns, level=-50);

            # Add number-of-loci annotation to graph
            n <- nbrOfLoci(rawCns, na.rm=TRUE);
            stext(text=sprintf("n=%d", n), side=4, pos=0, line=0, cex=0.8);

            # Plot raw CNs data points (and highlight outliers)
            # Generic function
            verbose && enter(verbose, "Adding data points");
            pointsRawCNs(fit, unit=unit, ...);
#            points(rawCns, xScale=1/10^unit, pch=20, col="black");
            verbose && exit(verbose);

            # Draw CNRs
            verbose && enter(verbose, "Adding segmentation results");
            cnRegions <- extractCopyNumberRegions(fit);
            verbose && print(verbose, cnRegions, level=-50);
            drawLevels(cnRegions, lwd=4, col="black", xScale=1/10^unit);
            verbose && exit(verbose);

            # Model-specific annotations (optional; class specific)
            verbose && enter(verbose, "Adding additional annotations");
            drawExtraAnnotations(fit);
            verbose && exit(verbose);

            # Add genotype call tracks, if available (optional; class specific)
            onFitAddGenotypeCalls(fit, callList=callList, arrayName=arrayName, resScale=10^unit*resScale, ...);
          });

          # Add chip-type annotation
          stext(chipType, side=4, pos=1, line=0, cex=0.8);

          verbose && exit(verbose);
        }, error = function(ex) {
          print(ex);
        }, finally = {
          par(opar);
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
}) # plot()




##############################################################################
# HISTORY:
# 2013-07-20
# o CLEANUP: Replaces an x11() with a dev.new().
# 2010-12-02
# o BUG FIX: plot() for CopyNumberSegmentationModel would throw exception
#   "Cannot infer number of bases in chromosome. No such chromosome: 25"
#   for chromosome 25.
# o CLEANUP: Dropped any usage getChromosomeLength().
# 2010-11-23
# o Added more verbose output to plot() of CopyNumberSegmentationModel.
# 2009-07-01
# o ROBUSTNESS/BUG FIX: Now plot() of CopyNumberSegmentationModel asserts
#   that the RColorBrewer package is avaiable at the very beginning.  This
#   will avoid generating image files where the data points are missing.
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
