###########################################################################/**
# @RdocClass GladModel
#
# @title "The GladModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Gain and Loss Analysis of DNA regions
#  (GLAD) model [1].
#  This class can model chip-effect estimates obtained from multiple
#  chip types, and not all samples have to be available on all chip types.
# }
#
# @synopsis
#
# \arguments{
#   \item{cesTuple}{A @see "CopyNumberDataSetTuple".}
#   \item{...}{Arguments passed to the constructor of
#              @see "CopyNumberSegmentationModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   Data from multiple chip types are combined "as is".  This is based
#   on the assumption that the relative chip effect estimates are
#   non-biased (or at the equally biased across chip types).
#   Note that in GLAD there is no way to down weight certain data points,
#   which is why we can control for differences in variance across
#   chip types.
# }
#
# \section{Benchmarking}{
#   In high-density copy numbers analysis, the most time consuming step
#   is fitting the GLAD model.  The complexity of the model grows
#   more than linearly (squared? exponentially?) with the number of data
#   points in the chromosome and sample being fitted.  This is why it
#   take much more than twice the time to fit two chip types together
#   than separately.
# }
#
# @author
#
# \seealso{
#  @see "CopyNumberSegmentationModel".
# }
#
# \references{
#  [1] Hupe P et al. \emph{Analysis of array CGH data: from signal ratio to
#      gain and loss of DNA regions}. Bioinformatics, 2004, 20, 3413-3422.\cr
# }
#*/###########################################################################
setConstructorS3("GladModel", function(cesTuple=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(cesTuple)) {
    requireWithMemory("GLAD") || throw("Package not loaded: GLAD");
  }

  extend(CopyNumberSegmentationModel(cesTuple=cesTuple, ...), "GladModel")
})

setMethodS3("getFitFunction", "GladModel", function(this, ...) {
  segmentByGLAD;
}, protected=TRUE)


setMethodS3("writeRegions", "GladModel", function(this, arrays=NULL, format=c("xls", "wig"), digits=3, ..., oneFile=TRUE, skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  arrays <- indexOf(this, arrays);

  # Argument 'format':
  format <- match.arg(format);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Setup
  fullname <- getFullName(this);
  arrayNames <- getNames(this);

  path <- getPath(this);
  path <- Arguments$getWritablePath(path);

  if (oneFile) {
    filename <- sprintf("%s,regions.%s", fullname, format);
    pathname <- filePath(path, filename);
    pathname <- Arguments$getWritablePathname(pathname);
    if (!skip && isFile(pathname)) {
      file.remove(pathname);
    }
  }

  res <- list();
  for (aa in seq_along(arrays)) {
    array <- arrays[aa];
    name <- arrayNames[array];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                               aa, name, length(arrays)));
    df <- getRegions(this, arrays=array, ..., verbose=less(verbose))[[1]];
    names(df) <- gsub("(Smoothing|mean)", "log2", names(df))
    verbose && str(verbose, df)

    if (nrow(df) > 0) {
      if (identical(format, "xls")) {
        # Append column with sample names
        df <- cbind(sample=name, df);
      } else if (identical(format, "wig")) {
        # Write a four column WIG/BED table
        df <- df[,c("Chromosome", "start", "stop", "log2")];

        # In the UCSC Genome Browser, the maximum length of one element
        # is 10,000,000 bases.  Chop up long regions in shorter contigs.
        verbose && enter(verbose, sprintf("Chopping up too long segment"));
        MAX.LENGTH = 10e6-1;
        start <- df[,"start"];
        stop <- df[,"stop"];
        len <- stop-start;
        tooLong <- which(len > MAX.LENGTH);
        if (length(tooLong) > 0) {
          dfXtra <- NULL;
          for (rr in tooLong) {
            x0 <- start[rr];
            while (x0 < stop[rr]) {
              x1 <- min(x0 + MAX.LENGTH, stop[rr]);
              df1 <- df[rr,];
              df1[,"start"] <- x0;
              df1[,"stop"] <- x1;
              dfXtra <- rbind(dfXtra, df1);
              x0 <- x1+1;
            }
          }
          df <- df[-tooLong,];
          df <- rbind(df, dfXtra);
          # Not needed anymore
          dfXtra <- NULL;
          row.names(df) <- seq_len(nrow(df));
        }
        verbose && exit(verbose);
        # Make sure the items are ordered correctly
        chrIdx <- as.integer(df[,"Chromosome"]);
        o <- order(chrIdx, df[,"start"]);
        df <- df[o,];

        # All chromosomes should have prefix 'chr'.
        chrIdx <- as.integer(df[,"Chromosome"]);
        ## df[chrIdx == 23,"Chromosome"] <- "X"; ## REMOVED 2007-03-15
        df[,"Chromosome"] <- paste("chr", df[,"Chromosome"], sep="");
      }

      # Apply digits
      for (cc in seq_len(ncol(df))) {
        value <- df[,cc];
        if (is.double(value)) {
          df[,cc] <- round(value, digits=digits);
        }
      }
    } # if (nrow(df) > 0)

    if (!oneFile) {
      savename <- name;
      filename <- sprintf("%s,regions.%s", savename, format);
      pathname <- filePath(path, filename);
      if (!oneFile && !skip && isFile(pathname))
        file.remove(pathname);
    }

    # Writing to file
    verbose && cat(verbose, "Pathname: ", pathname);
    if (identical(format, "xls")) {
      col.names <- (array == arrays[1]);
      suppressWarnings({
        write.table(df, file=pathname, sep="\t", col.names=col.names, row.names=FALSE, quote=FALSE, append=oneFile);
      })
    } else if (identical(format, "wig")) {
      # Write track control
      trackAttr <- c(type="wiggle_0");
      trackAttr <- c(trackAttr, name=sprintf("\"%s\"", name));
      trackAttr <- c(trackAttr, group="\"GLAD regions\"");
      trackAttr <- c(trackAttr, priority=array);
      trackAttr <- c(trackAttr, graphType="bar");
      trackAttr <- c(trackAttr, visibility="full");
      trackAttr <- c(trackAttr, maxHeightPixels="128:96:64");
      trackAttr <- c(trackAttr, yLineOnOff="on");
# HARD WIRED FOR NOW.  TO DO /hb 2006-11-27
col <- c("117,112,179", "231,41,138");
ylim <- c(-1,1);
      if (!is.null(col)) {
        trackAttr <- c(trackAttr, color=col[1], altColor=col[2]);
      }
      if (!is.null(ylim)) {
        trackAttr <- c(trackAttr, autoScale="off",
              viewLimits=sprintf("%.2f:%.2f ", ylim[1], ylim[2]));
      }
      trackAttr <- paste(names(trackAttr), trackAttr, sep="=");
      trackAttr <- paste(trackAttr, collapse=" ");
      trackAttr <- paste("track ", trackAttr, "\n", sep="");
      verbose && str(verbose, trackAttr);
      cat(file=pathname, trackAttr, append=oneFile);

      # Write data
      verbose && str(verbose, df);
      write.table(df, file=pathname, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE);
    }
    verbose && exit(verbose);
    res[[array]] <- df;
  }

  invisible(pathname);
})


##############################################################################
# HISTORY:
# 2012-10-21
# o ROBUSTNESS: Now using Arguments$getWritablePath() instead of mkdirs()
#   whereever applicable.  Soon, when Arguments$getWritablePath() will
#   assert that (i) the requested path exists/or created, and (ii) allowing
#   for small delay between creating the path and testing for it's existance
#   on slow file systems.
# 2009-05-16
# o Added getFitFunction().  Removed fitOne().
# 2007-09-04
# o Move plot() of GladModel to its own file.
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
