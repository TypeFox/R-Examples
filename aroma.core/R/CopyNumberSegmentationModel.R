##########################################################################/**
# @RdocClass CopyNumberSegmentationModel
#
# @title "The CopyNumberSegmentationModel class"
#
# \description{
#  @classhierarchy
#
#  This \emph{abstract} class represents a copy-number segmentation model.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to constructor
#      @see "CopyNumberChromosomalModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("CopyNumberSegmentationModel", function(...) {
  extend(CopyNumberChromosomalModel(...), "CopyNumberSegmentationModel");
})


setMethodS3("getAsteriskTags", "CopyNumberSegmentationModel", function(this, collapse=NULL, ..., tag=NULL) {
  if (is.null(tag)) {
    # Infer 'GLAD' from GladModel, 'CBS' from CbsModel,
    # 'HAARSEG' from HaarSegModel, and so on.
    tag <- class(this)[1];
    tag <- gsub("Model$", "", tag);
    tag <- toupper(tag);
  }

  tags <- tag;

  # Add class-specific tags
  if (isPaired(this))
    tags <- c(tags, "paired");

  # Collapse?
  tags <- Arguments$getTags(tags, collapse=collapse);

  tags;
}, protected=TRUE)



setMethodS3("getTags", "CopyNumberSegmentationModel", function(this, collapse=NULL, ...) {
  tags <- getTags(getSetTuple(this), collapse=collapse, ...);

  # Add model tags
  tags <- c(tags, this$.tags);

  # In case this$.tags is not already split
  tags <- strsplit(tags, split=",", fixed=TRUE);
  tags <- unlist(tags);

  # Update default tags
  asteriskTags <- paste(getAsteriskTags(this)[-1], collapse=",");
  if (length(asteriskTags) == 0)
    asteriskTags <- "";
  tags[tags == "*"] <- asteriskTags;

  tags <- Arguments$getTags(tags, collapse=NULL);

  # Get unique tags
  tags <- locallyUnique(tags);

  # Collapsed or split?
  tags <- Arguments$getTags(tags, collapse=collapse);

  tags;
})



setMethodS3("getFitFunction", "CopyNumberSegmentationModel", abstract=TRUE, protected=TRUE);



###########################################################################/**
# @RdocMethod fit
#
# @title "Fits the model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{A @vector of array indices specifying which arrays to
#    be considered.  If @NULL, all are processed.}
#   \item{chromosome}{A @vector of chromosomes indices specifying which
#     chromosomes to be considered.  If @NULL, all are processed.}
#   \item{force}{If @FALSE, the model will not be fitted again if it was
#     already fitted.}
#   \item{...}{Additional arguments passed to the segmentation method for
#     the @see "aroma.core::RawGenomicSignals".}
#   \item{.retResults}{If @TRUE, the segmentation fit structures are
#     returned for each fitted array and chromosome.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a named @list of named @lists.
# }
#
# \section{Additional arguments to the internal fit function}{
#   Arguments in \code{...} are passed down to the internal fit function,
#   which means that it is possible to fine tune even further.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fit", "CopyNumberSegmentationModel", function(this, arrays=NULL, chromosomes=getChromosomes(this), maxNAFraction=getMaxNAFraction(this), force=FALSE, ..., .retResults=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getTagsFromList <- function(files, ...) {
    keep <- (!sapply(files, FUN=is.null));
    files <- files[keep];

    # Get tags *common* across chip types
    tags <- lapply(files, FUN=function(file) {
      if (is.character(file)) return(file);
      getTags(file);
    });
    tags <- getCommonListElements(tags);
    tags <- unlist(tags, use.names=FALSE);
    # BEGIN: AFFX
    tags <- setdiff(tags, "chipEffects");
    # END: AFFX
    tags <- locallyUnique(tags);
    tags;
  } # getTagsFromList()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  if (identical(arrays, "fitted")) {
  } else {
    arrays <- indexOf(this, arrays);
  }

  allChromosomes <- getChromosomes(this);

  # Argument 'chromosomes':
  if (identical(chromosomes, "fitted")) {
  } else if (is.null(chromosomes)) {
    chromosomes <- getChromosomes(this);
  } else if (is.numeric(chromosomes)) {
    chromosomes <- Arguments$getChromosomes(chromosomes,
                                                range=range(allChromosomes));
##    chromosomes <- as.character(chromosomes);  ## TODO
##    chromosomes[chromosomes == "23"] <- "X";   ## TODO
    chromosomes <- intersect(chromosomes, allChromosomes);
  } else if (is.character(chromosomes)) {
    chromosomes <- Arguments$getChromosomes(chromosomes,
                                                range=range(allChromosomes));
##    chromosomes[chromosomes == "23"] <- "X";   ## TODO
    chromosomes <- intersect(chromosomes, getChromosomes(this));
  }

  # Argument 'maxNAFraction':
  if (!missing(maxNAFraction)) {
    msg <- sprintf("Argument 'maxNAFraction' to fit() of CopyNumberSegmentationModel is deprecated. Instead, specify when setting up the %s object.", class(this)[1]);
    warning(msg);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- getPath(this);
  path <- Arguments$getWritablePath(path);

  fitFcn <- getFitFunction(this, verbose=less(verbose, 50));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving tuple of reference data files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Note: Do *not* pass done 'force' to getReferenceFiles(), because then
  # it will calculate the average regardless of reference. /HB 2007-03-24
  refTuple <- getReferenceSetTuple(this, verbose=verbose);
  verbose && cat(verbose, "Using reference tuple:");
  verbose && print(verbose, refTuple);


  ## Call onFit() hook functions later?
  hookName <- "onFit.CopyNumberSegmentationModel"
  hasHooks <- (length(getHook(hookName)) > 0)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Array by array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- listenv();
  arrayNames <- getNames(this)[arrays];
  nbrOfArrays <- length(arrayNames);
  for (aa in seq_len(nbrOfArrays)) {
    array <- arrays[aa];
    arrayName <- arrayNames[aa];

    files <- getDataFileMatrix(this, array=array, verbose=less(verbose,5));

    # Extract the test and reference files
    ceList <- files[,"test", drop=FALSE];
    rfList <- files[,"reference", drop=FALSE];

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get tags for test sample
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get copy-number signal tags *common* across chip types
    ceTags <- getTagsFromList(ceList);
    verbose && cat(verbose, "Genomic-signal tags: ", paste(ceTags, collapse=","));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get tags for reference sample
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get copy-number signal tags *common* across chip types
    rfTags <- getTagsFromList(rfList);

    # HB 2007-02-19 To fix: Should the name and the tags of average files
    # be replaced?!? We do not get the same names if we specify the average
    # files explicitly or indirectly.

    # Add combined reference name
    names <- sapply(rfList, FUN=function(file) {
      if (is.character(file)) return(file);
      getName(file);
    });
    names <- mergeByCommonTails(names,"+");
    rfTags <- c(names, rfTags);
    rfTags <- unique(rfTags);
    if (length(rfTags) != 1L) {
      rfTags <- getChecksum(rfTags);
    }
    verbose && cat(verbose, "Reference tags: ", paste(rfTags, collapse=","));



    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Chromosome by chromosome
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    res[[arrayName]] %<=% {
      resArray <- list()
      for (chr in chromosomes) {
        verbose && enter(verbose,
                            sprintf("Array #%d ('%s') of %d on chromosome %s",
                                             aa, arrayName, nbrOfArrays, chr));

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Get pathname
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Add tags chrNN,<reference tags>
        tags <- c(ceTags, sprintf("chr%02d", chr));
        tags <- c(tags, rfTags);
        fullname <- paste(c(arrayName, tags), collapse=",");
        filename <- sprintf("%s.xdr", fullname);
        pathname <- filePath(path, filename);
        verbose && cat(verbose, "Pathname: ", pathname)

        # Already done?
        if (!force && isFile(pathname)) {
          verbose && cat(verbose, "Already done. Skipping.")
          fit <- loadObject(pathname)
        } else {
          # Time the fitting.
          startTime <- processTime();

          timers <- list(total=0, read=0, fit=0, write=0, gc=0);

          tTotal <- processTime();

          # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          # Get (x, M, stddev, chipType, unit) data from all chip types
          # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          tRead <- processTime();
          cn <- extractRawCopyNumbers(this, array=array, chromosome=chr, ...);
          timers$read <- timers$read + (processTime() - tRead);
          verbose && print(verbose, cn);


          # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          # Fit segmentation model
          # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          verbose && enter(verbose, "Calling model fit function");

          optArgs <- getOptionalArguments(this);
          verbose && cat(verbose, "Optional arguments (may be ignored/may give an error/warning):");
          verbose && str(verbose, optArgs);
          userArgs <- list(...);
          excl <- which(names(userArgs) == "maxNAFraction");
          if (length(excl) > 0L) userArgs <- userArgs[-excl];
          if (length(userArgs) > 0L) {
            verbose && cat(verbose, "User arguments (may be ignored/may give an error/warning):");
            verbose && str(verbose, userArgs);
          }
          args <- list(cn);
          args <- c(args, optArgs, userArgs);
          verbose && cat(verbose, "All arguments:");
          verbose && str(verbose, args);
          args <- c(args, list(...), list(verbose=less(verbose, 1)));
          tFit <- processTime();
          fit <- do.call(fitFcn, args);
          verbose && str(verbose, fit);
          timers$fit <- timers$fit + (processTime() - tFit);
          # Not needed anymore
          cn <- NULL;

          verbose && cat(verbose, "Class of fitted object: ", class(fit)[1]);
          verbose && printf(verbose, "Time to fit segmentation model: %.2fmin\n", timers$fit[3]/60);

          verbose && exit(verbose);


          verbose && enter(verbose, "Validate that it can be coerced");
          rawCns <- extractRawCopyNumbers(fit);
          nbrOfLoci <- nbrOfLoci(rawCns);
          verbose && print(verbose, rawCns);
          sigmaM <- estimateStandardDeviation(rawCns);
          verbose && printf(verbose, "Robust first-order standard deviation estimate: %g\n", sigmaM);
          cnRegions <- extractCopyNumberRegions(fit);
          verbose && print(verbose, cnRegions);
          # Not needed anymore
          rawCns <- cnRegions <- NULL;
          verbose && exit(verbose);


          # Garbage collection
          tGc <- processTime();
          gc <- gc();
          timers$gc <- timers$gc + (processTime() - tGc);
          verbose && print(verbose, gc);

          verbose && enter(verbose, "Saving to file");
          verbose && cat(verbose, "Pathname: ", pathname);
          tWrite <- processTime();
          saveObject(fit, file=pathname);
          timers$write <- timers$write + (processTime() - tWrite);
          verbose && exit(verbose);

          timers$total <- timers$total + (processTime() - tTotal);

          # Report time profiling
          totalTime <- processTime() - startTime;

          if (verbose) {
            t <- totalTime[3];
            printf(verbose, "Total time for chromosome %d: %.2fs == %.2fmin\n", chr, t, t/60);
            t <- totalTime[3]/nbrOfLoci;
            printf(verbose, "Total time per 1000 locus (with %d loci): %.2fs\n", nbrOfLoci, 1000*t);
            # Get distribution of what is spend where
            t <- lapply(timers, FUN=function(timer) unname(timer[3]));
            t <- unlist(t);
            t <- 100 * t / t["total"];
            printf(verbose, "Fraction of time spent on different tasks: Fitting: %.1f%%, Reading: %.1f%%, Writing: %.1f%%, Explicit garbage collection: %.1f%%\n", t["fit"], t["read"], t["write"], t["gc"]);
          }
        }

        ## Record segmentation results?
        if (.retResults) resArray[[chr]] <- fit

        ## Call onFit() hooks?
        if (hasHooks) {
          verbose && enter(verbose, sprintf("Calling %s() hooks", hookName))
          callHooks("onFit.CopyNumberSegmentationModel", fit=fit, chromosome=chr, fullname=fullname)
          verbose && exit(verbose)
        }

        # Not needed anymore
        fit <- NULL;

        verbose && exit(verbose);
      } # for (chr in ...)

      resArray
    } ## %<=%
  } # for (aa in ...)

  ## Resolve futures
  res <- as.list(res)

  invisible(res)
})



setMethodS3("getLog2Ratios", "CopyNumberSegmentationModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "getLog2Ratios()");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the regions for each of the fits (per array)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Obtaining fits (or fit if missing)");
  res <- suppressWarnings(local({
    ## Assume fit() has already been called;
    ## avoids void asynchroneous processing.
    oplan <- plan("eager")
    on.exit(plan(oplan))
    fit(this, ..., .retResults=TRUE, verbose=less(verbose,10))
  }))
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting regions from all fits");
  res <- lapply(res, FUN=function(arrayFits) {
    df <- NULL;
    # For each chromosome
    for (kk in seq_along(arrayFits)) {
      fit <- arrayFits[[kk]];
      if (!is.null(fit)) {
        verbose && enter(verbose, "Extracting regions for chromosome #", kk);
        suppressWarnings({
          df0 <- getRegions(fit, ...);
        })
        df <- rbind(df, df0);
        verbose && exit(verbose);
      }
    }
    rownames(df) <- seq_len(nrow(df));
    df;
  })
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}, private=TRUE) # getLog2Ratios()



setMethodS3("getRegions", "CopyNumberSegmentationModel", function(this, ..., url="ucsc", organism="Human", hgVersion="hg18", margin=0.1, flat=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'url':
  if (identical(url, "ucsc")) {
    # The UCSC browser accepts chromsomes either 'X' or 23.  In other words,
    # we can stick with integers to be more general.
    url <- function(chromosome, start, stop) {
      suppressWarnings({
        start <- as.double(start);
        stop <- as.double(stop);
      })
      uri <- "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=vertebrate";
      sprintf("%s&org=%s&db=%s&position=chr%s%%3A%.0f-%.0f", uri, organism, hgVersion, chromosome, start, stop);
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extracting regions from all fits");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the regions for each of the CN model fits (per array)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Obtaining CN model fits (or fit if missing)");
  res <- suppressWarnings(local({
    ## Assume fit() has already been called;
    ## avoids void asynchroneous processing.
    oplan <- plan("eager")
    on.exit(plan(oplan))
    fit(this, ..., .retResults=TRUE, verbose=less(verbose,10))
  }))
  verbose && exit(verbose);

  res <- lapply(res, FUN=function(arrayFits) {
    df <- NULL;
    # For each chromosome
    for (kk in seq_along(arrayFits)) {
      fit <- arrayFits[[kk]];
      if (!is.null(fit)) {
        verbose && enter(verbose, "Extracting regions for chromosome #", kk);
        suppressWarnings({
#          df0 <- getRegions(fit, ...);
          cnr <- extractCopyNumberRegions(fit, ...);
          df0 <- as.data.frame(cnr);
        })
        df <- rbind(df, df0);
        verbose && exit(verbose);
      }
    }
    rownames(df) <- seq_len(nrow(df));

    verbose && cat(verbose, "Extracted regions:");
    verbose && str(verbose, df);

    # Add URL?
    if (!is.null(url)) {
      chromosome <- df[,"chromosome"];
      start <- df[,"start"];
      stop <- df[,"stop"];
      m <- margin*abs(stop-start);
      start <- start-m;
      start[start < 0] <- 0;
      stop <- stop + m;
      urls <- character(nrow(df));
      for (rr in seq_along(urls)) {
        urls[rr] <- url(chromosome[rr], start[rr], stop[rr]);
      }
      df <- cbind(df, url=urls);
    }

    df;
  })
  verbose && exit(verbose);

  if (flat) {
    df <- NULL;
    for (kk in seq_along(res)) {
      df <- rbind(df, cbind(sample=names(res)[kk], res[[kk]]));
      res[[kk]] <- NA;
    }
    row.names(df) <- seq_len(nrow(df));
    res <- df;
  }

  res;
})


setMethodS3("writeRegions", "CopyNumberSegmentationModel", function(this, arrays=NULL, format=c("xls", "wig"), digits=3, ..., oneFile=TRUE, skip=TRUE, verbose=FALSE) {
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
  arrayNames <- getNames(this)[arrays];

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
    name <- arrayNames[aa];
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
        df <- df[,c("chromosome", "start", "stop", "log2")];

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
        chrIdx <- as.integer(df[,"chromosome"]);
        o <- order(chrIdx, df[,"start"]);
        df <- df[o,];

        # All chromosomes should have prefix 'chr'.
        chrIdx <- as.integer(df[,"chromosome"]);
        ## df[chrIdx == 23,"chromosome"] <- "X"; ## REMOVED 2007-03-15
        df[,"chromosome"] <- paste("chr", df[,"chromosome"], sep="");
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
      trackAttr <- c(trackAttr,
                     group=sprintf("\"%s regions\"",
                                getAsteriskTags(this, collapse=",")));
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
      write.table(df, file=pathname, sep="\t", col.names=FALSE,
                  row.names=FALSE, quote=FALSE, append=TRUE)
    }
    verbose && exit(verbose);
    res[[array]] <- df;
  }

  invisible(pathname);
})


setMethodS3("getFullNames", "CopyNumberSegmentationModel", function(this, ...) {
  # Get paired names
  names <- getNames(this, ...);

  # Get tags
  testTuple <- getSetTuple(this);
  fullnames <- getFullNames(testTuple);
  tags <- gsub("^[^,]*(|,)", "", fullnames);
  # Not needed anymore
  testTuple <- fullnames <- NULL;

  fullnames <- paste(names, tags, sep=",");
  fullnames <- gsub(",$", "", fullnames);

  fullnames;
})


##############################################################################
# HISTORY:
# 2011-02-07
# o Now fit() for CopyNumberSegmentationModel passes down argument
#   'maxNAFraction' to extractRawCopyNumbers().  This argument used to
#   work before aroma.core v1.3.4.
# 2010-10-25
# o Now fit() for CopyNumberSegmentationModel also passed the optional
#   arguments ('...') passed to the constructor function.  This makes it
#   possible to specify all arguments while initiating the model, e.g.
#   sm <- CbsModel(..., min.width=5, alpha=0.05).
# 2010-01-01
# o Added getFullNames() to CopyNumberSegmentationModel in order to be
#   backward compatible with previous versions.
# 2009-12-31
# o BUG FIX: After the recent updates, the getTags() methods of
#   CopyNumberSegmentationModel could give "Error in strsplit(tags,
#   split = ",") : non-character argument".
# 2009-11-22
# o CLEAN UP: Now extractRawCopyNumbers() is used; not old getRawCnData().
# 2009-11-16
# o Except from drop 'chipEffects' tags, the code of this class is completely
#   generic, that is, it does not assume Affymetrix data.  Note however,
#   that the code of super classes still assumes Affymetrix data.
# o CLEAN UP: Using getDataFileMatrix() instead of the old name
#   getMatrixChipEffectFiles().
# o CLEAN UP: Cleaning up code and comments so it is less specific to
#   Affymetrix data.
# 2009-05-16
# o Added generic getAsteriskTags() for CopyNumberSegmentationModel.
# o Classes extending CopyNumberSegmentationModel do no longer need to have
#   a fitOne() method.  Instead, they need to implement a static
#   getFitFunction() which should return a segmentByNnn() function for the
#   RawGenomicSignals class.
# 2008-03-10
# o Now verbose output of fit() reports the robust estimate of the log-ratios.
# 2008-02-27
# o BUG FIX: The URLs returned by getRegions() had broken positions.
# 2007-10-17
# o Most of the remain of this class was moved to superclass
#   CopyNumberChromosomalModel.
# o Remove methods that were moved up to the ChromosomalModel class.
# 2007-09-29
# o Now argument 'refTuple' works again.
# 2007-09-25
# o Now CopyNumberSegmentationModel extends ChromosomalModel.
# 2007-09-16
# o Now process() of CopyNumberSegmentationModel reports timing information
#   for each chromosome fitted.
# o Made a better job cleaning out non-needed objects in process().
# 2007-09-15
# o Added getGenome() and setGenome().  Now, getGenomeData() searches for a
#   matching <genome>(,<tags)*,chromosomes.txt in
#   annotationData/genomes/<genome>/. As a backup it also searches in the
#   corrsponding package directory(ies).
#   hgChromosomes.txt is no longer used.
# o BUG FIX: getGenomeData() was declared static, but it wasn't used that way.
# 2007-09-05
# o Was thinking to add a default asterisk tag to the output data set name.
#   However, although this works beautifully, the care has to be taken to
#   redesign ChromosomeExplorer, e.g. do you want to treat GLAD and CBS data
#   as totally different data sets?!? Possibly, because it is more consistent
#   with everything else, but for now we leave it as it since that works well.
# 2007-09-04
# o Now plot() is fully implemented CopyNumberSegmentationModel.  Subclasses
#   pretty much only have to implement pointsRawCNs() if wanted extra
#   features, e.g. colors and outliers.
# 2007-08-20
# o It was actually not much that was hardwired to the GLAD model.
#   Note that most methods, including fit(), are generic enough to be
#   defined in the superclass.  Methods that need to be implemented in
#   subclasses are: fitOne().
# o Created from GladModel.R.
##############################################################################
