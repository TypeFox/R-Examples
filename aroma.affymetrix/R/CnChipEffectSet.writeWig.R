setMethodS3("writeWig", "CnChipEffectSet", function(this, reference=NULL, arrays=1:length(this), chromosomes=c(1:22,"X"), na.rm=TRUE, oneFile=FALSE, gzip=TRUE, ylim=c(-1,1), digits=3, group="Copy numbers", col=NULL, smoothingWindow=5, verbose=FALSE, ...) {
  allChromosomes <- c(1:22,"X");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'reference':
  if (!is.null(reference)) {
    reference <- Arguments$getInstanceOf(reference, "CnChipEffects");

    if (reference$combineAlleles != this$combineAlleles) {
       throw("The reference chip effects are not compatible with the chip-effect set. One is combining the alleles the other is not.");
    }

    if (reference$mergeStrands != this$mergeStrands) {
       throw("The reference chip effects are not compatible with the chip-effect set. One is merging the strands the other is not.");
    }
  }

  # Argument 'arrays':
  arrays <- Arguments$getIndices(arrays, max=length(this));

  # Argument 'chromosomes':
  chromosomes <- Arguments$getCharacters(chromosomes);
  chromosomes <- intersect(chromosomes, allChromosomes);
  chrIdxs <- match(chromosomes, allChromosomes);
  chrsStr <- seqToHumanReadable(chrIdxs, collapse=";");

  # Argument 'oneFile':
  oneFile <- Arguments$getLogical(oneFile);

  # Argument 'gzip':
  gzip <- Arguments$getLogical(gzip);

  # Argument 'group':
  group <- Arguments$getCharacter(group, length=c(1,1));

  # Argument 'smoothingWindow':
  smoothingWindow <- Arguments$getInteger(smoothingWindow, range=c(0,16));

  # Argument 'col':
  if (is.null(col))
    col <- RColorBrewer::brewer.pal(5, "Dark2")[3:4];
  col <- col2rgb(col);
  col <- apply(col, MARGIN=2, FUN=paste, collapse=",");
  col <- rep(col, length.out=2);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  dataSetName <- getFullName(this);
  cdf <- getCdf(this);
  chipType <- getChipType(cdf, fullname=FALSE);
  arrayNames <- getNames(this);

  path <- filePath("glad", dataSetName, chipType);
  path <- Arguments$getWritablePath(path);

  pathnames <- c();
  con <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reference?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(reference)) {
    verbose && enter(verbose, "No reference specified. Calculating average chip effects");
    reference <- getAverageFile(this, verbose=less(verbose));
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Open output file?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (oneFile) {
    filename <- sprintf("%s,chr%s,log2.wig", getName(this), chrsStr);
    pathname <- filePath(path, filename);
    if (gzip) {
      pathname <- paste(pathname, ".gz", sep="");
      con <- gzfile(pathname, open="w", compression=9);
    } else {
      con <- file(pathname, open="w");
    }
    verbose && cat(verbose, "Pathname: ", pathname);
  }

  on.exit({
    if (!is.null(con)) {
      close(con);
      pathnames <- c(pathnames, pathname);
    }
  });

  res <- list();
  for (aa in arrays) {
    name <- arrayNames[aa];

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Open output file?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (!oneFile) {
      filename <- sprintf("%s,chr%s,log2.wig", name, chrsStr);
      pathname <- filePath(path, filename);
      if (gzip) {
        pathname <- paste(pathname, ".gz", sep="");
        con <- gzfile(pathname, open="w", compression=9);
      } else {
        con <- file(pathname, open="w");
      }
      verbose && cat(verbose, "Pathname: ", pathname);
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Write WIG track
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Writing data");

    # Write browser control
#     cat(file=con, sprintf("browser position chr%d\n", chrIdx));

    # Write track control
    cat(file=con, "track type=wiggle_0 ");
    cat(file=con, sprintf("name=\"%s\" ", name));
    cat(file=con, sprintf("group=\"%s\" ", group));
    cat(file=con, sprintf("priority=%d ", aa));
    cat(file=con, "graphType=bar visibility=full maxHeightPixels=128:96:64 ");
    cat(file=con, "yLineOnOff=on ");
    cat(file=con, sprintf("color=%s altColor=%s ", col[1], col[2]));
    if (!is.null(ylim))
      cat(file=con, sprintf("autoScale=off viewLimits=%.2f:%.2f ", ylim[1], ylim[2]));
    if (smoothingWindow > 0)
      cat(file=con, sprintf("windowingFunction=mean smoothingWindow=%d", smoothingWindow));
    cat(file=con, "\n");

    ce <- this[[aa]];

    for (chr in chromosomes) {
      verbose && enter(verbose,
                             sprintf("Array %s (#%d of %d) on chromosome %s",
                                             name, aa, length(arrays), chr));

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Get the chromosomal position and log-relative chip effects
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Extracting data");
      xm <- getXAM(ce, other=reference, chromosome=chr, verbose=less(verbose))[,c("x","M")];
      xm[,"M"] <- round(xm[,"M"], digits=digits);
      # Remove cases with missing values?
      if (na.rm) {
        xm <- xm[!is.na(xm[,"M"]),];
        xm <- xm[!is.na(xm[,"x"]),];
      }
      # Order by chromosomal position
      xm <- xm[order(xm[,"x"]),];
      verbose && exit(verbose);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Write data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      cat(file=con, sprintf("variableStep\tchrom=chr%s\n", chr));
      write.table(xm, file=con, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE);

      cat(file=con, "\n");

      verbose && exit(verbose);
    }

    # Flush and close file
    if (!oneFile) {
      close(con);
      con <- NULL;
      pathnames <- c(pathnames, pathname);
    }
  }

  invisible(pathnames);
}) # writeWig()



##############################################################################
# HISTORY:
# 2007-06-11
# o writeWig() of CnChipEffectSet used non-existing 'ces' (instead of 'this').
# 2006-11-27
# o BUG FIX: writeWig() would write NAs, but the UCSC Genome Browser does not
#   accept missing values.
# 2006-11-24
# o Added track arguments 'group' and 'priority' too. Now the arrays are
#   ordered in the UCSC Genome Browser as they are ordered in the data set.
# 2006-11-23
# o Added writeWig().
# o Created.
##############################################################################
