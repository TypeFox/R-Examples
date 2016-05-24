setMethodS3("shiftTCN", "PairedPSCBS", function(fit, shift, update=TRUE, ...) {
  # Argument 'shift':
  shift <- Arguments$getDouble(shift, disallow=c("NA", "NaN", "Inf"));

  data <- getLocusData(fit);
  data$CT <- data$CT + shift;
  fit$data <- data;
  # Not needed anymore
  data <- NULL;

  if (update) {
    fit <- updateMeans(fit, ...);
  }

  fit;
}, protected=TRUE)


setMethodS3("bootstrapCIs", "PairedPSCBS", function(fit, ...) {
  # ...
}, private=TRUE)


###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod extractTCNAndDHs
#
# @title "Extract TCN and DH mean levels per segment"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to \code{getSegments()}.}
# }
#
# \value{
#   Returns a @data.frame.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "extractMinorMajorCNs".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractTCNAndDHs", "PairedPSCBS", function(fit, ...) {
  segs <- getSegments(fit, ...);
  stopifnot(!is.null(segs));

  data <- segs[,c("tcnMean", "dhMean", "tcnNbrOfLoci", "dhNbrOfLoci"), drop=FALSE];
  data;
}, protected=TRUE)



###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod extractMinorMajorCNs
# @aliasmethod extractC1C2
#
# @title "Extract minor and major copy-number mean levels per segment"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @data.frame.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "extractTCNAndDHs"
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractMinorMajorCNs", "PairedPSCBS", function(fit, ...) {
  data <- extractTCNAndDHs(fit, ...);

  gamma <- data[,1L];
  rho <- data[,2L];
  C1 <- 1/2*(1-rho)*gamma;
  C2 <- gamma - C1;

  data[,1L] <- C1;
  data[,2L] <- C2;
  colnames(data)[1:2] <- c("C1", "C2");

  # Swap (C1,C2)?
  segs <- getSegments(fit, ...);
  flipped <- segs$c1c2Swap;
  if (!is.null(flipped)) {
    idxs <- which(flipped);
    if (length(idxs) > 0L) {
      data[idxs,1:2] <- data[idxs,2:1];
    }
  }

  data;
}, protected=TRUE)


setMethodS3("extractC1C2", "PairedPSCBS", function(...) {
  extractMinorMajorCNs(...);
}, protected=TRUE)


setMethodS3("extractCNs", "PairedPSCBS", function(fit, splitters=TRUE, ...) {
  data <- extractC1C2(fit, splitters=splitters, ...);
  data[,c("C1", "C2"), drop=FALSE];
})


setMethodS3("extractDeltaC1C2", "PairedPSCBS", function(...) {
  xy <- extractC1C2(..., splitters=TRUE, addGaps=TRUE);
  X <- xy[,1:2,drop=FALSE];
  dX <- colDiffs(X);
  dX;
}, protected=TRUE)



setMethodS3("postsegmentTCN", "PairedPSCBS", function(fit, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Post-segmenting TCNs");

  flavor <- fit$params$flavor;
  if (!force && regexpr("&", flavor, fixed=TRUE) != -1) {
    verbose && cat(verbose, "Nothing to do. Already postsegmentTCN:ed: ", flavor);
    verbose && exit(verbose);
    return(fit);
  }

  joinSegments <- fit$params$joinSegments;
  if (!joinSegments) {
    throw("Postsegmentation of TCNs is only implemented for the case when joinSegments=TRUE: ", joinSegments);
  }


  # Get mean estimators
  estList <- getMeanEstimators(fit, "tcn");
  avgTCN <- estList$tcn;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit);

  segs <- getSegments(fit);
  keep <- is.finite(segs$chromosome);
  segs <- segs[keep,,drop=FALSE];
  tcnSegRows <- fit$tcnSegRows[keep,,drop=FALSE];
  dhSegRows <- fit$dhSegRows[keep,,drop=FALSE];

  # Sanity check
  stopifnot(nrow(dhSegRows) == nrow(tcnSegRows));
  stopifnot(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE));
#  stopifnot(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE));
  stopifnot(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE));
  stopifnot(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE));


  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  chromosome <- data$chromosome;
  x <- data$x;
  CT <- data$CT;
  muN <- data$muN;
  rho <- data$rho;
  hasDH <- !is.null(rho)
  if (hasDH) {
    isHet <- !is.na(rho)
    isSnp <- isHet
  } else {
    isSnp <- !is.na(muN)
    isHet <- (isSnp & (muN == 1/2))
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update the TCN segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chromosomes <- getChromosomes(fit);
  nbrOfChromosomes <- length(chromosomes);
  verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes);
  verbose && print(verbose, chromosomes);

  for (cc in seq(length=nbrOfChromosomes)) {
    chr <- chromosomes[cc];
    chrTag <- sprintf("chr%02d", chr);
    verbose && enter(verbose, sprintf("Chromosome %d ('%s') of %d", cc, chrTag, nbrOfChromosomes));
    rows <- which(is.element(segs[["chromosome"]], chr));
    verbose && cat(verbose, "Rows:");
    verbose && print(verbose, rows);

    segsCC <- segs[rows,,drop=FALSE];
    tcnSegRowsCC <- tcnSegRows[rows,,drop=FALSE];
    dhSegRowsCC <- dhSegRows[rows,,drop=FALSE];
    nbrOfSegmentsCC <- nrow(segsCC);
    verbose && cat(verbose, "Number of segments: ", nbrOfSegmentsCC);

    tcnIds <- sort(unique(segsCC[["tcnId"]]));
    I <- length(tcnIds);
    for (ii in seq(length=I)) {
      tcnId <- tcnIds[ii];
      verbose && enter(verbose, sprintf("TCN segment #%d ('%s') of %d", ii, tcnId, I));

      rowsII <- which(segsCC[["tcnId"]] == tcnId);
      J <- length(rowsII);
      # Nothing todo?
      if (!force && J == 1) {
        verbose && cat(verbose, "Nothing todo. Only one DH segmentation. Skipping.");
        verbose && exit(verbose);
        next;
      }

      verbose && cat(verbose, "Rows:");
      verbose && print(verbose, rowsII);
      segsII <- segsCC[rowsII,,drop=FALSE];

      tcnSegRowsII <- tcnSegRowsCC[rowsII,,drop=FALSE];
      dhSegRowsII <- dhSegRowsCC[rowsII,,drop=FALSE];

      verbose && cat(verbose, "TCN & DH segRows before:");
      verbose && print(verbose, cbind(tcn=tcnSegRowsII, dh=dhSegRowsII));

      segRowsRange <- range(c(tcnSegRowsII, dhSegRowsII), na.rm=TRUE);
      verbose && printf(verbose, "Range [%d,%d]\n",
                                    segRowsRange[1], segRowsRange[2]);

      tcnSegRowsIIBefore <- tcnSegRowsII;
      nbrOfTCNsBefore <- segsII[1,"tcnNbrOfLoci"];
      # Sanity check
      stopifnot(diff(segRowsRange)+1L == nbrOfTCNsBefore);

      for (jj in seq(length=J)) {
        verbose && enter(verbose, sprintf("DH segment #%d of %d", jj, J));
        seg <- segsII[jj,,drop=FALSE];
        tcnSegRow <- unlist(tcnSegRowsII[jj,,drop=FALSE], use.names=FALSE);
        dhSegRow <- unlist(dhSegRowsII[jj,,drop=FALSE], use.names=FALSE);
        # Sanity check
        stopifnot(all(is.na(tcnSegRow)) || (tcnSegRow[1] <= tcnSegRow[2]));
        stopifnot(all(is.na(dhSegRow)) || (dhSegRow[1] <= dhSegRow[2]));

        # Sanity check
        idxsTCN <- tcnSegRow[1]:tcnSegRow[2];
        nbrOfTCNs <- sum(!is.na(CT[idxsTCN]));
        stopifnot(nbrOfTCNs == nbrOfTCNsBefore);

        if (joinSegments) {
          # (a) The TCN segment should have identical (start,end) boundaries as the DH region
          xStart <- seg[["dhStart"]];
          xEnd <- seg[["dhEnd"]];
          verbose && printf(verbose, "[xStart,xEnd] = [%.0f,%.0f]\n", xStart, xEnd);
          stopifnot(xStart <= xEnd);

          # (b) Identify units
          units <- which(chromosome == chr & xStart <= x & x <= xEnd);

          # (c) Drop units that are outside both the TCN and DH segments
          keep <- (segRowsRange[1] <= units & units <= segRowsRange[2]);
          units <- units[keep];

          tcnSegRow <- range(units);
          verbose && printf(verbose, "[idxStart,idxEnd] = [%d,%d]\n", tcnSegRow[1], tcnSegRow[2]);
          verbose && cat(verbose, "Number of TCN loci: ", length(units));

          # (c) Adjust for missing values
          keep <- which(!is.na(CT[units]));
          units <- units[keep];

          # (d) Adjust for DH boundaries
          if (jj > 1L) {
            minIdx <- tcnSegRowsII[jj-1L,2L, drop=TRUE];
            units <- units[units > minIdx];
          }
          if (jj < J) {
            maxIdx <- dhSegRowsII[jj+1L,1L, drop=TRUE];
            units <- units[units < maxIdx];
          }

          if (jj == J) {
#           maxIdx <- dhSegRowsII[jj+1L,1L, drop=TRUE];
#           units <- units[units < maxIdx];
          }

          tcnSegRow <- range(units);
          verbose && printf(verbose, "[idxStart,idxEnd] = [%d,%d]\n", tcnSegRow[1], tcnSegRow[2]);
          verbose && cat(verbose, "Number of non-missing TCN loci: ", length(units));
        } else {
          throw("Not implemented yet.")  # /HB 2010-12-02
        } # if (joinSegments)

        gamma <- avgTCN(CT[units]);
        # Sanity check
        stopifnot(length(units) == 0 || !is.na(gamma));

        # Update the segment boundaries, estimates and counts
        seg[["tcnStart"]] <- xStart;
        seg[["tcnEnd"]] <- xEnd;
        seg[["tcnMean"]] <- gamma;
        seg[["tcnNbrOfLoci"]] <- length(units);
        seg[["tcnNbrOfSNPs"]] <- sum(isSnp[units]);
        seg[["tcnNbrOfHets"]] <- sum(isHet[units]);

        # Sanity check
        stopifnot(nrow(seg) == length(jj));

        segsII[jj,] <- seg;
        tcnSegRowsII[jj,] <- tcnSegRow;

        verbose && exit(verbose);
      } # for (jj ...)

      # Sanity check
      stopifnot(nrow(segsII) == length(rowsII));

      verbose && cat(verbose, "TCN & DH segRows afterward:");
      verbose && print(verbose, cbind(tcn=tcnSegRowsII, dh=dhSegRowsII));

##print(segsII);

      # Sanity check
      nbrOfTCNsAfter <- sum(segsII[,"tcnNbrOfLoci"], na.rm=TRUE);
      verbose && cat(verbose, "Number of TCNs before: ", nbrOfTCNsBefore);
      verbose && cat(verbose, "Number of TCNs after: ", nbrOfTCNsAfter);
      stopifnot(nbrOfTCNsAfter >= nbrOfTCNsBefore);

      # Sanity check
      stopifnot(nrow(dhSegRowsII) == nrow(tcnSegRowsII));
      stopifnot(all(tcnSegRowsII[,1] <= tcnSegRowsII[,2], na.rm=TRUE));
      stopifnot(all(tcnSegRowsII[-nrow(tcnSegRowsII),2] < tcnSegRowsII[-1,1], na.rm=TRUE));
      stopifnot(all(dhSegRowsII[,1] <= dhSegRowsII[,2], na.rm=TRUE));
      stopifnot(all(dhSegRowsII[-nrow(dhSegRowsII),2] < dhSegRowsII[-1,1], na.rm=TRUE));

      segsCC[rowsII,] <- segsII;
      tcnSegRowsCC[rowsII,] <- tcnSegRowsII;

      # Not needed anymore
      rowsII <- segsII <- NULL;
      verbose && exit(verbose);
    } # for (ii ...)

    # Sanity check
    stopifnot(nrow(segsCC) == length(rows));

    # Sanity check
    stopifnot(nrow(dhSegRowsCC) == nrow(tcnSegRowsCC));
    stopifnot(all(tcnSegRowsCC[,1] <= tcnSegRowsCC[,2], na.rm=TRUE));
####################
if (!all(tcnSegRowsCC[-nrow(tcnSegRowsCC),2] < tcnSegRowsCC[-1,1], na.rm=TRUE)) {

  aa <- tcnSegRowsCC[-nrow(tcnSegRowsCC),2];
  bb <- tcnSegRowsCC[-1,1];
  delta <- bb - aa;
  dd <- cbind(aa, bb, delta=delta);
  print(dd);
  dd <- subset(dd, delta == 0);
  print(dd);
  row <- dd[,1L,drop=TRUE];
  print(row);
  rr <- row + -10:10;
  dd <- data[rr,];
  rownames(dd) <- rr;
  print(dd);
print(tcnSegRowsII);
}
####################

    stopifnot(all(tcnSegRowsCC[-nrow(tcnSegRowsCC),2] < tcnSegRowsCC[-1,1], na.rm=TRUE));
    stopifnot(all(dhSegRowsCC[,1] <= dhSegRowsCC[,2], na.rm=TRUE));
    stopifnot(all(dhSegRowsCC[-nrow(dhSegRowsCC),2] < dhSegRowsCC[-1,1], na.rm=TRUE));

    segs[rows,] <- segsCC;
    tcnSegRows[rows,] <- tcnSegRowsCC;

    # Not needed anymore
    rows <- segsCC <- NULL;
    verbose && exit(verbose);
  } # for (cc ...)

  # Sanity check
  stopifnot(nrow(dhSegRows) == nrow(tcnSegRows));
  stopifnot(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE));
  stopifnot(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE));
  stopifnot(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE));
  stopifnot(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE));

  verbose && enter(verbose, "Update (C1,C2) per segment");
  # Append (C1,C2) estimates
  tcn <- segs$tcnMean;
  dh <- segs$dhMean;
  C1 <- 1/2*(1-dh)*tcn;
  C2 <- tcn - C1;
  segs$c1Mean <- C1;
  segs$c2Mean <- C2;
  verbose && exit(verbose);

  # Return results
  keep <- which(is.finite(fit$output$chromosome));
  fitS <- fit;
  fitS$data <- data;
  fitS$output[keep,] <- segs;
  fitS$tcnSegRows[keep,] <- tcnSegRows;

  # Sanity check
  tcnSegRows <- fitS$tcnSegRows;
  dhSegRows <- fitS$dhSegRows;
  stopifnot(nrow(dhSegRows) == nrow(tcnSegRows));
  stopifnot(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE));
  stopifnot(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE));
  stopifnot(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE));
  stopifnot(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE));

  # Update 'flavor'
  fitS$params$flavor <- gsub(",", "&", flavor, fixed=TRUE);

  verbose && exit(verbose);

  fitS;
}, protected=TRUE) # postsegmentTCN()




############################################################################
# HISTORY:
# 2013-08-15
# o Made extractMinorMajorCNs() for PairedPSCBS acknowledge the
#   'c1c2Swap' field.
# 2013-01-15
# o Now postsegmentTCN() uses the params$avgTCN estimator, iff given.
# 2012-09-21
# o ROBUSTNESS: Now extractDeltaC1C2() for PairedPSCBS makes sure to
#   retrieve segments with NA splitters between chromosomes and gaps.
# 2012-09-13
# o Added shiftTCN() for PairedPSCBS.
# 2012-01-21
# o CLEANUP: Removed left-over debug output in postsegmentTCN().
# 2012-01-09
# o Minor correction of a verbose message in postsegmentTCN().
# 2011-10-16
# o Added extractCNs().
# 2011-10-14
# o Now extractTCNAndDHs() passes '...' to getSegments().
# 2011-10-02
# o CLEANUP: Dropped empty callSegments() for PairedPSCBS.
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-04-08
# o BUG FIX: postsegmentTCN() for PairedPSCBS could generate an invalid
#   'tcnSegRows' matrix, where the indices for two consecutive segments
#   would overlap, which is invalid.
# 2011-04-05
# o BUG FIX: estimateHighDHQuantileAtAB() for PairedPSCBS would throw
#   an error on an undefined 'trim' if verbose output was used.
# 2011-02-17
# o Added arguments 'robust' and 'trim' to estimateMeanForDH().
# 2011-02-03
# o Added argument 'tauTCN' to estimateMeanForDH().
# 2011-01-27
# o Added flavor="DHskew" to estimateTauAB().
# o Added flavor="DH" to estimateTauAB() to estimate from DH instead
#   of hBAF.  As argued by the equations in the comments, these two
#   approaches gives virtually the same results.  The advantage with the
#   DH approach is that it requires one less degree of freedom.
# o Added estimateMeanForDH().
# 2011-01-18
# o BUG FIX: 'tcnSegRows' and 'dhSegRows' where not updated by
#   extractByRegions() for PairedPSCBS.
# 2011-01-14
# o Added estimateTauAB() for estimating the DeltaAB parameter.
# o Added estimateStdDevForHeterozygousBAF() for PairedPSCBS.
# o BUG FIX: extractByRegions() did not handle the case where multiple loci
#   at the same position are split up in two different segments.
# 2011-01-12
# o Added extractByRegions() and extractByRegion() for PairedPSCBS.
# o Now postsegmentTCN(..., force=TRUE) for PairedPSCBS also updates
#   the TCN estimates even for segments where the DH segmentation did
#   not find any additional change points.
# 2010-12-02
# o Now postsegmentTCN() assert that total number of TCN loci before
#   and after is the same.
# o Now postsegmentTCN() assert that joinSegment is TRUE.
# 2010-12-01
# o Now postsegmentTCN() checks if it is already postsegmented.
# 2010-11-30
# o TODO: postsegmentTCN() does not make sure of 'dhLociToExclude'. Why?
# o Now postsegmentTCN() recognizes the new 'tcnLociToExclude'.
# 2010-11-28
# o BUG FIX: postsegmentTCN() did not handle loci with the same positions
#   and that are split in two different segments.  It also did not exclude
#   loci with missing values.
# 2010-11-21
# o Adjusted postsegmentTCN() such that the updated TCN segment boundaries
#   are the maximum of the DH segment and the support by the loci.  This
#   means that postsegmentTCN() will work as expected both when signals
#   where segmented with 'joinSegments' being TRUE or FALSE.
# 2010-10-25
# o Now subsetByDhSegments() for PairedPSCBS handles the rare case when
#   markers with the same positions are split in two different segments.
# o Renamed subsetBySegments() for PairedPSCBS to subsetByDhSegments().
# 2010-09-26
# o Now subsetBySegments() for PairedPSCBS handles multiple chromosomes.
# o Now postsegmentTCN() PairedPSCBS handles multiple chromosomes.
# 2010-09-21
# o Added postsegmentTCN() for PairedPSCBS.
# 2010-09-19
# o BUG FIX: plot() used non-defined nbrOfLoci; now length(x).
# 2010-09-15
# o Added subsetBySegments().
# o Added linesC1C2() and arrowsC1C2().
# o Now the default 'cex' for pointsC1C2() corresponds to 'dh.num.mark'.
# o Now extractTotalAndDH() also returns 'dh.num.mark'.
# 2010-09-08
# o Added argument 'add=FALSE' to plot().
# o Added plotC1C2().
# o Added extractTotalAndDH() and extractMinorMajorCNs().
# 2010-09-04
# o Added drawLevels() for PairedPSCBS.
# o Added as.data.frame() and print() for PairedPSCBS.
# 2010-09-03
# o Added plot() for PairedPSCBS.
# o Created.
############################################################################
