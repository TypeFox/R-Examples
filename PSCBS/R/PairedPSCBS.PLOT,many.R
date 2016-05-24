#   \item{chromosomes}{An optional @numeric @vector specifying which
#     chromosomes to plot.}
#
#   \item{seed}{An (optional) @integer specifying the random seed to be
#     set before subsampling.  The random seed is
#     set to its original state when exiting.  If @NULL, it is not set.}
#
#   \item{verbose}{See @see "R.utils::Verbose".}
#
setMethodS3("plotTracksManyChromosomes", "PairedPSCBS", function(fit, chromosomes=getChromosomes(fit), tracks=NULL, scatter="*", calls=if (callLoci || length(chromosomes) == 1L) ".*" else NULL, callLoci=FALSE, callThresholds=TRUE, boundaries=TRUE, knownSegments=FALSE, quantiles=c(0.05,0.95), seed=0xBEEF, pch=".", Clim=c(0,3*ploidy(fit)), Blim=c(0,1), xScale=1e-6, xlabTicks=if (length(chromosomes) == 1L) "[pos]" else "[chr]", ..., subset=if (length(chromosomes) > 1L) 0.1 else NULL, add=FALSE, subplots=!add && (length(tracks) > 1L), oma=c(0,0,2,0), mar=c(2,5,1,3)+0.1, onBegin=NULL, onEnd=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  attachGH <- function(gh, envir=parent.frame()) {
    if (!is.list(gh)) return();
    if (is.null(gh$track)) return();
    if (!is.null(value <- gh$track)) assign("track", value, envir=envir);
    if (!is.null(value <- gh$subtracks)) assign("trackT", value, envir=envir);
    if (!is.null(value <- gh$scatter$col)) assign("colS", value, envir=envir);
    if (!is.null(value <- gh$scatter$pch)) assign("pchT", value, envir=envir);
    if (!is.null(value <- gh$level$col)) assign("colL", value, envir=envir);
    if (!is.null(value <- gh$cis$col)) assign("colC", value, envir=envir);
  } # attachGH()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Graphical styles
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  opts <- list(
    scatter = list(pch=".", col=c("#aaaaaa")),
    callScatter = list(col=c("#aaaaaa", LOSS="blue", GAIN="red", LOH="purple")),
    smoothScatter = list(pch=".", col=c("#666666")),
    level = list(lty=1L, col=c("black", tcn="purple", c1="blue", c2="red", dh="orange")),
    callLevel = list(lty=1L, col=c("#666666")),
    knownSegment = list(lty=1L, col=c("#aaaaaa"))
  );

  getOptionValue <- function(option, what, track, ...) {
    values <- opts[[option]][[what]];
    value <- values[track];
    if (is.na(value)) value <- values[1L];
    unname(value);
  } # getOptionValue()

  getScatterColor <- function(track, ...) {
    getOptionValue("scatter", "col", track, ...);
  } # getScatterColor()

  getLevelColor <- function(track, ...) {
    getOptionValue("level", "col", track, ...);
  } # getLevelColor()

  getCallScatterColor <- function(track, ...) {
    getOptionValue("callScatter", "col", track, ...);
  } # getCallScatterColor()

  getCallLevelColor <- function(track, ...) {
    getOptionValue("callLevel", "col", track, ...);
  } # getCallLevelColor()

  getCallLevelLty <- function(track, ...) {
    getOptionValue("callLevel", "lty", track, ...);
  } # getCallLevelColor()

  getCIColor <- function(track, ...) {
    getLevelColor(track, ...);
  } # getLevelColor()

  getKnownSegmentColor <- function(track, ...) {
    getOptionValue("knownSegment", "col", track, ...);
  } # getKnownSegmentColor()

  getKnownSegmentLty <- function(track, ...) {
    getOptionValue("knownSegment", "lty", track, ...);
  } # getKnownSegmentColor()

  drawXLabelTicks <- function() {
    if (identical(xlabTicks, "[chr]")) {
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
    } else if (identical(xlabTicks, "[pos]")) {
      axis(side=1L);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fit':

  # Argument 'chromosomes':
  if (!is.null(chromosomes)) {
    disallow <- c("NaN", "Inf");
    chromosomes <- Arguments$getIntegers(chromosomes, range=c(0,Inf), disallow=disallow);
    stopifnot(is.element(chromosomes, getChromosomes(fit)));
  }

  # Argument 'tracks':
  knownTracks <- c("tcn", "dh", "tcn,c1,c2", "c1,c2", "c1", "c2",
                   "betaN", "betaT", "betaTN");
  defaultTracks <- knownTracks[1:3];
  if (is.null(tracks)) {
    tracks <- defaultTracks;
  } else {
    tracks <- match.arg(tracks, choices=knownTracks, several.ok=TRUE);
    tracks <- unique(tracks);
  }

  # Argument 'scatter':
  if (!is.null(scatter)) {
    scatter <- Arguments$getCharacter(scatter);
    if (scatter == "*") {
      scatter <- tracks;
    } else {
      scatterT <- strsplit(scatter, split=",", fixed=TRUE);
      tracksT <- strsplit(tracks, split=",", fixed=TRUE);
      stopifnot(all(is.element(scatterT, tracksT)));
      # Not needed anymore
      scatterT <- tracksT <- NULL;
    }
  }

  # Argument 'calls':
  if (!is.null(calls)) {
    calls <- sapply(calls, FUN=Arguments$getRegularExpression);
  }

  # Argument 'callLoci':
  callLoci <- Arguments$getLogical(callLoci);

  # Argument 'callThresholds':
  callThresholds <- Arguments$getLogical(callThresholds);

  # Argument 'boundaries':
  boundaries <- Arguments$getLogical(boundaries);

  # Argument 'knownSegments':
  knownSegments <- Arguments$getLogical(knownSegments);

  # Argument 'add':
  add <- Arguments$getLogical(add);

  # Argument 'Clim' & 'Blim':
  if (!add) {
    Clim <- Arguments$getNumerics(Clim, length=c(2L,2L),
                                        disallow=c("Inf", "NA", "NaN"));
    Blim <- Arguments$getNumerics(Blim, length=c(2L,2L),
                                        disallow=c("Inf", "NA", "NaN"));
  }

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf));

  # Argument 'xlabTicks':
  if (!is.null(xlabTicks)) {
    xlabTicks <- Arguments$getCharacter(xlabTicks);
  }

  # Argument 'subset':
  if (!is.null(subset)) {
    subset <- Arguments$getDouble(subset);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Plotting PSCN tracks");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset by chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(chromosomes)) {
    verbose && enter(verbose, "Plotting a subset of the chromosomes");
    fit <- extractChromosomes(fit, chromosomes=chromosomes, verbose=verbose);
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit <- tileChromosomes(fit, verbose=verbose);
  verbose && str(verbose, fit);

  # Extract the input data
  data <- getLocusData(fit);
  if (is.null(data)) {
    throw("Cannot plot segmentation results. No input data available.");
  }

  # Extract the segmentation
  segs <- as.data.frame(fit);

  # Identify available calls
  callData <- NULL;
  if (!is.null(calls) || callThresholds) {
    verbose && enter(verbose, "Identifying calls");

    pattern <- "Call$";
    allCallColumns <- grep(pattern, colnames(segs), value=TRUE);
    allCallLabels <- toupper(gsub(pattern, "", allCallColumns));
    verbose && cat(verbose, "Call columns:");
    verbose && print(verbose, allCallColumns);

    if (!is.null(calls)) {
      callColumns <- allCallColumns;
      if (length(callColumns) > 0L) {
        keep <- sapply(calls, FUN=function(pattern) {
          (regexpr(pattern, callColumns) != -1L);
        });
        if (is.matrix(keep)) {
          keep <- rowAnys(keep);
        }
        callColumns <- callColumns[keep];
        callLabels <- allCallLabels[keep];

        # Annotate individual loci by calls?
        if (callLoci) {
          callData <- extractCallsByLocus(fit, verbose=less(verbose,5));
        }
      }
      verbose && cat(verbose, "Call to be annotated:");
      verbose && print(verbose, callColumns);
    }

    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset of the loci?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(subset)) {
    # (a) Set and unset the random seed
    if (!is.null(seed)) {
      randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
      on.exit(randomSeed("reset"), add=TRUE)
      verbose && printf(verbose, "Random seed temporarily set (seed=%d)\n", seed)
    }

    # (b) Subset
    n <- nrow(data);
    keep <- sample(n, size=subset*n);
    data <- data[keep,];
    if (!is.null(callData)) {
      callData <- callData[keep,];
    }
  }

  # To please R CMD check
  CT <- rho <- muN <- betaT <- betaN <- betaTN <- rho <- NULL;
  rm(list=c("CT", "rho", "muN", "betaT", "betaN", "betaTN"));
  attachLocally(data);
  # Calculate (C1,C2)
  C1 <- 1/2*(1-rho)*CT;
  C2 <- CT - C1;

  # BACKWARD COMPATIBILITY:
  # If 'rho' is not available, recalculate it from tumor BAFs.
  # NOTE: This should throw an error in the future. /HB 2013-10-25
  if (is.null(data$rho)) {
    isSnp <- (!is.na(betaTN) & !is.na(muN));
    isHet <- isSnp & (muN == 1/2);
    rho <- rep(NA_real_, length=nbrOfLoci);
    rho[isHet] <- 2*abs(betaTN[isHet]-1/2);
    warning(sprintf("Locus-level DH signals ('rho') were not available in the %s object and therefore recalculated from the TumorBoost-normalized tumor BAFs ('betaTN').", class(fit)[1L]));
  }

  x <- xScale * x;
  vs <- xScale * fit$chromosomeStats[,1:2,drop=FALSE];
  mids <- (vs[,1]+vs[,2])/2;

  nbrOfLoci <- length(x);
  chrTags <- sprintf("Chr%02d", chromosomes);

  if (subplots) {
    subplots(length(tracks), ncol=1L);
    par(oma=oma, mar=mar);
  }

  pchT <- if (!is.null(scatter)) { pch } else { NA };
  xlim <- range(x, na.rm=TRUE);
  xlab <- "Genomic position";

  # Graphical handle
  gh <- list(fit=fit);
  gh$xScale <- xScale;
  gh$xlim <- xlim;
  gh$xlab <- xlab;
  if (!is.null(callData)) {
    gh$callsByLocus <- callData;
  }

  for (tt in seq(along=tracks)) {
    track <- tracks[tt];
    verbose && enter(verbose, sprintf("Track #%d ('%s') of %d",
                                             tt, track, length(tracks)));

    # Get graphical style parameters.
    tracksT <- unlist(strsplit(track, split=",", fixed=TRUE));
    colS <- sapply(tracksT, FUN=getScatterColor);
    colL <- sapply(tracksT, FUN=getLevelColor);
    colC <- sapply(tracksT, FUN=getCIColor);

    # Color scatter plot according to calls?
    if (!is.null(calls) && callLoci && length(callColumns) > 0L) {
      colsT <- rep(colS[1L], times=nrow(callData));
      for (cc in seq(along=callColumns)) {
        callColumn <- callColumns[cc];
        callLabel <- callLabels[cc];
        verbose && enter(verbose, sprintf("Call #%d ('%s') of %d",
                                      cc, callLabel, length(callColumns)));

        verbose && cat(verbose, "Column: ", callColumn);

        skip <- TRUE;
        if (regexpr("tcn", track) != -1L) {
          skip <- !is.element(callLabel, c("LOSS", "NTCN", "GAIN", "LOH"));
        } else if (track == "dh") {
          skip <- !is.element(callLabel, c("AB", "LOH"));
        }
        if (skip) {
          verbose && exit(verbose);
          next;
        }

        callsCC <- callData[[callColumn]];
        idxs <- which(callsCC);
        # Nothing to do?
        if (length(idxs) == 0L) {
          verbose && exit(verbose);
          next;
        }

        callCol <- getCallScatterColor(callLabel);

        colsT[idxs] <- callCol;
      } # for (cc in ...)

      colS <- colsT;
    } # if (!is.null(calls))


    # Assign graphical-handle parameters
    gh$track <- track;
    gh$subtracks <- tracksT;
    gh$scatter <- list(col=colS, pch=pchT);
    gh$level <- list(col=colL);
    gh$cis <- list(col=colC);


    if (track == "tcn") {
      plot(NA, xlim=xlim, ylim=Clim, xlab=xlab, ylab="TCN", axes=FALSE);
      if (!is.null(onBegin)) attachGH(onBegin(gh=gh));
      if (!is.na(pchT)) {
        points(x, CT, pch=pchT, col=colS);
      }
      drawConfidenceBands(fit, what="tcn", quantiles=quantiles, col=colC["tcn"], xScale=xScale);
      drawLevels(fit, what="tcn", col=colL, xScale=xScale);
    }

    if (is.element(track, c("tcn,c1,c2", "c1,c2", "c1", "c2"))) {
      tracksT <- unlist(strsplit(track, split=",", fixed=TRUE));
      plot(NA, xlim=xlim, ylim=Clim, xlab=xlab, ylab="C1, C2, TCN", axes=FALSE);
      if (!is.null(onBegin)) attachGH(onBegin(gh=gh));

      # Draw scatter for TCN or C1 and C2.
      if (!is.na(pchT)) {
        if (is.element("tcn", tracksT)) {
          points(x, CT, pch=pchT, col=colS);
        } else {
          if (is.element("c1", tracksT)) {
            points(x, C1, pch=pchT, col=colS);
          }
          if (is.element("c2", tracksT)) {
            points(x, C2, pch=pchT, col=colS);
          }
        }
      }

      # Draw confidence bands for TCN, C1, C2.
      if (is.element("tcn", tracksT)) {
        drawConfidenceBands(fit, what="tcn", quantiles=quantiles, col=colC["tcn"], xScale=xScale);
      }
      if (is.element("c2", tracksT)) {
        drawConfidenceBands(fit, what="c2", quantiles=quantiles, col=colC["c2"], xScale=xScale);
      }
      if (is.element("c1", tracksT)) {
        drawConfidenceBands(fit, what="c1", quantiles=quantiles, col=colC["c1"], xScale=xScale);
      }

      # Draw segment means for TCN, C1, C2.
      if (is.element("tcn", tracksT)) {
        drawLevels(fit, what="tcn", col=colL["tcn"], xScale=xScale);
      }
      if (is.element("c2", tracksT)) {
        drawLevels(fit, what="c2", col=colL["c2"], xScale=xScale);
      }
      if (is.element("tcn", tracksT)) {
        # In case C2 overlaps with TCN
        drawLevels(fit, what="tcn", col=colL["tcn"], lty="22", xScale=xScale);
      }
      # In case C1 overlaps with C2
      if (is.element("c1", tracksT)) {
        drawLevels(fit, what="c1", col=colL["c1"], xScale=xScale);
        if (is.element("c2", tracksT)) {
          drawLevels(fit, what="c2", col=colL["c2"], lty="22", xScale=xScale);
        }
        if (is.element("tcn", tracksT)) {
          drawLevels(fit, what="tcn", col=colL["tcn"], lty="22", xScale=xScale);
        }
      }
    }

    if (track == "betaN") {
      plot(NA, xlim=xlim, ylim=Blim, xlab=xlab, ylab="BAF_N", axes=FALSE);
      if (!is.null(onBegin)) attachGH(onBegin(gh=gh));
      if (!is.na(pchT)) {
        points(x, betaN, pch=pchT, col=colS);
      }
    }

    if (track == "betaT") {
      plot(NA, xlim=xlim, ylim=Blim, xlab=xlab, ylab="BAF_T", axes=FALSE);
      if (!is.null(onBegin)) attachGH(onBegin(gh=gh));
      if (!is.na(pchT)) {
        points(x, betaT, pch=pchT, col=colS);
      }
    }

    if (track == "betaTN") {
      plot(NA, xlim=xlim, ylim=Blim, xlab=xlab, ylab="BAF_TN", axes=FALSE);
      if (!is.null(onBegin)) attachGH(onBegin(gh=gh));
      if (!is.na(pchT)) {
        points(x, betaTN, pch=pchT, col=colS);
      }
    }

    if (track == "dh") {
      plot(NA, xlim=xlim, ylim=Blim, xlab=xlab, ylab="DH", axes=FALSE);
      if (!is.null(onBegin)) attachGH(onBegin(gh=gh));
      if (!is.na(pchT)) {
        points(x, rho, pch=pchT, col=colS);
      }
      drawConfidenceBands(fit, what="dh", quantiles=quantiles, col=colC["dh"], xScale=xScale);
      drawLevels(fit, what="dh", col=colL["dh"], xScale=xScale);
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # For each panel of tracks, annotate segments with calls?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (!is.null(calls) && !callLoci && length(callColumns) > 0L) {
      for (cc in seq(along=callColumns)) {
        callColumn <- callColumns[cc];
        callLabel <- callLabels[cc];
        verbose && enter(verbose, sprintf("Call #%d ('%s') of %d",
                                      cc, callLabel, length(callColumns)));

        verbose && cat(verbose, "Column: ", callColumn);

        segsT <- segs[,c("dhStart", "dhEnd", callColumn)];
        isCalled <- which(segsT[[callColumn]]);
        segsT <- segsT[isCalled,1:2,drop=FALSE];
        verbose && printf(verbose, "Number of segments called %s: %d\n",
                                                  callLabel, nrow(segsT));
        segsT <- xScale * segsT;

        verbose && str(verbose, segsT);

        side <- 2*((cc+1) %% 2) + 1;
        # For each segment called...
        for (ss in seq(length=nrow(segsT))) {
          x0 <- segsT[ss,1,drop=TRUE];
          x1 <- segsT[ss,2,drop=TRUE];
          abline(v=c(x0,x1), lty=3, col="gray");
          xMid <- (x0+x1)/2;
          mtext(side=side, at=xMid, line=-1, cex=0.7, col="#666666", callLabel);
        } # for (ss in ...)
        verbose && exit(verbose);
      } # for (cc in ...)
    } # if (!is.null(calls))


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # For each panel of tracks, annotate with call thresholds?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (callThresholds) {
      # Add call parameter estimates, e.g. deltaAB
      colCL <- sapply(tracksT, FUN=getCallLevelColor);
      ltyCL <- sapply(tracksT, FUN=getCallLevelLty);

      trackT <- track;
      for (cc in seq(along=allCallColumns)) {
        callColumn <- allCallColumns[cc];
        callLabel <- allCallLabels[cc];

        h <- NULL;
        if (callLabel == "AB") {
          if (track == "dh") {
            h <- fit$params$deltaAB;
            label <- expression(Delta[AB]);
          }
        } else if (callLabel == "LOH") {
          if (regexpr("c1", track) != -1L) {
            h <- fit$params$deltaLowC1;
            label <- expression(Delta[LOH]);
            trackT <- "c1";
          }
        } else if (callLabel == "NTCN") {
          if (track == "tcn") {
            h <- fit$params$ntcnRange;
            label <- c(expression(Delta[-NTCN]), expression(Delta[+NTCN]));
          }
        }

        if (!is.null(h)) {
          abline(h=h, lty=ltyCL[trackT], lwd=2, col=colCL[trackT]);
          for (ss in 1:2) {
            side <- c(2,4)[ss];
            adj <- c(1.2,-0.2)[ss];
            mtext(side=side, at=h, label, adj=adj, las=2, xpd=TRUE);
          }
        }
      } # for (cc in ...)
    } # if (callThresholds)

    drawXLabelTicks();
    if (boundaries) {
      abline(v=vs, lty=1, lwd=2);
    }
    if (knownSegments) {
      colT <- getKnownSegmentColor();
      ltyT <- getKnownSegmentLty();
      drawKnownSegments(fit, col=colT, lty=ltyT);
    }
    axis(side=2); box();
    if (!is.null(onEnd)) onEnd(gh=gh);

    verbose && exit(verbose);
  } # for (tt ...)

  verbose && exit(verbose);

  invisible(gh);
}, private=TRUE) # plotTracksManyChromosomes()



############################################################################
# HISTORY:
# 2013-10-28
# o Now plotTracksManyChromosomes() also supports
#   tracks=c("c1,c2", "c1", "c2").
# 2013-10-25
# o Now plotTracksManyChromosomes() uses the locus data field 'rho'
#   when plotting DH locus-level data.  It only recalculates it from
#   the tumor BAFs if the DH signals are not available - if so a
#   warning is generated.
# 2013-10-20
# o BUG FIX: plotTracksManyChromosomes() for PairedPSCBS would use
#   Blim=Clim, regardless of what argument 'Blim' was.
# 2013-04-13
# o Added argument 'boundaries' to plotTracksManyChromosomes().
# 2013-04-11
# o BUG FIX: plotTracksManyChromosomes(..., callLoci=TRUE) would color
#   loci incorrectly if more than one chromosome are plotted.
# 2013-04-05
# o Now plotTracks() passes more information to onBegin(gh)/onEnd(gh)
#   hooks via the graphical handle object, cf. str(gh).
# 2013-03-21
# o Added argument 'knownSegments' to plotTracksManyChromosomes().
# o Generalized plotTracksManyChromosomes() a bit such that it can be
#   used for a single chromosome as well.  All col and lty annotations
#   are now specified at the very top of the function.
# 2013-03-18
# o Now plotTracksManyChromosomes() draws AB and LOH call thresholds.
# 2012-09-23
# o Now plotTracks() [and plotTracksManyChromosomes()] draws segment levels
#   in TCN-C2-C1 order, and then goes back and draws C2 and TCN with dashed
#   lines, just in case C1 is overlapping C2 and C2 is overlapping TCN.
# 2012-09-21
# o ROBUSTNESS: Now drawChangePointsC1C2() and arrowsC1C2() for PairedPSCBS
#   makes sure to retrieve segments with NA splitters between chromosomes
#   and gaps.
# 2012-02-29
# o BUG FIX: plotTracks(..., add=TRUE) for PairedPSCBS would add TCNs
#   when BAFs and DHs where intended.
# 2012-02-22
# o BUG FIX: Argument 'calls' of plotTracks() for PairedPSCBS was ignored
#   if more than one chromosome was plotted.
# 2011-12-03
# o Added drawChangePointsC1C2() for PairedPSCBS.
# o Added drawChangePoints() for PairedPSCBS.
# 2011-11-12
# o Added argument col="#00000033" to plotC1C2() and linesC1C2().
# o Added argument 'oma' and 'mar' to plotTracksManyChromosomes() for
#   PairedPSCBS for setting graphical parameters when 'add' == FALSE.
# 2011-09-30
# o GENERALIZATION: Now drawLevels() for PairedPSCBS allows for drawing
#   segmentation results in 'betaT' space.
# 2011-07-10
# o BUG FIX: tileChromosomes() for PairedPSCBS was still assuming the
#   old naming convention of column names.
# o ROBUSTNESS: Fixed partial argument matchings in arrowsC1C2() and
#   arrowsDeltaC1C2() for PairedPSCBS.
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-01-19
# o Added argument 'subplots'.
# 2011-01-18
# o DOCUMENTATION: Documented more plotTracks() arguments for PairedPSCBS.
# o Now plotTracks(..., add=TRUE) for PairedPSCBS plots to the current
#   figure/panel.
# o Now plotTracks(..., add=FALSE) for PairedPSCBS only sets up subplots
#   if argument 'tracks' specifies more than one panel.
# o Added argument 'col' to plotTracks() for PairedPSCBS.
# 2011-01-12
# o Added argument 'changepoints' to plotTracks() for PairedPSCBS for
#   highlightning change-point locations as vertical lines.
# 2010-12-01
# o Now using a default 'seed' for plotTracksManyChromosomes() such
#   that the scatter data in the plots are reproducible by default.
# 2010-11-27
# o BUG FIX: plotTracksManyChromosomes() would annotate called regions
#   incorrectly.
# o Added more verbouse output to plotTracksManyChromosomes().
# o Added missing argument 'verbose' to plotTracksManyChromosomes().
# o plotTracksManyChromosomes() gained argument 'scatter'.
# o REPRODUCIBILITY: plotTracksManyChromosomes() for PairedPSCBS gained
#   argument 'seed', because if 'subset' is specified then a random
#   subset of the data points are displayed.
# 2010-11-26
# o Added optional argument 'chromosomes' to plotTracks() to plot a
#   subset of all chromosomes.
# o Now the default confidence intervals for plotTracks() is (0.05,0.95),
#   if existing.
# 2010-11-22
# o ROBUSTNESS: Now drawConfidenceBands() of PairedPSCBS silently does
#   nothing if the requested bootstrap quantiles are available.
# o Added argument 'calls' to plotTracks() and plotTracksManyChromosomes()
#   for highlighing called regions.
# 2010-11-21
# o Now plotTracks() supports tracks "tcn,c1", "tcn,c2" and "c1,c2" too.
# o Added argument 'xlim' to plotTracks() making it possible to zoom in.
# o Now plotTracks() and plotTracksManyChromosomes() draws confidence
#   bands, iff argument quantiles is given.
# o Added drawConfidenceBands() for PairedPSCBS.
# 2010-11-09
# o Added argument 'cex=1' to plotTracks().
# o BUG FIX: It was not possible to plot BAF tracks with plotTracks().
# 2010-10-20
# o Added arguments 'onBegin' and 'onEnd' to plotTracksManyChromosomes().
# 2010-10-18
# o Now plotTracks() can plot whole-genome data.
# o Now plotTracks() utilized plotTracksManyChromosomes() if there is
#   more than one chromosome.
# o Added internal plotTracksManyChromosomes().
# o Added internal tileChromosomes().
# 2010-10-03
# o Now the default is that plotTracks() for PairedPSCBS generated three
#   panels: (1) TCN, (2) DH, and (3) C1+C2+TCN.
# o Added plotTracks() to be more explicit than just plot().
# o Added argument 'xScale' to plot() for PairedPSCBS.
# o Now plot() for PairedPSCBS adds a chromosome tag.
# 2010-09-21
# o Added argument 'what' to plot() for PairedPSCBS.
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
