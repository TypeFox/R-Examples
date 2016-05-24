setMethodS3("plotTracks2", "PairedPSCBS", function(x, panels=NULL, calls=".*", pch=".", col=NULL, cex=1, lwd=2, changepoints=FALSE, grid=FALSE, quantiles=c(0.05,0.95), xlim=NULL, Clim=c(0,3*ploidy(x)), Blim=c(0,1), xScale=1e-6, ..., add=FALSE, subplots=!add && (length(panels) > 1), verbose=FALSE) {

  # To please R CMD check
  fit <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fit':
  if (nbrOfChromosomes(fit) > 1) {
    throw("Multiple chromosomes detected. Not yet implemented.");
  }

  # Argument 'panels':
  panels <- unique(panels);
  panelsOrg <- panels;
  panels <- gsub("[-*]", "", panelsOrg);
  knownPanels <- c("tcn", "tcn,c1,c2", "tcn,c1", "tcn,c2", "c1,c2", "c1", "c2", "dh", "betaN", "betaT", "betaTN");
  panels <- match.arg(panels, choices=knownPanels, several.ok=TRUE);
  panels <- panelsOrg;

  # Argument 'calls':
  if (!is.null(calls)) {
    calls <- sapply(calls, FUN=Arguments$getRegularExpression);
  }

  # Argument 'changepoints':
  changepoints <- Arguments$getLogical(changepoints);

  # Argument 'grid':
  grid <- Arguments$getLogical(grid);

  # Argument 'xlim':
  if (!is.null(xlim)) {
    xlim <- Arguments$getNumerics(xlim, length=c(2,2));
  }

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf));

  # Argument 'add':
  add <- Arguments$getLogical(add);

  # Argument 'Clim':
  if (!add) {
    Clim <- Arguments$getNumerics(Clim, length=c(2L,2L),
                                        disallow=c("Inf", "NA", "NaN"));
  }

  # Argument 'subplots':
  subplots <- Arguments$getLogical(subplots);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Plotting PSCN panels");

  # Extract the input data
  data <- getLocusData(fit);
  if (is.null(data)) {
    throw("Cannot plot segmentation results. No input data available.");
  }

  chromosomes <- getChromosomes(fit);
  chromosome <- chromosomes[1];
  x <- data$x;
  CT <- data$CT;
  rho <- data$rho;
  betaT <- data$betaT;
  betaN <- data$betaN;
  betaTN <- data$betaTN;
  muN <- data$muN;
  rho <- data$rho
  hasDH <- !is.null(rho)
  if (hasDH) {
    isHet <- !is.na(rho)
    isSnp <- isHet
  } else {
    isSnp <- (!is.na(betaTN) & !is.na(muN))
    isHet <- isSnp & (muN == 1/2)
  }
  nbrOfLoci <- length(x)

  # BACKWARD COMPATIBILITY:
  # If 'rho' is not available, recalculate it from tumor BAFs.
  # NOTE: This should throw an error in the future. /HB 2013-10-25
  if (!hasDH) {
    rho <- rep(NA_real_, length=nbrOfLoci);
    rho[isHet] <- 2*abs(betaTN[isHet]-1/2);
    warning(sprintf("Locus-level DH signals ('rho') were not available in the %s object and therefore recalculated from the TumorBoost-normalized tumor BAFs ('betaTN').", class(fit)[1L]));
  }

  # Extract the segmentation
  segs <- as.data.frame(fit);

  # Identify available calls
  if (!is.null(calls)) {
    verbose && enter(verbose, "Identifying calls");

    pattern <- "Call$";
    callColumns <- grep(pattern, colnames(segs), value=TRUE);
    if (length(callColumns) > 0) {
      keep <- sapply(calls, FUN=function(pattern) {
        (regexpr(pattern, callColumns) != -1);
      });
      if (is.matrix(keep)) {
        keep <- rowAnys(keep);
      }
      callColumns <- callColumns[keep];
      callLabels <- gsub(pattern, "", callColumns);
      callLabels <- toupper(callLabels);
    }
    verbose && cat(verbose, "Call columns:");
    verbose && print(verbose, callColumns);

    verbose && exit(verbose);
  } else {
    callColumns <- NULL;
  }

  if (chromosome != 0) {
    chrTag <- sprintf("Chr%02d", chromosome);
  } else {
    chrTag <- "";
  }

  if (xScale != 1) {
    x <- xScale * x;
    if (is.null(xlim)) {
      xlim <- range(x, na.rm=TRUE);
    } else {
      xlim <- xScale * xlim;
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Number of panels?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfPanels <- length(panels);
  if (subplots) {
    subplots(nbrOfPanels, ncol=1);
    par(mar=c(1,4,1,2)+1);
  }

  # Color loci by genotype
  if (hasDH) {
    colMu <- c("gray", "black")[!is.na(rho) + 1]
  } else {
    colMu <- c("gray", "black")[(muN == 1/2) + 1]
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For each panel...
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (pp in seq(length=nbrOfPanels)) {
    panel <- panels[pp];
    verbose && enter(verbose, sprintf("Panel #%d ('%s') of %d",
                                             pp, panel, length(panels)));


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Setup empty plot
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (!add) {
      verbose && enter(verbose, "Creating empty plot");

      tracks <- strsplit(panel, split=",", fixed=TRUE)[[1]];
      tracks <- gsub("[-*]", "", tracks);

      # Defaults
      ylim <- Clim;
      ylab <- paste(toupper(tracks), collapse=", ");

      if (any(is.element(c("betaN", "betaT", "betaTN", "dh"), tracks))) {
        ylim <- Blim;
      }

      verbose && cat(verbose, "ylim:");
      verbose && print(verbose, ylim);

      plot(NA, xlim=xlim, ylim=ylim, ylab=ylab);

      # Geometrical annotations
      stext(side=3, pos=1, chrTag);
      if (grid) {
        abline(h=seq(from=0, to=ylim[2], by=2), lty=3, col="gray");
        abline(h=0, lty=1, col="black");
      }

      verbose && exit(verbose);
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Scatter tracks
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tracks <- strsplit(panel, split=",", fixed=TRUE)[[1]];
    keep <- grep("-", tracks, fixed=TRUE, invert=TRUE);
    tracks <- tracks[keep];

    verbose && cat(verbose, "Number tracks with scatter: ", length(tracks));
    verbose && cat(verbose, "Tracks with scatter:");
    verbose && print(verbose, tracks);
    nbrOfTracks <- length(tracks);

    for (tt in seq(length=nbrOfTracks)) {
      track <- tracks[tt];
      verbose && enter(verbose, sprintf("Scatter track #%d ('%s') of %d",
                                                tt, track, nbrOfTracks));
      track <- gsub("[-*]", "", track);

      # Defaults
      colT <- if (is.null(col)) "black" else col;

      if (track == "tcn") {
        y <- CT;
      } else if (track == "c1") {
        y <- 1/2*(1-rho)*CT;
      } else if (track == "c2") {
        y <- 1/2*(1+rho)*CT;
      } else if (track == "betaN") {
        y <- betaN;
        colT <- if (is.null(col)) colMu else col;
        ylab <- expression(BAF[N]);
        ylim <- Blim;
      } else if (track == "betaT") {
        y <- betaT;
        colT <- if (is.null(col)) colMu else col;
        ylab <- expression(BAF[T]);
        ylim <- Blim;
      } else if (track == "betaTN") {
        y <- betaTN;
        colT <- if (is.null(col)) colMu else col;
        ylab <- expression(BAF[T]^"*");
        ylim <- Blim;
      } else if (track == "dh") {
        y <- rho;
        colT <- if (is.null(col)) colMu[isHet] else col;
        ylab <- "DH";
        ylim <- Blim;
      } else {
        y <- NULL;
      }

      # Nothing to do?
      if (!is.null(y)) {
        points(x, y, pch=pch, cex=cex, col=colT);
      }

      verbose && exit(verbose);
    } # for (tt ...)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Segment tracks
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tracks <- strsplit(panel, split=",", fixed=TRUE)[[1]];
    keep <- grep("*", tracks, fixed=TRUE);
    tracks <- gsub("[-*]", "", tracks[keep]);

    # Keep only supported tracks
    tracksWithLevels <- c("tcn", "betaTN", "c1", "c2", "dh");
    stopifnot(all(is.element(tracks, tracksWithLevels)));
    tracks <- intersect(tracks, tracksWithLevels);

    verbose && cat(verbose, "Number tracks with levels: ", length(tracks));
    verbose && cat(verbose, "Tracks with levels:");
    verbose && print(verbose, tracks);
    nbrOfTracks <- length(tracks);

    for (tt in seq(length=nbrOfTracks)) {
      track <- tracks[tt];
      verbose && enter(verbose, sprintf("Level track #%d ('%s') of %d",
                                                 tt, track, nbrOfTracks));

      if (track == "tcn") {
        colT <- "purple";
      } else if (track == "c1") {
        colT <- "blue";
      } else if (track == "c2") {
        colT <- "red";
      } else if (track == "betaTN") {
        colT <- "orange";
      } else if (track == "dh") {
        colT <- "orange";
      } else {
        colT <- if (is.null(col)) "black" else col;
      }

      # Nothing to do?
      if (track != "betaTN") {
        drawConfidenceBands(fit, what=track, quantiles=quantiles,

                        col=colT, xScale=xScale);
      }
      drawLevels(fit, what=track, col=colT, lwd=lwd, xScale=xScale);

      verbose && exit(verbose);
    } # for (tt ...)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Draw change points?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (changepoints) {
      segs <- as.data.frame(fit);
      xStarts <- segs[,"tcnStart"];
      xEnds <- segs[,"tcnEnd"];
      xs <- sort(unique(c(xStarts, xEnds)));
      abline(v=xScale*xs, lty=1, col="gray");
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # For each panel of tracks, annotate with calls?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (length(callColumns) > 0) {
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
    } # if (length(callColumns) > 0)

    verbose && exit(verbose);
  } # for (pp ...)

  verbose && exit(verbose);

  invisible();
}, protected=TRUE) # plotTracks2()



############################################################################
# HISTORY:
# 2011-09-30
# o BUG FIX: plotTracks2(..., panels="dh") gave an error due to a
#   forgotten assigment.
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-01-19
# o Added plotTracks2().  Completely rewritten plotTracks().
# o Created.
############################################################################
