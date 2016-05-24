###########################################################################/**
# @set "class=CBS"
# @RdocMethod plotTracks
#
# @title "Plots copy numbers along the genome"
#
# \description{
#  @get "title" for one or more chromosomes.
#  Each type of track is plotted in its own panel.
# }
#
# @synopsis
#
# \arguments{
#   \item{x}{A result object returned by @see "segmentByCBS".}
#   \item{pch}{The type of points to use.}
#   \item{Clim}{The range of copy numbers.}
#   \item{xScale}{The scale factor used for genomic positions.}
#   \item{...}{Not used.}
#   \item{add}{If @TRUE, the panels plotted are added to the existing plot,
#     otherwise a new plot is created.}
# }
#
# \value{
#   Returns nothing.
# }
#
# @author "HB"
#
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("plotTracks", "CBS", function(x, scatter=TRUE, pch=20, col="gray", meanCol="purple", cex=1, grid=FALSE, Clim="auto", xScale=1e-6, Clab="auto", ..., byIndex=FALSE, mar=NULL, add=FALSE) {
  # To please R CMD check
  fit <- x;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'add':
  add <- Arguments$getLogical(add);

  # Argument 'Clim':
  if (identical(Clim, "auto")) {
    signalType <- getSignalType(fit);
    ploidy <- ploidy(fit);
    Clim <- switch(signalType,
      "log2ratio" = c(-2,2) + c(-1,1)*ploidy/2,
      "ratio"     = c(0,3*ploidy),
      NULL
    );
##  NOTE: Don't understand why, but with this 'R CMD build' gives:
##    "Error: processing vignette 'CBS.tex.rsp' failed with diagnostics:
##     Failed to infer argument 'Clim' due to an unknown signalType(): NA"
##  /HB 2013-10-14
##  if (!add && is.null(Clim)) {
##    throw("Failed to infer argument 'Clim' due to an unknown signalType(): ", signalType);
##  }
  } else if (!add) {
    Clim <- Arguments$getNumerics(Clim, length=c(2L,2L),
                                        disallow=c("Inf", "NA", "NaN"));
  }

  if (identical(Clab, "auto")) {
    signalType <- getSignalType(fit);
    Clab <- switch(signalType,
      "log2ratio" = "log2 CN ratio",
      "ratio"     = "CN ratio",
      NULL
    );
  }

  # Argument 'fit':
  if (nbrOfChromosomes(fit) > 1L) {
    res <- plotTracksManyChromosomes(fit, scatter=scatter, pch=pch, col=col, cex=cex, meanCol=meanCol, Clim=Clim, xScale=xScale, Clab=Clab, ..., byIndex=byIndex, mar=mar, add=add);
    return(invisible(res));
  }

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf));

  # Extract the input data
  data <- getLocusData(fit);
  if (is.null(data)) {
    throw("Cannot plot segmentation results. No input data available.");
  }

  chromosomes <- getChromosomes(fit);
  chromosome <- chromosomes[1];
  x <- data$x;
  CT <- data[,3];
  nbrOfLoci <- length(x);

  # Extract the segmentation
  segs <- getSegments(fit);


  if (chromosome != 0) {
    chrTag <- sprintf("Chr%02d", chromosome);
  } else {
    chrTag <- "";
  }

  if (xScale != 1) {
    x <- xScale * x;
  }

  if (!add && !is.null(mar)) {
    par(mar=mar);
  }


  pchT <- if (scatter) { pch } else { NA };

  plot(x, CT, pch=pchT, cex=cex, col=col, ..., ylim=Clim, ylab=Clab);
  stext(side=3, pos=1, chrTag);
  if (grid) {
    yrange <- par("usr")[3:4];
    yrange[1] <- floor(yrange[1]);
    yrange[2] <- ceiling(yrange[2]);
    abline(h=seq(from=yrange[1], to=yrange[2], by=2), lty=3, col="gray");
    abline(h=0, lty=1, col="black");
  }
  drawLevels(fit, col=meanCol, xScale=xScale);

  invisible();
}) # plotTracks()


setMethodS3("plot", "CBS", function(x, ...) {
  plotTracks(x, ...);
}, protected=TRUE)


setMethodS3("drawLevels", "CBS", function(fit, col="purple", xScale=1e-6, byIndex=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- tileChromosomes(fit);

  # Get segmentation results
  segs <- as.data.frame(fitT);

  # Extract subset of segments
  fields <- c("start", "end", "mean");
  segs <- segs[,fields, drop=FALSE];
  segs <- unique(segs);

  # Reuse drawLevels() for the DNAcopy class
  colnames(segs) <- c("loc.start", "loc.end", "seg.mean");
  dummy <- list(output=segs);
  class(dummy) <- "DNAcopy";
  drawLevels(dummy, col=col, xScale=xScale, ...);
}, protected=TRUE)


setMethodS3("highlightCalls", "CBS", function(fit, pch=20, callCols=c(loss="red", gain="green", "amplification"="blue"), lwd=3, meanCol="purple", ..., xScale=1e-6, byIndex=FALSE, verbose=FALSE) {
  segs <- getSegments(fit, splitter=FALSE);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify segment calls
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  callFields <- grep("Call$", colnames(segs));
  callTypes <- gsub("Call$", "", colnames(segs)[callFields]);
  nbrOfCalls <- length(callFields);

  # Nothing todo?
  if (nbrOfCalls == 0L) {
    return();
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- tileChromosomes(fit);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Highlight threshold levels
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- fit$params$callGainsAndLosses;
  abline(h=params$muR, col="gray", lty=3);
  abline(h=params$tauLoss, col=callCols["loss"], lty=3);
  abline(h=params$tauGain, col=callCols["gain"], lty=3);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Highlight gains and losses
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dataT <- getLocusData(fitT);
  segsT <- getSegments(fitT, splitter=FALSE);
  chr <- dataT[,"chromosome"];
  x <- dataT[,"x"];
  y <- dataT[,3];
  nbrOfLoci <- nbrOfLoci(fitT);
  nbrOfSegments <- nbrOfSegments(fitT);
  # Not needed anymore
  dataT <- NULL;

  # For each non-neutral segment
  for (ss in seq(length=nbrOfSegments)) {
    seg <- segsT[ss,];

    for (tt in seq(along=callTypes)) {
      field <- callFields[tt];
      type <- callTypes[tt];

      # Called?
      call <- seg[[field]];
      if (isTRUE(call)) {
        col <- callCols[type];
        idxs <- which(chr == seg$chromosome & seg$start <= x & x <= seg$end);
        idxs <- Arguments$getIndices(idxs, max=nbrOfLoci);

        if (byIndex) {
          xs <- idxs;
        } else {
          xs <- x[idxs] * xScale;
        }
        ys <- y[idxs];
        points(xs, ys, pch=pch, col=col, ...);
        xx <- range(xs, na.rm=TRUE);
        yy <- rep(seg$mean, times=2);
        lines(xx, yy, lwd=lwd, col=meanCol);
      }
    } # for (tt ...)
  } # for (ss ...)
}, protected=TRUE) # highlightCalls()



setMethodS3("highlightLocusCalls", "CBS", function(fit, callPchs=c(negOutlier=25, posOutlier=24), callCols=c(negOutlier="blue", posOutlier="blue"), ..., xScale=1e-6, byIndex=FALSE, verbose=FALSE) {
  data <- getLocusData(fit);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify segment calls
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  callFields <- grep("Call$", colnames(data));
  callTypes <- gsub("Call$", "", colnames(data)[callFields]);
  nbrOfCalls <- length(callFields);

  # Nothing todo?
  if (nbrOfCalls == 0) {
    return();
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- tileChromosomes(fit);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Highlight gains and losses
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dataT <- getLocusData(fitT);
  chr <- dataT[,"chromosome"];
  x <- dataT[,"x"];
  y <- dataT[,3];
  nbrOfLoci <- nbrOfLoci(fitT);

  # For each non-neutral segment
  for (tt in seq(along=callTypes)) {
    field <- callFields[tt];
    type <- callTypes[tt];

    isCalled <- dataT[[field]];
    idxs <- which(isCalled);

    if (length(idxs) == 0L) {
      next;
    }

    if (byIndex) {
      xs <- idxs;
    } else {
      xs <- x[idxs] * xScale;
    }
    ys <- y[idxs];
    pch <- callPchs[type];
    col <- callCols[type];
    points(xs, ys, pch=pch, col=col, ...);
  } # for (tt ...)
}, protected=TRUE) # highlightLocusCalls()




setMethodS3("drawChromosomes", "CBS", function(x, lty=3, xScale=1e-6, ..., byIndex=FALSE, verbose=FALSE) {
  # To please R CMD check
  fit <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fit':

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- tileChromosomes(fit);
  # Sanity check
  stopifnot(!is.null(fitT$chromosomeStats));

  chrStats <- fitT$chromosomeStats;
  chrStats <- chrStats[-nrow(chrStats),,drop=FALSE];
  chrRanges <- as.matrix(chrStats[,c("start","end")]);
  vs <- xScale * chrRanges;
  mids <- (vs[,1]+vs[,2])/2;
  chromosomes <- getChromosomes(fitT);
  chrLabels <- sprintf("%02d", chromosomes);
  side <- rep(c(1,3), length.out=length(chrLabels));
  mtext(text=chrLabels, side=side, at=mids, line=0.1, cex=0.7*par("cex"));
  abline(v=vs, lty=lty);
}, protected=TRUE) # drawChromosomes()



setMethodS3("drawCentromeres", "CBS", function(fit, genomeData, what=c("start", "end"), xScale=1e-6, col="gray", lty=3, ..., byIndex=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'genomeData':
  stopifnot(inherits(genomeData, "data.frame"));
  stopifnot(is.element("chromosome", colnames(genomeData)));
  stopifnot(is.element("centroStart", colnames(genomeData)));
  stopifnot(is.element("centroEnd", colnames(genomeData)));

  # Calculate the midpoints of the centromeres
  colnames(genomeData) <- tolower(gsub("centro", "", colnames(genomeData)));
  genomeData$mid <- (genomeData$start + genomeData$end) / 2;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- tileChromosomes(fit);
  # Sanity check
  stopifnot(!is.null(fitT$chromosomeStats));

  chrStats <- fitT$chromosomeStats;
  offsets <- chrStats[,"start"] - chrStats[1,"start"];

  # Centroid locations in the tiled space
  offsetsT <- offsets[seq(length=nrow(genomeData))];

  xx <- genomeData[,what,drop=FALSE];
  xx <- as.matrix(xx);
  xx <- offsetsT + xx;

  ats <- xScale * xx;
  for (cc in seq(length=ncol(xx))) {
    abline(v=ats[,cc], col=col, lty=lty, ...);
  }

  invisible(ats);
}, protected=TRUE) # drawCentromeres()


setMethodS3("highlightArmCalls", "CBS", function(fit, genomeData, minFraction=0.95, callCols=c("loss"="red", "gain"="green"),   xScale=1e-6, ...) {
  # To please/trick R CMD check
  chromosome <- x <- NULL; rm(list=c("chromosome", "x"));

  callStats <- callArms(fit, genomeData=genomeData, minFraction=minFraction);

  callTypes <- grep("Fraction", colnames(callStats), value=TRUE);
  callTypes <- gsub("Fraction", "", callTypes);

  callTypes <- intersect(callTypes, c("loss", "gain"));

  # Adjust (start, end)
  offsets <- getChromosomeOffsets(fit);
  offsets <- offsets[callStats[,"chromosome"]];
  callStats[,c("start","end")] <- offsets + callStats[,c("start","end")];

  nbrOfRegions <- nrow(callStats);

  # Nothing todo?
  if (nbrOfRegions == 0) {
    return(invisible(callStats));
  }


  usr <- par("usr");
  dy <- diff(usr[3:4]);
  yy <- usr[3]+c(0,0.05*dy);
  abline(h=usr[3]+0.95*0.05*dy, lty=1, col="gray");

  xx <- callStats[,c("start", "end")];
  xx <- as.matrix(xx);
  xx <- xx * xScale;

  for (type in callTypes) {
    col <- callCols[type];
    keyA <- sprintf("%sFraction", type);
    keyB <- sprintf("%sCall", type);
    for (kk in seq(length=nbrOfRegions)) {
      xs <- xx[kk,];
      score <- callStats[kk, keyA];
      if (is.finite(score) && score > 0) {
        ys <- rep(yy[1]+callStats[kk, keyA]*0.05*dy, times=2);
        lines(x=xs, y=ys, col=col);
        call <- callStats[kk, keyB];
        if (call) {
          rect(xs[1], yy[1], xs[2], yy[2], col=col, border=NA);
        }
      }
    }
  } # for (type ...)

  invisible(callStats);
}, protected=TRUE); # highlightArmCalls()



############################################################################
# HISTORY:
# 2013-11-23
# o Added "dummy" 'byIndex' argument to drawLevels() for CBS to avoid
#   that argument being passed down lines().
# 2013-10-14
# o Now plotTracks() for CBS gives a more informative error if 'Clim'
#   is invalid or "auto" and could not be inferred due to an unknown
#   or unset signal type.
# 2013-04-18
# o Now drawLevels() also works for multiple chromosomes.
# 2011-12-06
# o Now plotTracks() for CBS always returns an invisible object.
# 2011-12-03
# o Added drawChangePoints() for CBS.
# 2011-10-23
# o BUG FIX: highlightArmCalls() for CBS did not handle empty chromosomes.
# 2011-10-08
# o Added drawChromosomes() for CBS.
# 2011-10-07
# o Added highlightArmCalls() for CBS.
# 2011-10-06
# o Now drawCentromeres() for CBS can also plot start and stop.
# o ROBUSTNESS: Now plotTracksManyChromosomes() extract (start,end)
#   information by names (no longer assuming certain indices).
# 2011-09-07
# o Added highlightLocusCalls().
# 2011-09-06
# o Added highlightCalls().
# o Added getChromosomeOffsets().
# 2011-09-01
# o BUG FIX: plotTracksManyChromosomes() for CBS gave an error because
#   internal variable 'CT' was not defined.
# o BUG FIX: tileChromosomes() for CBS not identify the chromosomes of
#   the loci, and hence generated corrupt/missing values while tiling.
# 2010-11-19
# o Created from PairedPSCBS.R.
############################################################################
