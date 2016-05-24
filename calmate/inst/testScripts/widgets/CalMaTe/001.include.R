###########################################################################
# 001.include.R
#
# Author: Henrik Bengtsson
# 
# Description:
# This script defines functions and other objects that are needed by
# all scripts in this directory.
###########################################################################
library("aroma.cn");
library("calmate");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
loadSets <- function(..., verbose=FALSE) {
  dsT <- AromaUnitTotalCnBinarySet$byName(..., verbose=verbose);
  dsB <- AromaUnitFracBCnBinarySet$byName(..., verbose=verbose);
  list(total=dsT, fracB=dsB);
} # loadSets()



extractSignals <- function(dsList, sampleName, ..., reference=c("auto", "none", "median"), verbose=FALSE) {
  # Argument 'reference':
  reference <- match.arg(reference);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extracting total CN signals");
  
  verbose && enter(verbose, "Extracting sample of interest");
  verbose && cat(verbose, "Sample name: ", sampleName);
  
  idxT <- indexOf(dsList$total, sampleName);
  dfT <- getFile(dsList$total, idxT);
  idxB <- indexOf(dsList$fracB, sampleName);
  dfB <- getFile(dsList$fracB, idxB);
  
  verbose && print(verbose, list(tumor=dfT, fracB=dfB));
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Extracting TCNs");
  tcn <- extractRawCopyNumbers(dfT, logBase=NULL, ...);
  verbose && print(verbose, tcn);
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting BAFs");
  baf <- extractRawAlleleBFractions(dfB, ...);
  verbose && print(verbose, baf);
  verbose && exit(verbose);

  if (reference == "auto") {
    verbose && enter(verbose, "Inferring from data set name if CN ratios needs to be calculated");
    hasRatios <- hasTag(dsList$total, "CMTN");
    verbose && cat(verbose, "Has 'CMTN' tag: ", hasRatios);
    reference <- ifelse(hasRatios, "none", "median");
    verbose && cat(verbose, "Inferred argument 'reference': ", reference);
    verbose && exit(verbose);
  }
  
  if (reference == "median") {
    dfTR <- getAverageFile(dsList$total, verbose=less(verbose,5));
    
    verbose && enter(verbose, "Extracting TCNs for reference pool");
    tcnR <- extractRawCopyNumbers(dfTR, logBase=NULL, ...);
    verbose && print(verbose, tcnR);
    verbose && exit(verbose);
    
    tcn <- divideBy(tcn, tcnR);
    tcn$y <- 2*tcn$y;
  }
  
  res <- list(tcn=tcn, baf=baf);
  verbose && exit(verbose);
  
  res;
} # extractSignals()


plotBvsB <- function(beta, col=NULL, pch=19, Blim=c(-0.2,+1.2), xlab="Normal BAF", ylab="Tumor BAF", dataSet=NULL, chipType=NULL, sampleName=NULL, chrTag=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
  par(mar=c(5.8,7,3,3)+0.1, mgp=c(4,1.4,0), cex=0.8, cex.lab=3.8, cex.axis=2.8);
  plot(NA, xlim=Blim, ylim=Blim, xlab=xlab, ylab=ylab, axes=FALSE);
  for (ss in 1:2) {
    axis(side=ss, at=c(0,1/2,1), label=c("0","1/2","1"));
  }
  box();
  abline(a=0, b=1, lty=3, lwd=2);

  # Text annotations
  if (!is.null(chrTag)) {
    stext(side=3, pos=0, chrTag,  cex=1.5);
  } else if (!is.null(segName) & !is.null(segTag)) {
    stext(side=3, pos=0, sprintf("%s @ %s", segName, segTag), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName,  cex=1.5);
  }

  points(beta, pch=pch, col=col, cex=1);

  if (is.element("dens", anTags)) {
    for (dd in 1:2) {
      d <- density(beta[,dd], adjust=0.4, from=Blim[1], to=Blim[2], na.rm=TRUE);
      draw(d, lwd=3, col="#666666", side=dd, 
                       height=0.05, scale="relative", xpd=FALSE);
    }
  }
} # plotBvsB()



plotASCN <- function(ascn, col=NULL, pch=19, Clim=c(-0.4,4), xlab=expression(C[A]), ylab=expression(C[B]), dataSet=NULL, chipType=NULL, sampleName=NULL, chrTag=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
#  par(mar=c(4,5,3,3)+0.1, cex=0.8, cex.lab=2.4, cex.axis=2.2);
  par(mar=c(5.8,7,3,3)+0.1, mgp=c(4,1.4,0), cex=0.8, cex.lab=3.8, cex.axis=2.8);
  plot(NA, xlim=Clim, ylim=Clim, xlab=xlab, ylab=ylab, axes=FALSE);
  for (ss in 1:2) {
    axis(side=ss, at=0:round(max(Clim)+1L));
  }
  box();
  abline(a=0, b=1, lty=3, lwd=2);

  # Text annotations
  if (!is.null(chrTag)) {
    stext(side=3, pos=0, chrTag,  cex=1.5);
  } else if (!is.null(segName) & !is.null(segTag)) {
    stext(side=3, pos=0, sprintf("%s @ %s", segName, segTag), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName,  cex=1.5);
  }

  points(ascn, pch=pch, col=col, cex=1);

  if (is.element("dens", anTags)) {
    for (dd in 1:2) {
      d <- density(ascn[,dd], adjust=0.4, from=Clim[1], to=Clim[2], na.rm=TRUE);
      draw(d, lwd=3, col="#666666", side=dd, 
                       height=0.05, scale="relative", xpd=FALSE);
    }
  }
} # plotASCN()



plotTrackC1C2 <- function(ascn, x, muN, pch=".", lwd=3, xlim=NULL, Clim=c(-0.4,3.4), xlab="Position (Mb)", ylab=expression(list(C[1],C[2])), dataSet=NULL, chipType=NULL, sampleName=NULL, chrTag=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
  if (is.null(xlim)) {
    xlim <- range(x, na.rm=TRUE);
    dx <- diff(xlim);
    xlim[1] <- xlim[1] - 0.05*dx;
    xlim[2] <- xlim[2] + 0.05*dx;
  }

#  par(mar=c(4,5,3,3)+0.1, cex=0.8, cex.lab=2.4, cex.axis=2.2);
#  par(mar=c(5.8,7,3,3)+0.1, mgp=c(4,1.4,0), cex=0.8, cex.lab=3.8, cex.axis=2.8);
  par(mar=c(2.5,2.5,1.1,1)+0.1, tcl=-0.3, mgp=c(1.4,0.4,0), cex=2);
  plot(NA, xlim=xlim, ylim=Clim, xlab=xlab, ylab=ylab, axes=FALSE);
  axis(side=1);
  axis(side=2, at=seq(from=0, to=as.integer(Clim[2]+1L)));
  box();

  # Text annotations
  if (!is.null(chrTag)) {
    stext(side=3, pos=0, chrTag,  cex=1.5);
  } else if (!is.null(segName) & !is.null(segTag)) {
    stext(side=3, pos=0, sprintf("%s @ %s", segName, segTag), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
#    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName,  cex=1.5);
  }

  hets <- which(muN == 1/2);
  xH <- x[hets];
  c1c2 <- ascn[hets,];
  rr <- which(c1c2[,2] < c1c2[,1]);
  c1c2[rr,] <- c1c2[rr,2:1];

  c1 <- c1c2[,1];
  c2 <- c1c2[,2];

  points(xH, c1, pch=pch, cex=1, col="blue");
  points(xH, c2, pch=pch, cex=1, col="red");

  if (!is.null(segName)) {
    abline(h=median(c1,na.rm=TRUE), lwd=lwd, col="blue");
    abline(h=median(c2,na.rm=TRUE), lwd=lwd, col="red");

    if (is.element("dens", anTags)) {
      d <- density(c1, adjust=0.8, na.rm=TRUE);
      draw(d, lwd=3, col="blue", side=2, height=0.05, scale="relative", xpd=FALSE);
  
      d <- density(c2, adjust=0.8, na.rm=TRUE);
      draw(d, lwd=3, col="red", side=2, height=0.05, scale="relative", xpd=FALSE);
    }
  }

  by <- NULL;
  if (is.element("1Mb", anTags)) {
    by <- 1;
  } else if (is.element("2Mb", anTags)) {
    by <- 2;
  } else if (is.element("5Mb", anTags)) {
    by <- 5;
  }

  if (!is.null(by)) {
    c1s <- binnedSmoothing(c1, xH, by=by);
    c2s <- binnedSmoothing(c2, xH, by=by);
    xOut <- attr(c1s, "xOut");
    lwd <- 4;
    lines(xOut, c1s, lwd=lwd+2, col="white");
    lines(xOut, c2s, lwd=lwd+2, col="white");
    lines(xOut, c1s, lwd=lwd, col="blue");
    lines(xOut, c2s, lwd=lwd, col="red");
  }
} # plotTrackC1C2()


plotTrackTCN <- function(C, x, col=NULL, pch=".", xlim=NULL, Clim=c(-0.4,5.4), ylab="Total CN", xlab="Position (Mb)", dataSet=NULL, chipType=NULL, sampleName=NULL, chrTag=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
  par(mar=c(2.5,2.5,1.1,1)+0.1, tcl=-0.3, mgp=c(1.4,0.4,0), cex=2);
  plot(x, C, pch=pch, col=col, xlim=xlim, ylim=Clim, xlab=xlab, ylab=ylab, axes=FALSE);
  axis(side=1);
  axis(side=2, at=seq(from=0, to=as.integer(Clim[2]+1L)));
  box();

  # Text annotations
  if (!is.null(chrTag)) {
    stext(side=3, pos=0, chrTag,  cex=1.5);
  } else if (!is.null(segName) & !is.null(segTag)) {
    stext(side=3, pos=0, sprintf("%s @ %s", segName, segTag), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
#    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName, cex=1.5, line=-0.2);
  }
} # plotTrackTCN()



plotTrackBAF <- function(B, x, col=NULL, pch=".", xlim=NULL, Blim=c(-0.2,+1.2), ylab="BAF", xlab="Position (Mb)", dataSet=NULL, chipType=NULL, sampleName=NULL, chrTag=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
  par(mar=c(2.5,2.5,1.1,1)+0.1, tcl=-0.3, mgp=c(1.4,0.4,0), cex=2);
  plot(x, B, pch=pch, col=col, xlim=xlim, ylim=Blim, xlab=xlab, ylab=ylab, axes=FALSE);
  axis(side=1);
  axis(side=2, at=c(0, 1/2, 1), label=c(0, "1/2", 1));
  box();

  # Text annotations
  if (!is.null(chrTag)) {
    stext(side=3, pos=0, chrTag,  cex=1.5);
  } else if (!is.null(segName) & !is.null(segTag)) {
    stext(side=3, pos=0, sprintf("%s @ %s", segName, segTag), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
#    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName, cex=1.5, line=-0.2);
  }
} # plotTrackBAF()



plotTrackDH <- function(B, x, muN, col=NULL, pch=".", xlim=NULL, Blim=c(0,+1.2), ylab="DH", xlab="Position (Mb)", dataSet=NULL, chipType=NULL, sampleName=NULL, chrTag=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
  par(mar=c(2.5,2.5,1.1,1)+0.1, tcl=-0.3, mgp=c(1.4,0.4,0), cex=2);
  DH <- 2*abs(B - 1/2);
  DH[muN != 1/2] <- NA;
  plot(x, DH, pch=pch, col=col, xlim=xlim, ylim=Blim, xlab=xlab, ylab=ylab, axes=FALSE);
  axis(side=1);
  axis(side=2, at=c(0, 1), label=c(0, 1));
  box();

  # Text annotations
  if (!is.null(chrTag)) {
    stext(side=3, pos=0, chrTag,  cex=1.5);
  } else if (!is.null(segName) & !is.null(segTag)) {
    stext(side=3, pos=0, sprintf("%s @ %s", segName, segTag), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
#    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName, cex=1.5, line=-0.2);
  }
} # plotTrackDH()




extractCACB <- function(dsList, units, ..., reference=c("auto", "none", "median"), samplesLast=FALSE) {
  # Argument 'reference':
  reference <- match.arg(reference);

  C <- extractMatrix(dsList$total, units=units, ...);
  B <- extractMatrix(dsList$fracB, units=units, ...);
  stopifnot(identical(dim(B), dim(C)));

  if (reference == "auto") {
    verbose && enter(verbose, "Inferring from data set name if CN ratios needs to be calculated");
    hasRatios <- hasTag(dsList$total, "CMTN");
    verbose && cat(verbose, "Has 'CMTN' tag: ", hasRatios);
    reference <- ifelse(hasRatios, "none", "median");
    verbose && cat(verbose, "Inferred argument 'reference': ", reference);
    verbose && exit(verbose);
  }

  if (reference == "median") {
    muC <- matrixStats::rowMedians(C, na.rm=TRUE);
    C <- 2 * C / muC;
  }

  CB <- B*C;
  CA <- C - CB;
  data <- array(c(CA, CB), dim=c(dim(C), 2));
  dimnames(data)[[3]] <- c("A", "B");

  if (samplesLast) {
    data <- aperm(data, perm=c(1,3,2));
  }

  data;
} # extractCACB()


plotMultiArrayCACB <- function(cacb, pch=19, cex=2, ..., Clim=c(-0.2,3.4), xlab=expression(C[A]), ylab=expression(C[B]), dataSet=NULL, chipType=NULL, unitName=NULL, tagsT=NULL) {
  par(mar=c(5.8,7,3,3)+0.1, mgp=c(4,1.4,0), cex=0.8, cex.lab=3.8, cex.axis=2.8);
  plot(NA, xlim=Clim, ylim=Clim, xlab=xlab, ylab=ylab, axes=FALSE);
  for (ss in 1:2) {
    axis(side=ss, at=0:round(max(Clim)+1L));
  }
  box();
  lines(x=c(0,2), y=c(2,0), lty=2, lwd=2, col="#999999");
  lines(x=c(0,Clim[2]), y=c(0,Clim[2]), lty=2, lwd=2, col="#999999");
  lines(x=c(0,0), y=c(0,Clim[2]), lty=2, lwd=2, col="#999999");
  lines(x=c(0,Clim[2]), y=c(0,0), lty=2, lwd=2, col="#999999");

  if (!is.null(unitName)) {
    stext(side=3, pos=0, sprintf("Unit: %s", unitName), cex=1.5);
  }

  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }

  if (!is.null(dataSet) & !is.null(chipType)) {
    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }

  stext(side=4, pos=1, sprintf("n=%d", nrow(cacb)), cex=1.5);

  points(cacb, pch=pch, cex=cex, ...);
} # plotMultiArrayCACB()

##########################################################################
# HISTORY:
# 2011-03-09 [HB]
# o Added option "auto" for argument 'reference' of extractSignals().
# o Created.
###########################################################################
