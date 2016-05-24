library("aroma.cn");
library("PSCBS");
library("R.devices");
library("R.menu");
verbose <- Arguments$getVerbose(-10);

# Local functions
deShearC1C2 <- deShearC1C2_20120922;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("doPlots", "PairedPSCBS", function(fit, sampleName=NULL, tags=NULL, ...) {
  # Argument 'sampleName':
  if (is.null(sampleName)) {
    sampleName <- sampleName(fit);
  }
  stopifnot(!is.null(sampleName));

  nCPsTag <- sprintf("#CPs=%d", nbrOfChangePoints(fit));
  toPNG(sampleName, tags=c("(C1,C2)", nCPsTag, tags), width=800, {
    plotC1C2Grid(fit);
    linesC1C2(fit);
    stext(side=3, pos=0, sampleName);
    stext(side=3, pos=1, nCPsTag);
    stext(side=4, pos=0, dataSet, cex=0.7);
    stext(side=4, pos=1, chipType, cex=0.7);
  });


  toPNG(sampleName, tags=c("tracks", nCPsTag, tags), width=1200, aspectRatio=0.25, {
    plotTracks(fit, tracks="tcn,c1,c2");
    stext(side=4, pos=0, sampleName);
    stext(side=4, pos=1, nCPsTag);
  });
}) # doPlots()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup Paired PSCBS segmentation data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "pscbsData";
path <- Arguments$getReadablePath(rootPath);

dataSets <- list.files(rootPath);
if (length(dataSets) > 1) {
 dataSet <- textMenu(dataSets, value=TRUE);
} else {
 dataSet <- dataSets[1];
}

path <- file.path(rootPath, dataSet);
path <- Arguments$getReadablePath(path);
chipTypes <- list.files(path);
if (length(chipTypes) > 1) {
 chipType <- textMenu(chipTypes, value=TRUE);
} else {
 chipType <- chipTypes[1];
}

ds <- PairedPSCBSFileSet$byName(dataSet, chipType=chipType);
print(ds);
dsName <- getName(ds);

if (length(ds) == 0) {
 throw("No PairedPSCBS data file found.")
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select tumor-normal pair
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (length(ds) > 1) {
 ii <- textMenu(getNames(ds));
} else {
 ii <- 1L;
}

if (!exists("fit") || !inherits(fit, "PairedPSCBS")) {
  df <- getFile(ds, ii);
  fit <- loadObject(df);
  sampleName <- getName(df);
  rm(segList, fitList);
}

fit0 <- fit;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Configure report
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
figPath <- file.path("figures", dataSet);
options("devEval/args/path"=figPath);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot (C1,C2)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
doPlots(fit);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Prune change points using dynamic programming
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("segList", mode="list")) {
  segList <- seqOfSegmentsByDP(fit, verbose=-10);
  modelFit <- attr(segList, "modelFit");
  modelFit$seqOfSegmentsByDP <- NULL;
  str(modelFit);
}


toPNG(sampleName, tags=c("DP", "RSEvsCPs"), width=800, aspectRatio=0.7, {
  plot(modelFit$nbrOfChangePoints, modelFit$rse,
       xlab="Number of change points", ylab="RSE");
  stext(side=3, pos=0, sampleName);
  stext(side=4, pos=0, dataSet, cex=0.7);
  stext(side=4, pos=1, chipType, cex=0.7);
});


nbrOfCPs <- c(100, 50, 25)[1:2];
if (!exists("fitList", mode="list")) {
  fitList <- list();
}
for (kk in seq(along=nbrOfCPs)) {
  key <- sprintf("nbrOfCPs=%d", nbrOfCPs[kk]);

  verbose && enter(verbose, sprintf("Change point set #%d ('%s') of %d", kk, key, length(nbrOfCPs)));

  verbose && cat(verbose, "Number of change points: ", nbrOfCPs[kk]);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Pruning CPs via dynamic programming
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitT <- fitList[[key]];
  if (is.null(fitT)) {
    verbose && enter(verbose, "Resegmenting");
    knownSegments <- segList[[nbrOfCPs[kk]+1L]];
    fitT <- resegment(fit, knownSegments=knownSegments, undoTCN=+Inf, undoDH=+Inf);
    fitList[[key]] <- fitT;
    verbose && exit(verbose);
  }
  sampleName(fitT) <- sampleName(fit);
  fitDP <- fitT;
  doPlots(fitDP);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Deshear
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitD <- deShearC1C2(fitDP);
  doPlots(fitD, tags="deShear");

  nCPsTag <- sprintf("#CPs=%d", nbrOfChangePoints(fitD));
  toPNG(sampleName, tags=c("cpCallDensity", nCPsTag, "deShear"), width=800, aspectRatio=0.5, {
    debug <- fitD$modelFit$debug;
    d <- debug$cpAngleDensity;
    pfp <- debug$pfp;
    expected <- attr(pfp, "expected");
    par(mar=c(5,4,2,2));
    plot(d, lwd=2, main="");
    abline(v=expected);
    text(x=expected, y=par("usr")[4], names(expected), adj=c(0.5,-0.5), cex=1.5, xpd=TRUE);

    # Annotate called peaks
    idxs <- match(pfp$call, expected);
    text(x=pfp$x, y=pfp$density, names(expected)[idxs], adj=c(0.5,-0.5), cex=1.5, col="blue");

    stext(side=4, pos=0, sampleName);
    stext(side=4, pos=1, nCPsTag);
  });

  verbose && exit(verbose);
} # for (kk ...)
