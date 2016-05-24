library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE13372,testset";
chipType <- "GenomeWideSNP_6,Full";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, verbose=verbose);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Spatial intensity plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ae <- ArrayExplorer(csR);
setColorMaps(ae, "sqrt,yellow");
print(ae);

process(ae, verbose=verbose);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Spatial probe log-ratio plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cfR <- getAverageFile(csR, verbose=verbose);
reporter <- SpatialReporter(csR, reference=cfR);
addColorMap(reporter, "log2center,rainbow");
process(reporter, zrange=c(-2,2), verbose=verbose);

ylab <- expression(log[2](y/y[R]));
for (ii in seq_along(csR)) {
  df <- csR[[ii]];

  toPNG(getFullName(df), tags="spatial,rowMedians", width=300, aspectRatio=8/3, {
    plotMargins(reporter, array=ii, ylim=c(-1,1)*0.2, ylab=ylab, margins="rows", rotate=90);
  });

  toPNG(getFullName(df), tags="spatial,colMedians", width=800, aspectRatio=3/8, {
    plotMargins(reporter, array=ii, ylim=c(-1,1)*0.2, ylab=ylab, margins="columns", rotate=0);
  });
} # for (ii ...)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Spatial residual plots test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
plm <- res$plm;
print(plm);

# Calculate PLM residuals
rs <- calculateResidualSet(plm, verbose=verbose);
print(rs);

ae <- ArrayExplorer(rs);
setColorMaps(ae, c("log2,log2neg,rainbow", "log2,log2pos,rainbow"));
print(ae);

process(ae, interleaved="auto", verbose=verbose);
