library("aroma.affymetrix");
verbose <- Verbose(threshold=-4, timestamp=TRUE);

dataSet <- "FusionSDK_Test3";
chipType <- "Test3";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);

cdf <- getCdf(csR);
print(cdf);

ae <- ArrayExplorer(csR);
setColorMaps(ae, "log2,yellow");
process(ae, verbose=verbose);
