library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "GSE9890";
chipType <- "HG-U133_Plus_2";
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);

res <- doGCRMA(csR, drop=FALSE, verbose=verbose);
print(res);

ces <- doGCRMA(csR, verbose=verbose);
print(ces);
