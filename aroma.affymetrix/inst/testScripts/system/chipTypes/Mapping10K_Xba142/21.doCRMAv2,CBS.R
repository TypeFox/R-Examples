library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
dsT <- doCRMAv2(dataSet, chipType=chipType, verbose=verbose);
print(dsT);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CBS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
segB <- CbsModel(dsT);
print(segB);

# Try to segment
fit(segB, arrays=1:2, chromosomes=19, verbose=verbose);
