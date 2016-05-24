library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
dsNList <- doASCRMAv2(dataSet, chipType=chipType, verbose=verbose);
print(dsNList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CBS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsT <- dsNList$total;
seg <- CbsModel(dsT);
print(seg);

# Try to segment
fit(seg, arrays=c(1,2), chromosomes=19, verbose=verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ChromsomeExplorer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ce <- ChromosomeExplorer(seg, zooms=2^(0:5));
print(ce);
process(ce, arrays=1:2, chromosomes=c(1:2,19,23), verbose=verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CBS with a common reference
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Use the robust average of the first five arrays as a reference
dsR <- dsT[1:5];
dfR <- getAverageFile(dsR);
print(dfR);

seg <- CbsModel(dsT, ref=dfR, tags=c(getTags(dsT), "commonRef"));
print(seg);

# Try to segment
fit(seg, arrays=c(1,2), chromosomes=19, verbose=verbose);

# Write to file
pathname <- writeRegions(seg, arrays=1, oneFile=FALSE, verbose=verbose);
print(pathname);
