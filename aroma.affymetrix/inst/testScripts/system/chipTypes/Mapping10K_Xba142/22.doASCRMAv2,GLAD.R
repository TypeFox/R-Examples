library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
dsNList <- doASCRMAv2(dataSet, chipType=chipType, verbose=verbose);
print(dsNList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CBS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
seg <- GladModel(dsNList$total);
print(seg);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ChromsomeExplorer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ce <- ChromosomeExplorer(seg, zooms=2^(0:5));
print(ce);
process(ce, arrays=1:2, chromosomes=c(1:2,19,23), verbose=verbose);
