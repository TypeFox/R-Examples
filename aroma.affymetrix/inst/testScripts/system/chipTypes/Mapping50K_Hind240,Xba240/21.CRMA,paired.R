library("aroma.affymetrix")
log <- Arguments$getVerbose(-4, timestamp=TRUE);


dataSet <- "HapMap,CEU,testset";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");
#chipTypes <- chipTypes[2];

# Expected sample names
sampleNames <- c("NA06985", "NA06991", "NA06993", 
                 "NA06994", "NA07000", "NA07019");
tags <- "ACC,-XY,RMA,+300,A+B,FLN,-XY";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesList <- list();
for (chipType in chipTypes) {
  ces <- CnChipEffectSet$byName(dataSet, tags=tags, chipType=chipType);
  print(ces);
  stopifnot(identical(getNames(ces), sampleNames));
  cesList[[chipType]] <- ces;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired GLAD model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired CN anlysis, by setting up some fake (test,control) pairs
testList <- lapply(cesList, FUN=extract, 1:3);
refList <- lapply(cesList, FUN=extract, 4:6);
glad <- GladModel(testList, refList);
print(glad);

fit(glad, arrays=1, chromosomes=19, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ChromosomeExplorer test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ce <- ChromosomeExplorer(glad);
print(ce);
process(ce, arrays=1:2, chromosomes=c(19,22), verbose=log);
## process(ce, verbose=log);
