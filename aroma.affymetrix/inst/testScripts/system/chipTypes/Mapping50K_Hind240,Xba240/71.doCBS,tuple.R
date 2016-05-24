library("aroma.affymetrix")

verbose <- Verbose(threshold=-4, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# doCBS() with explicit data set tuple
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "HapMap,CEU,testset";
tags <- "ACC,-XY,RMA,+300,A+B,FLN,-XY";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

nbrOfSets <- length(chipTypes);
dsList <- vector("list", nbrOfSets);
for (kk in seq_len(nbrOfSets)) {
  chipType <- chipTypes[kk];
  ds <- CnChipEffectSet$byName(dataSet, tags=tags, chipType=chipType,
                              mergeStrands=TRUE, combineAlleles=TRUE);
  dsList[[kk]] <- ds;
}
print(dsList);

dsTuple <- as.CopyNumberDataSetTuple(dsList);
print(dsTuple);

res <- doCBS(dsTuple, arrays=1:2, chromosomes=c(19,21), verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# doCBS() with data set tuple names
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "HapMap270,100K,CEU,testSet";
tags <- "ACC,-XY,RMA,+300,A+B,FLN,-XY";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

res <- doCBS(dataSet, tags=tags, chipTypes=chipTypes,
             arrays=1:2, chromosomes=c(19,21), verbose=verbose);
print(res);
