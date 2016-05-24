library("aroma.affymetrix");

verbose <- Verbose(threshold=-4, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A tuple
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "HapMap,CEU,testset";
chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

dsList <- lapply(chipTypes, FUN=function(chipType) {
  AffymetrixCelSet$byName(dataSet, chipType=chipType);
});
print(dsList);
dsTuple <- AffymetrixCelSetTuple(dsList);
dsTuple <- GenericDataFileSetList(dsList);
print(dsTuple);
print(getNames(dsTuple));  # HMM... == NULL?!?


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Same tuple with fullname translators
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (ds in getSets(dsTuple)) {
  fnts <- getAromaFullNameTranslatorSet(ds);
  appendFullNamesTranslator(ds, fnts);
}
print(dsTuple);
print(getNames(dsTuple));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Attributes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds <- dsList[[1]];
n23 <- sapply(ds, getAttribute, "n23");
n24 <- sapply(ds, getAttribute, "n24");
isFemale <- (n23 == 2 & n24 == 0);
dsXX <- ds[isFemale];
print(dsXX);
