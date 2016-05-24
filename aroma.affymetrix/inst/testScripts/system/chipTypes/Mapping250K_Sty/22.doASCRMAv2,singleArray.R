library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE12019,testset";
chipType <- "Mapping250K_Sty";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (a) CRMAv2 - six arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsNList <- doASCRMAv2(csR, verbose=verbose);
print(dsNList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (b) CRMAv2 - single array
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup a single array
csR1 <- csR[1];

# Rename data set in order to not pick up existing results
setFullName(csR1, sprintf("%s,singleArray", getFullName(csR1)));

dsNList1 <- doASCRMAv2(csR1, verbose=verbose);
print(dsNList1);


# Sanity checks
dsNList0 <- lapply(dsNList, FUN=extract, 1);
for (key in names(dsNList1)) {
  print(key);
  dsN0 <- dsNList0[[key]];
  dsN1 <- dsNList1[[key]];
  stopifnot(getFullNames(dsN1) == getFullNames(dsN0));
  data0 <- extractMatrix(dsN0);
  data1 <- extractMatrix(dsN1);
  rho <- cor(data1, data0, use="complete.obs");
  print(rho);
  res <- all.equal(data1, data0);
  print(res);

  # Sanity check
  stopifnot(all.equal(data1, data0));
} # for (key ...)
