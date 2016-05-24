library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE13372,testset";
chipType <- "GenomeWideSNP_6,Full";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (a) AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
dsNList0 <- lapply(dsNList, FUN=`[`, 1);
for (key in names(dsNList1)) {
  verbose && enter(verbose, key);

  dsN0 <- dsNList0[[key]];
  verbose && print(verbose, dsN0);
  dsN1 <- dsNList1[[key]];
  verbose && print(verbose, dsN1);

  stopifnot(getFullNames(dsN1) == getFullNames(dsN0));
  data0 <- extractMatrix(dsN0);
  data1 <- extractMatrix(dsN1);

  rho <- cor(data1, data0, use="complete.obs");
  verbose && print(verbose, rho);

  res <- all.equal(data1, data0);
  verbose && print(verbose, res);

  # Sanity check
  stopifnot(all.equal(data1, data0));

  verbose && exit(verbose);
} # for (key ...)
