##########################################################################
# AS-CRMAv2 and Paired PSCBS
##########################################################################
library("aroma.affymetrix");
library("aroma.cn");  # PairedPscbsModel
stopifnot(packageVersion("aroma.cn") >= "1.2.4");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE12702";
chipType <- "Mapping250K_Nsp";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsNList <- doASCRMAv2(csR, verbose=verbose);
print(dsNList);

dsN <- exportAromaUnitPscnBinarySet(dsNList);
print(dsN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CalMaTe
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cmt <- CalMaTeCalibration(dsNList);
print(cmt);

dsCList <- process(cmt, verbose=verbose);
print(dsCList);

dsC <- exportAromaUnitPscnBinarySet(dsCList);
print(dsC);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Group samples by name and type
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AD HOC: For now, just hardwire the path.
path <- file.path("testScripts/complete/dataSets", dataSet);
db <- TabularTextFile(sprintf("%s,samples.txt", dataSet), path=path);
setColumnNameTranslator(db, function(names, ...) {
  names <- gsub("id", "fixed", names);
  names <- gsub("fullname", "replacement", names);
  names;
});
df <- readDataFrame(db, colClasses=c("*"="character"));
setFullNamesTranslator(dsC, df);

# Identify unique sample names
sampleNames <- unique(getNames(dsC));

dsList <- lapply(sampleNames, FUN=function(sampleName) {
  ds <- dsC[sampleName];
  lapply(c(T="T", N="N"), FUN=function(type) {
    ds[sapply(ds, hasTag, type)];
  });
});
names(dsList) <- sampleNames;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract (single) tumor-normal pairs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dfTList <- lapply(dsList, FUN=function(dsList) { dsList$T[[1]] });
dsT <- newInstance(dsList[[1]]$T, dfTList);
dfNList <- lapply(dsList, FUN=function(dsList) { dsList$N[[1]] });
dsN <- newInstance(dsList[[1]]$T, dfNList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segment tumor-normal pairs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sm <- PairedPscbsModel(dsT=dsT, dsN=dsN, gapMinLength=2e6);
print(sm);

res <- fit(sm, verbose=verbose);
print(res);

sms <- getOutputDataSet(sm);
print(sms);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Call segments to be in ROH, AB and LOH.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
caller <- PairedPscbsCaller(sms);
print(caller);
scs <- process(caller, verbose=verbose);
print(scs);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Generate report (just to check)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setOption("PSCBS::reports/pscnSegmentationTransitions", TRUE);

# Generate reports for tumor-normal pairs
for (ii in 1:min(length(scs),5)) {
  df <- scs[[ii]];
  fit <- loadObject(df);
  pathname <- report(fit, studyName=getFullName(dsT), verbose=verbose);
  print(pathname);
} # for (ii ...)
