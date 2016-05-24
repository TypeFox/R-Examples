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
dataSet <- "GSE34754";
chipType <- "Mapping250K_Nsp";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsNList <- doASCRMAv2(csR);
print(dsNList);

dsN <- exportAromaUnitPscnBinarySet(dsNList);
print(dsN);


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
setFullNamesTranslator(dsN, df);

# Identify unique sample names
sampleNames <- unique(getNames(dsN));

dsList <- lapply(sampleNames, FUN=function(sampleName) {
  ds <- dsN[sampleName];
  lapply(c(T="T", N="N"), FUN=function(type) {
    ds[sapply(ds, hasTag, type)];
  });
});
names(dsList) <- sampleNames;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract tumor-normal pairs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
R <- 5; # Number of replicates
dsList2 <- list();
for (key in names(dsList)) {
  dsListKK <- dsList[[key]];
  dsT <- dsListKK$T;
  dsN <- dsListKK$N;

  # (a) Replicated tumor vs same normal
  dsTa <- dsT[1:min(length(dsT),R)];
  dsNa <- dsN[rep(1L, times=length(dsTa))];

  # (b) Same tumor vs replicated normals
  dsNb <- dsN[1:min(length(dsN),R)];
  dsTb <- dsT[rep(1L, times=length(dsNb))];

  # (c) Merge
  dsT <- append(dsTa, dsTb);
  dsN <- append(dsNa, dsNb);

  # (d) Drop duplicated tumor-normal pairs
  nT <- getFullNames(dsT);
  nN <- getFullNames(dsN);
  nP <- paste(nT, nN, sep="_vs_");
  dups <- duplicated(nP);
  dsT <- dsT[!dups];
  dsN <- dsN[!dups];
  stopifnot(length(dsT) == length(dsN));

  dsList2[[key]] <- list(T=dsT, N=dsN);
} # for kk ...)

dsT <- Reduce(append, lapply(dsList2, FUN=function(x) x$T));
dsN <- Reduce(append, lapply(dsList2, FUN=function(x) x$N));


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
## caller <- PairedPscbsCaller(sms);
## print(caller);
## sms <- process(caller, verbose=verbose);
## print(sms);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Generate report (just to check)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setOption("PSCBS::reports/pscnSegmentationTransitions", TRUE);

# Generate reports for tumor-normal pairs
for (ii in 1:min(length(sms),5)) {
  df <- sms[[ii]];
  fit <- loadObject(df);
  sampleName <- getFullName(df);
  pathname <- report(fit, sampleName=sampleName, studyName=getFullName(dsT), verbose=verbose);
  print(pathname);
} # for (ii ...)
