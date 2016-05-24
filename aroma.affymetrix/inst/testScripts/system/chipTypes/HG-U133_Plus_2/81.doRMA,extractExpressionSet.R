library("aroma.affymetrix");
library("Biobase")
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "GSE9890";
chipType <- "HG-U133_Plus_2";
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);

ces <- doRMA(csR, verbose=verbose);
print(ces);

eset <- extractExpressionSet(ces, verbose=verbose);
print(eset);
print(sampleNames(eset));

# Sanity checks
stopifnot(identical(sampleNames(eset), getNames(ces)));
stopifnot(identical(featureNames(eset), getUnitNames(getCdf(ces))));


# Add 'annotation' slot
eset <- extractExpressionSet(ces, annotationPkg="PDInfo", verbose=verbose);
print(eset);
print(annotation(eset));
stopifnot(annotation(eset) == "pd.hg.u133.plus.2");

eset <- extractExpressionSet(ces, annotationPkg="cdf", verbose=verbose);
print(eset);
print(annotation(eset));
stopifnot(annotation(eset) == "hgu133plus2cdf");

eset <- extractExpressionSet(ces, annotationPkg="ChipDb", verbose=verbose);
print(eset);
print(annotation(eset));
stopifnot(annotation(eset) == "hgu133plus2");
