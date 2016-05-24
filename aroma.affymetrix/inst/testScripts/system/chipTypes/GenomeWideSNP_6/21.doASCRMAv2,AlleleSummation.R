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
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Summarize theta = thetaA+thetaB
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesN <- res$cesN;
as <- AlleleSummation(cesN);
print(as);
cesT <- process(as, verbose=verbose);
print(cesT);

Y <- extractMatrix(cesT, verbose=verbose);
print(head(Y));


# Validation
dsNList <- res$dsNList;
Y0 <- extractMatrix(dsNList$total, verbose=verbose);
print(head(Y0));

# Sanity check
stopifnot(all.equal(Y, Y0));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cbs <- CbsModel(cesT);
print(cbs);

ce <- ChromosomeExplorer(cbs, zooms=2^(0:5));
process(ce, arrays=1:2, chromosomes=18:23, verbose=verbose);
