library("aroma.affymetrix");
verbose <- Verbose(threshold=-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "Affymetrix-HapMap,CEU,set1";
chipType <- "GenomeWideSNP_5";

cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full,r2");
print(cdf);

csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);

# Split up into names and tags
setFullNamesTranslator(csR, function(names, ...) {
  gsub("^(NA[0-9]+)_(.*)", "\\1,\\2", names);
});

# Drop duplicates
keep <- !duplicated(getNames(csR));
csR <- csR[keep];

# Keep only certain samples
sampleNames <- c("NA06985", "NA06991", "NA06993",
                 "NA07019", "NA07022", "NA07056");
csR <- csR[sampleNames];

print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsNList <- res$dsNList;

seg <- CbsModel(dsNList$total);
print(seg);

ce <- ChromosomeExplorer(seg);
print(ce);
process(ce, arrays=1, chromosomes=c(19,23), verbose=verbose);
