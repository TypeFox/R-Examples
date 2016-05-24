library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-3, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "Affymetrix-HeartBrain";
chipType <- "HuEx-1_0-st-v2";
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="coreR3,A20071112,EP");
print(cdf);

# Setup CEL set using the core CDF.
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process only cerebellum and heart
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
types <- c("cerebellum", "heart");
csR <- csR[indexOf(csR, patterns=types)];

setFullName(csR, sprintf("%s,%s", dataSet, paste(types, collapse="+")));
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# FIRMA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doFIRMA(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# PLMs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plmList <- list(
  merge   = res$plm,                               # all exons together
  noMerge = ExonRmaPlm(res$csN, mergeGroups=FALSE) # each exon separately
);
print(plmList);

# Fit the per-exon PLM
fit(plmList$noMerge, verbose=verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# PLM residuals (also calculated by FIRMA)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
resList <- lapply(plmList, FUN=calculateResidualSet, verbose=verbose);
print(resList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# PLM weights
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
weightList <- lapply(plmList, FUN=calculateWeights, verbose=verbose);
print(weightList);
