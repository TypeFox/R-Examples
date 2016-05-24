library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE20403";
chipType <- "MoEx-1_0-st-v1";

cdf <- AffymetrixCdfFile$byChipType(chipType, tags="coreR1,A20080718,MR");
print(cdf);

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
csR <- csR[1:6];
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# gcRMA-style background correction
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Currently, you must use the standard CDF file.
cdf <- getCdf(csR);
cdfS <- AffymetrixCdfFile$byChipType(getChipType(cdf, fullname=FALSE));
setCdf(csR, cdfS);
bc <- GcRmaBackgroundCorrection(csR, type="affinities");
print(bc);
csB <- process(bc, verbose=verbose);
print(csB);
# Now, use the custom CDF in what follows
setCdf(csB, cdf);
print(csB);
