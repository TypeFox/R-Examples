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
# Background correction and normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdfS <- AffymetrixCdfFile$byChipType(getChipType(cdf, fullname=FALSE));
setCdf(csR, cdfS);
bc <- GcRmaBackgroundCorrection(csR, type="affinities");
print(bc);

csB <- process(bc, verbose=verbose);
print(csB);
setCdf(csB, cdf);
