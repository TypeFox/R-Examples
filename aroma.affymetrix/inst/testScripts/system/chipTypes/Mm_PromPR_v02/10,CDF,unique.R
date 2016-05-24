library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

chipType <- "Mm_PromPR_v02";

cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

# Get the "unique" CDF, which is generated if missing
cdfU <- getUniqueCdf(cdf, verbose=more(verbose, 60));
print(cdfU);

# Note that the "unique" CDF has the same units as the main CDF.
stopifnot(nbrOfUnits(cdfU) == nbrOfUnits(cdf));
