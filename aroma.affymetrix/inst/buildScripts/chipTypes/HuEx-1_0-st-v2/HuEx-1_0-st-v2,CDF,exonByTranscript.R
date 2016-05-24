if (interactive()) savehistory();
library("aroma.affymetrix");
verbose <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chipType <- "HuEx-1_0-st-v2";
types <- c("core", "extended", "full");
type <- textMenu(types, value=TRUE);
tags <- "HB20110917";

cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

pattern <- "na31.hg19.probeset.csv$";
csv <- AffymetrixNetAffxCsvFile$byChipType(chipType, pattern=pattern);
print(csv);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Build exon-by-transcript CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdfT <- createExonByTranscriptCdf(cdf, csv=csv, type=type, 
                                  tags=c("*", tags), verbose=verbose);
print(cdfT);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Validation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cells <- getCellIndices(cdfT, unlist=TRUE, useNames=FALSE);
if (anyDuplicated(cells)) {
  throw("Detected cell indices that occurs more than once in the CDF: ", 
                                                      getPathname(cdfT));
}
