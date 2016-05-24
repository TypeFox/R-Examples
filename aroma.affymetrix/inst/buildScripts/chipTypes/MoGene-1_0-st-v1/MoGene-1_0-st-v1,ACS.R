if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);

email <- sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chipType <- "MoGene-1_0-st-v1";

footer <- list(
  createdOn = format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE),
  createdBy = list(
    fullname = "Henrik Bengtsson", 
    email = email
  ),
  srcFiles = list()
);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="r3");
print(cdf);

ptb <- AffymetrixProbeTabFile$byChipType(chipType);
print(ptb);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allocate aroma cell sequence file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acs <- AromaCellSequenceFile$allocateFromCdf(cdf, tags="*,HB20110325");
print(acs);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Import data from the Affymetrix probe-tab file (contains only PMs)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
importFrom(acs, ptb, keepSequenceLengths=25, verbose=log);

footer <- readFooter(acs);
footer$createdBy = list(
  fullname = "Henrik Bengtsson", 
  email = sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org")
);
writeFooter(acs, footer);
