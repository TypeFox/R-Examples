if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chipType <- "Canine_2";

footer <- list(
  createdOn = format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE),
  createdBy = list(
    fullname = "Henrik Bengtsson", 
    email = sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org")
  ),
  srcFiles = list()
);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

ptb <- AffymetrixProbeTabFile$byChipType(chipType);
print(ptb);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allocate aroma cell sequence file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acs <- AromaCellSequenceFile$allocateFromCdf(cdf, tags="HB20100918");
print(acs);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Import data from the Affymetrix probe-tab file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
importFrom(acs, ptb, verbose=log);

# Infer MM from PM sequences?  Will give an error if no MMs
inferMmFromPm(acs, cdf=cdf, verbose=log);
