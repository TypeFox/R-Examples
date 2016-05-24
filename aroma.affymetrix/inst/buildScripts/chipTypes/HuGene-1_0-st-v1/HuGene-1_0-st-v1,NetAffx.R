if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chipType <- "HuGene-1_0-st-v1";

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
tags <- ".na31.hg19.transcript";
#tags <- ".na31.hg19.probeset";
pattern <- sprintf("^%s%s([.]|_)(annot[.]|)(csv|CSV)$", chipType, tags);
nax <- AffymetrixNetAffxCsvFile$byChipType(chipType, tags=tags, pattern=pattern);
nbrOfLines(nax); # Count number of lines (once!)
print(nax);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Read data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- readGeneAssignments(nax, verbose=-10);

