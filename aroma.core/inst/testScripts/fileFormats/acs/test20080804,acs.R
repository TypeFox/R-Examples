library("aroma.core");
library("aroma.affymetrix");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Function to draw 'n' random sequences of length 'size'.
rseqs <- function(n=1, size=25) {
  seqs <- sample(c("A", "C", "G", "T"), size=n*size, replace=TRUE);
  dim(seqs) <- c(n, size);
  seqs <- apply(seqs, MARGIN=1, paste, collapse="");
  seqs;
} # rseqs()


log <- Arguments$getVerbose(-20, timestamp=TRUE);

chipType <- "Mapping50K_Hind240";
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allocate ACS file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tags <- c("*", "forTestingOnly");
acs <- AromaCellSequenceFile$allocateFromCdf(cdf, tags=tags, overwrite=TRUE);
print(acs);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Write and read some random sequences
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
n <- 10;
cells <- sample(nbrOfCells(acs), size=n, replace=FALSE);
seqs <- rseqs(n=n);

# A newly allocated ACS file contains all missing sequences
seqs0 <- readSequences(acs, cells=cells);
print(seqs0);
stopifnot(all(is.na(seqs0)));

# Write the test sequences to file
updateSequences(acs, cells=cells, seqs=seqs);

# Assert the correctness of the writing and reading
seqs1 <- readSequences(acs, cells=cells);
print(seqs1);
stopifnot(identical(seqs1, seqs));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Update the ACS file footer with info about the source file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fake that we got the sequences from a source file
pathname <- tempfile();
cat(file=pathname, seqs);
srcFile <- GenericDataFile(pathname);

# Add info about the source file to the ACS footer
footer <- readFooter(acs);
footer$srcFile = list(
  fullname = getFullName(srcFile),
  filesize = getFileSize(srcFile),
  checksum = getChecksum(srcFile)
);
writeFooter(acs, footer);
print(acs);

# Validate written footer
footer2 <- readFooter(acs);
footer$srcFile$filesize <- as.character(footer$srcFile$filesize);
stopifnot(identical(footer2, footer));

file.remove(getPathname(acs));
