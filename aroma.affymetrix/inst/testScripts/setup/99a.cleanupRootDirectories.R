library("R.utils");
if (interactive()) library("R.menu");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);


verbose && enter(verbose, "Removing all generated data");

verbose && cat(verbose, "This will *not* remove annotation data or raw data sets.");

# Root directories to be removed, iff existing
paths <- c(
  "probeData/",
  "plmData/",
  "qcData/",
  "firmaData/",
  "cbsData/",
  "gladData/",
  "haarData/",
  "totalAndFracBData/",
  "totalAndFracBData,txt/",
  "oligoData/",
  "crlmmData/",
  "Data/",
  "rawCnData/",
  "smoothCnData/",
  "reports/",
  "figures/"
);

# Identify which exists
keep <- sapply(paths, FUN=isDirectory);
paths <- paths[keep];

if (length(paths) > 0) {
  if (interactive()) {
    paths <- selectMenu(paths, selected=TRUE, title="Root directories to be removed");
    verbose && cat(verbose, "Root directories to be removed:");
    verbose && cat(verbose, paste(paths, collapse="\n"));
  }
  
  
  for (kk in seq_along(paths)) {
    path <- paths[kk];
    verbose && enter(verbose, sprintf("Root directory #%s ('%s') of %s", kk, path, length(paths)));
    removeDirectory(path, recursive=TRUE);
    verbose && exit(verbose);
  } # for (kk ...)
} else {
  verbose && cat(verbose, "No root directories available. Skipping!");
}

verbose && exit(verbose);
