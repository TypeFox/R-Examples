############################################################################
# Author: Henrik Bengtsson
# Created on: 2010-05-26
# Last updated: 2010-05-26
#
# USAGE:
# path <- "downloadScripts/chipTypes/Mapping50K_Hind240,Xba240";
# filename <- "HapMap,testSet,100K.R";
# pathname <- system.file(path, filename, package="aroma.affymetrix");
# source(pathname);
#
# DOWNLOADS:
# rawData/
#  HapMap,testSet/
#   Mapping50K_Hind240/
#    CEU_NA06985_HIND.CEL, CEU_NA06991_HIND.CEL, CEU_NA06993_HIND.CEL,
#    CEU_NA06994_HIND.CEL, CEU_NA07000_HIND.CEL, CEU_NA07019_HIND.CEL
#   Mapping50K_Xba240/
#    CEU_NA06985_XBA.CEL, CEU_NA06991_XBA.CEL, CEU_NA06993_XBA.CEL,
#    CEU_NA06994_XBA.CEL, CEU_NA07000_XBA.CEL, CEU_NA07019_XBA.CEL
############################################################################
library("R.filesets");
verbose <- Arguments$getVerbose(-5, timestamp=TRUE);

urlPath <- "http://hapmap.ncbi.nlm.nih.gov/downloads/raw_data/affy100k/";

files <- list(
  "Mapping50K_Hind240" = c(
    "CEU_NA06985_HIND.CEL.gz", 
    "CEU_NA06991_HIND.CEL.gz", 
    "CEU_NA06993_HIND.CEL.gz", 
    "CEU_NA06994_HIND.CEL.gz", 
    "CEU_NA07000_HIND.CEL.gz", 
    "CEU_NA07019_HIND.CEL.gz"
  ),
  "Mapping50K_Xba240" = c(
    "CEU_NA06985_XBA.CEL.gz", 
    "CEU_NA06991_XBA.CEL.gz", 
    "CEU_NA06993_XBA.CEL.gz", 
    "CEU_NA06994_XBA.CEL.gz", 
    "CEU_NA07000_XBA.CEL.gz", 
    "CEU_NA07019_XBA.CEL.gz"
  )
);

rootPath <- "rawData";
dataSet <- "HapMap,testSet";
chipTypes <- names(files);

dsList <- list();
for (kk in seq(along=chipTypes)) {
  chipType <- chipTypes[kk];
  verbose && enter(verbose, sprintf("Chip type #%d ('%s') of %d", 
                                         kk, chipType, length(chipTypes)));

  path <- file.path(rootPath, dataSet, chipType);
  filenames <- files[[chipType]];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Download only files if missing (with or without *.gz)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathnames <- file.path(path, filenames);
  pathnames <- sapply(pathnames, FUN=filePath, expandLinks="any");

  # Gunzip already downloaded *.gz files
  keep <- sapply(pathnames, FUN=isFile);
  ds <- GenericDataFileSet(lapply(pathnames[keep], FUN=GenericDataFile));
  sapply(ds, gunzip, verbose=TRUE);

  # See what is missing
  pathnames <- gsub(".gz", "", pathnames, fixed=TRUE);
  keep <- !sapply(pathnames, FUN=isFile);
  filenames <- filenames[keep];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Download and gunzip
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (length(filenames) > 0) {
    # Download
    urls <- paste(urlPath, filenames, sep="");
    pathnamesZ <- sapply(urls, FUN=downloadFile, path=path, verbose=TRUE);
    verbose && print(verbose, pathnamesZ);

    # Gunzip the data files
    ds <- GenericDataFileSet(lapply(pathnamesZ, FUN=GenericDataFile));
    sapply(ds, gunzip, verbose=TRUE);
    verbose && cat(verbose, "Downloaded:");
    verbose && print(verbose, ds);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup data set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- GenericDataFileSet(lapply(pathnames, FUN=GenericDataFile));
  verbose && cat(verbose, "Complete file set:");
  verbose && print(verbose, ds);

  dsList[[chipType]] <- ds;

  verbose && exit(verbose);
} # for (chipType ...)

verbose && print(verbose, dsList);


############################################################################
# HISTORY:
# 2010-05-26
# o Created.
############################################################################
