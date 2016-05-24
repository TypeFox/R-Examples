library("aroma.affymetrix");

##########################################################################
# Data set:
# FusionSDK_Test3/
#  Test3/
##########################################################################
dataSet <- "FusionSDK_Test3";
chipType <- "Test3";

path <- file.path("rawData", dataSet, chipType);
path <- Arguments$getWritablePath(path);
ds <- AffymetrixCelSet$byPath(path);
if (nbrOfFiles(ds) == 0) {
  pathT <- system.file(path, "2.Calvin", package="AffymetrixDataTestFiles");
  pathnames <- list.files(pathT, pattern="[.]CEL$", full.names=TRUE);
  for (ii in seq_along(pathnames)) {
    pathname <- pathnames[ii];
    pathnameD <- file.path(path, basename(pathname));
    copyFile(pathname, pathnameD);
  }
  ds <- AffymetrixCelSet$byPath(path);
}
print(ds);
## AffymetrixCelSet:
## Name: FusionSDK_Test3
## Tags:
## Path: rawData/FusionSDK_Test3/Test3
## Platform: Affymetrix
## Chip type: Test3
## Number of arrays: 2
## Names: Test3-1-121502, Test3-2-121502 [2]
## Time period: 2001-08-16 17:28:31 -- 2001-08-16 17:32:06
## Total file size: 0.31MB
## RAM: 0.01MB
