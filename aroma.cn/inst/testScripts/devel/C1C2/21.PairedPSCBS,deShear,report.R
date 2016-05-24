library("aroma.cn");
library("PSCBS");
library("R.devices");
library("R.menu");
verbose <- Arguments$getVerbose(-10);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup Paired PSCBS segmentation data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "pscbsData";
path <- Arguments$getReadablePath(rootPath);

dataSets <- list.files(rootPath);
if (length(dataSets) > 1) {
 dataSet <- textMenu(dataSets, value=TRUE);
} else {
 dataSet <- dataSets[1];
}

path <- file.path(rootPath, dataSet);
path <- Arguments$getReadablePath(path);
chipTypes <- list.files(path);
if (length(chipTypes) > 1) {
 chipType <- textMenu(chipTypes, value=TRUE);
} else {
 chipType <- chipTypes[1];
}

ds <- PairedPSCBSFileSet$byName(dataSet, chipType=chipType);
print(ds);
dsName <- getName(ds);

if (length(ds) == 0) {
 throw("No PairedPSCBS data file found.")
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select tumor-normal pair
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (length(ds) > 1) {
 ii <- textMenu(getFullNames(ds));
} else {
 ii <- 1L;
}

df <- getFile(ds, ii);
fit <- loadObject(df);
sampleName <- getFullName(df);
sampleName <- gsub(",PairedPSCBS", "", sampleName, fixed=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# REPORT
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tryCatch({
  report(fit, rspTags="deShear", sampleName=sampleName, studyName=dataSet, force=TRUE);
}, error=function(ex) {
  .lastError <<- ex;
});

