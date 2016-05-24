library("aroma.core");
log <- Arguments$getVerbose(-20, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allocate ACF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "FooDataSet";
platform <- "Generic";
chipType <- "FooBar";
nbrOfRows <- 402;
path <- filePath("callData", dataSet, chipType);

for (kk in 1:6) {
  filename <- sprintf("foo%02d.acf", kk);
  acf <- AromaUnitCallFile$allocate(filename, path=path, nbrOfRows=nbrOfRows, platform=platform, chipType=chipType);
  acf2 <- AromaUnitCallFile(filename, path=path);
  stopifnot(equals(acf2, acf));
}

acs <- AromaUnitCallSet$byName(dataSet, chipType=chipType);
