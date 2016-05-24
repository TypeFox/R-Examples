library("aroma.core");
library("aroma.affymetrix");
log <- Arguments$getVerbose(-20, timestamp=TRUE);


rootPath <- "callData/";
rootPath <- Arguments$getWritablePath(rootPath);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
unf <- AffymetrixCdfFile$byChipType("Mapping10K_Xba142");
print(unf);

dataSet <- "MyDataSet";
chipType <- getChipType(unf, fullname=FALSE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allocate call file(s)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path(rootPath, dataSet, chipType);

kk <- 1;
filename <- sprintf("foo%02d,genotypes.acf", kk);
acf <- AromaUnitGenotypeCallFile$allocateFromUnitNamesFile(unf, filename=filename, path=path, footer=list(comment="Hello world!"), overwrite=TRUE);
print(acf);

ftr <- readFooter(acf);
str(ftr);
ftr$comment2 <- "Hello moon!";
print(ftr);

writeFooter(acf, ftr);
print(acf);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup call set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acs <- AromaUnitGenotypeCallSet$byName(dataSet, chipType=chipType);
print(acs);

acf <- acs[[1]];
print(acf);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Update call file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acf[1:10,1] <- 1;
acf[1:10,2] <- 1;
acf[10,2] <- 2;
print(extractCalls(acf, units=1:12));
print(extractGenotypes(acf, units=1:12));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allocate confidence score files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (kk in seq_along(acs)) {
  acf <- acs[[kk]];
  fullname <- getFullName(acf);
  fullname <- gsub(",genotypes", ",confidenceScores", fullname);

  filename <- sprintf("%s.acf", fullname);
  pathname <- filePath(getPath(acf), filename);
  pathname <- Arguments$getWritablePathname(pathname, mustNotExist=FALSE);

  asf <- AromaUnitSignalBinaryFile$allocateFromUnitNamesFile(unf, filename=pathname, types="double", size=4, signed=TRUE, footer=list(comment="Hello world!"), overwrite=TRUE);
  asf[,1] <- NA_integer_;
}

ass <- AromaUnitSignalBinarySet$byName(dataSet, chipType=chipType, pattern=",confidenceScores.acf$", paths=rootPath);
print(ass);
asf <- ass[[1]];
print(asf);
asf[1:10,1] <- seq(from=0, to=1, length.out=10);
print(asf[1:12,1]);
