library(R.utils);

if (FALSE) {
  ## Define CN regions:
  ## source("05.defineCopyNumberSegments.R")
  source("trunk/inst/testScripts/system/preprocessing/GSE29172/05.defineCopyNumberSegments.R")
  str(regDat)
}

regPath <- "cnRegionData";
regPath <- Arguments$getReadablePath(regPath);

dataSet <- "GSE29172,ASCRMAv2"
chipType <- "GenomeWideSNP_6"

path <- file.path(regPath, dataSet, chipType);
path <- Arguments$getReadablePath(path);

sampleName <- "H1395vsBL1395"
pattern <- sprintf("%s,([0-9]+),\\(([0-9]),([0-9])\\).xdr", sampleName)
filenames <- list.files(path, pattern=pattern)
pcts <- unique(gsub(pattern, "\\1", filenames))

savPath <- file.path("extdata", "cnData", chipType)
savPath <- Arguments$getWritablePath(savPath);
  
types <- regDat[["type"]]
for (pct in pcts) {
  print(pct)
  pattern <- sprintf("H1395vsBL1395,%s,(.*).xdr", pct)
  filenames <- list.files(path, pattern=pattern)
  types <- gsub(pattern, "\\1", filenames)
  
  pathnames <- file.path(path, filenames)
  names(pathnames) <- types
  datList <- lapply(pathnames, loadObject)
  dat <- NULL
  for (dd in names(datList)) {
    datDD <- datList[[dd]]
    datDD$region <- dd
    dat <- rbind(dat, datDD)
  }
  str(dat)

  filename <- sprintf("%s,%s,%s,cnRegions.xdr", dataSet, sampleName, pct)
  pathname <- file.path(savPath, filename)
  saveObject(dat, file=pathname)
}


