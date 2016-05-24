library("aroma.affymetrix")
ovars <- ls(all.names=TRUE)
oplan <- future::plan()

message("*** OpticalBackgroundCorrection ...")

dataSet <- "GSE9890"
chipType <- "HG-U133_Plus_2"
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType)
csR <- csR[1:3]
print(csR)

cdf <- getCdf(csR)
acs <- getAromaCellSequenceFile(cdf)
print(acs)

strategies <- c("lazy", "eager")
if (future::supportsMulticore()) strategies <- c(strategies, "multicore")
if (packageVersion("future") > "0.10.9") strategies <- c(strategies, "multisession")
if (require("async")) {
  strategies <- c(strategies, "batchjobs")
  async::backend("local")
}

checksum <- NULL

for (strategy in strategies) {
  message(sprintf("*** Using %s futures ...", sQuote(strategy)))

  future::plan(strategy)
  tags <- c("*", strategy)

  bg <- OpticalBackgroundCorrection(csR, tags=tags)
  print(bg)
  csB <- process(bg, verbose=verbose)
  print(csB)

  csBz <- getChecksumFileSet(csB)
  print(csBz[[1]])
  checksumT <- readChecksum(csBz[[1]])
  if (is.null(checksum)) checksum <- checksumT
  stopifnot(identical(checksumT, checksum))

  message(sprintf("*** Using %s futures ... DONE", sQuote(strategy)))
}

message("*** OpticalBackgroundCorrection ... DONE")

## CLEANUP
future::plan(oplan)
rm(list=setdiff(ls(all.names=TRUE), ovars))
