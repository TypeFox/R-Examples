library("aroma.affymetrix")

ovars <- ls(all.names=TRUE)
oplan <- future::plan()

## Setup dataset
dataset <- "GSE8605"
chipType <- "Mapping10K_Xba142"

csR <- AffymetrixCelSet$byName(dataset, chipType=chipType)
csR <- csR[1:2]
print(csR)

## Create unique CDF
cdf <- getCdf(csR)
print(cdf)
cdfU <- getUniqueCdf(cdf)
print(cdfU)

checksum <- NULL

strategies <- c("lazy", "eager")
if (future::supportsMulticore()) strategies <- c(strategies, "multicore")
if (packageVersion("future") > "0.10.9") strategies <- c(strategies, "multisession")
if (require("async")) {
  strategies <- c(strategies, "batchjobs")
  async::backend("local")
}

for (strategy in strategies) {
  message(sprintf("*** Using %s futures ...", sQuote(strategy)))

  future::plan(strategy)
  tags <- c("*", strategy)

  ## (a) Process a single array
  csU1 <- convertToUnique(csR[1], tags=c(tags, "one-array"), verbose=verbose)
  print(csU1)
  csU1z <- getChecksumFileSet(csU1)
  print(csU1z[[1]])

  ## Compare file checksum to previous runs
  checksumT <- readChecksum(csU1z[[1]])
  if (is.null(checksum)) checksum <- checksumT
  ## FIXME: File checksums does not work for comparison since
  ## the CEL header has a timestamp of the creation time.
##  stopifnot(identical(checksumT, checksum))


  ## (b) Process two arrays
  csU <- convertToUnique(csR, tags=tags, verbose=verbose)
  print(csU)
  csUz <- getChecksumFileSet(csU)
  print(csU[[1]])
  res <- equals(csU1z[[1]], csUz[[1]])
  ## FIXME: File checksums does not work for comparison since
  ## the CEL header has a timestamp of the creation time.
##  if (!res) throw(res)  ## FIXME

  message(sprintf("*** Using %s futures ... DONE", sQuote(strategy)))
}


## CLEANUP
future::plan(oplan)
rm(list=setdiff(ls(all.names=TRUE), ovars))
