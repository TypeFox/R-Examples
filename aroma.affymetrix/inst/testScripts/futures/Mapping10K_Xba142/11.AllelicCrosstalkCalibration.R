library("aroma.affymetrix")

ovars <- ls(all.names=TRUE)
oplan <- future::plan()

## Setup dataset
dataset <- "GSE8605"
chipType <- "Mapping10K_Xba142"

csR <- AffymetrixCelSet$byName(dataset, chipType=chipType)
csR <- csR[1:2]
print(csR)

checksum <- NULL

strategies <- list("lazy", "eager")
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
  acc <- AllelicCrosstalkCalibration(csR[1], model="CRMAv2", tags=c(tags, "one-array"))
  print(acc)
  csC1 <- process(acc, verbose=-10)
  print(csC1)
  csC1z <- getChecksumFileSet(csC1)
  print(csC1z[[1]])

  ## Compare file checksum to previous runs
  checksumT <- readChecksum(csC1z[[1]])
  if (is.null(checksum)) checksum <- checksumT
  stopifnot(identical(checksumT, checksum))


  ## (b) Process two arrays
  acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2", tags=tags)
  print(acc)
  csC <- process(acc, verbose=-10)
  print(csC)
  csCz <- getChecksumFileSet(csC)
  print(csCz[[1]])
  res <- equals(csC1z[[1]], csCz[[1]])
  if (!res) throw(res)

  message(sprintf("*** Using %s futures ... DONE", sQuote(strategy)))
}


## CLEANUP
future::plan(oplan)
rm(list=setdiff(ls(all.names=TRUE), ovars))
