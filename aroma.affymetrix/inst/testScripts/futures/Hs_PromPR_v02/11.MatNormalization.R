library("aroma.affymetrix")
ovars <- ls(all.names=TRUE)
oplan <- future::plan()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE24546,testset"
chipType <- "Hs_PromPR_v02"

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType)
csR <- csR[1:2]
print(csR)

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
  mn <- MatNormalization(csR[1], tags=c(tags, "one-array"))
  print(mn)
  csM1 <- process(mn, verbose=-10)
  print(csM1)
  csM1z <- getChecksumFileSet(csM1)
  print(csM1z[[1]])
  checksumT <- readChecksum(csM1z[[1]])
  if (is.null(checksum)) checksum <- checksumT
  stopifnot(identical(checksumT, checksum))

  ## (b) Process two arrays
  mn <- MatNormalization(csR, tags=tags)
  print(mn)
  csM <- process(mn, verbose=-10)
  print(csM)
  csMz <- getChecksumFileSet(csM)
  print(csMz[[1]])
  res <- equals(csM1z[[1]], csMz[[1]])
  if (!res) throw(res)

  message(sprintf("*** Using %s futures ... DONE", sQuote(strategy)))
}

## CLEANUP
future::plan(oplan)
rm(list=setdiff(ls(all.names=TRUE), ovars))
