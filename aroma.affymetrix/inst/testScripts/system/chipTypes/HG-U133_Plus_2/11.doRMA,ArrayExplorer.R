library("aroma.affymetrix")
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

dataSet <- "GSE9890"
chipType <- "HG-U133_Plus_2"
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType)

res <- doRMA(csR, drop=FALSE, verbose=verbose)
print(res)

plm <- res$plm
rs <- calculateResidualSet(plm, verbose=verbose)
print(rs)

ae <- ArrayExplorer(rs)
setColorMaps(ae, c("log2,log2neg,rainbow", "log2,log2pos,rainbow"))
print(ae)

process(ae, arrays=1:2, interleaved="auto", verbose=verbose)
