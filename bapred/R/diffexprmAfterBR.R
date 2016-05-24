diffexprmAfterBR <-
function(x, xbrlist, y, batch, batchessuited) {

  pvalbatches <- sapply(levels(batch)[batchessuited], function(x1) apply(x[batch==x1,], 2, function(x2) fuzzywilcox(x2, y[batch==x1])))

  constantvarlist <- lapply(xbrlist, function(y) which(apply(y, 2, function(yvar) length(unique(yvar))==1)))
  pvalleavebatchout <- mapply(function(x1, xb, cvar) {
    if(length(cvar)>0) {
      pvals <- rep(NA, ncol(x1))
      pvals[-cvar] <- apply(x1[,-cvar], 2, function(x2) fuzzywilcox(x2, y[batch!=xb]))
      pvals[cvar] <- 1
    }
    else
      pvals <- apply(x1, 2, function(x2) fuzzywilcox(x2, y[batch!=xb]))
    return(pvals)
  }, xbrlist, levels(batch)[batchessuited], constantvarlist)

  signbatches <- apply(pvalbatches, 2, function(y) y <= sort(y)[floor(ncol(x)*0.05)])
  signleavebatchout <- apply(pvalleavebatchout, 2, function(y) y <= sort(y)[floor(ncol(x)*0.05)])

  freqcommon <- apply(signbatches + signleavebatchout, 2, function(y) sum(y==2))/floor(ncol(x)*0.05)
  batchtab <- sapply(1:length(levels(batch)), function(y) sum(batch==y))

  sum(batchtab[batchessuited]*freqcommon)/sum(batchtab[batchessuited])

}
