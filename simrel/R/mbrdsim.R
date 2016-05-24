mbrdsim <-function (simlist, fraction) 
{
  nlev <- unlist(lapply(simlist, length))
  l2lev <- log2(nlev)
  repnames <- rep(names(l2lev), times = l2lev)
  ext <- unlist(sapply(l2lev, FUN = function(x) {
    paste(".", 1:x, sep = "")
  }))
  bitnames <- paste(repnames, ext, sep = "")
  res <- mbrd(l2lev, fraction = fraction, fnames2 = names(simlist), 
              fnames1 = bitnames)
  runDesign <- res$Design
  for (i in 1:dim(runDesign)[2]) {
    runDesign[, i] <- simlist[[i]][res$Design[, i]]
  }
  res$Design <- runDesign
  res
}
