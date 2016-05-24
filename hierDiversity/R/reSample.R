reSample <-
function(dat, cluster, replace) {
  if (nrow(dat) == 1 || all(replace == FALSE)) 
    return(dat)
  cls <- sample(unique(dat[[cluster[1]]]), replace = replace[1])
  sub <- vector("list", length(cls))
  ind <- 1
  for (b in cls) {
    x <- dat[dat[[cluster[[1]]]] == b, ]
    if (sum(cls == b) > 1) {
      for (B in cluster) {
        x[[B]] <- paste(x[[B]], ind, sep = ".")
      }
    }
    sub[[ind]] <- x
    ind <- ind + 1
  }
  if (length(cluster) > 1) 
    sub <- lapply(sub, reSample, cluster = cluster[-1], 
      replace = replace[-1])
  do.call(rbind, sub)
}
