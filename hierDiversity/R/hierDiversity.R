hierDiversity <-
function(dat, group, replace = NULL, 
  reps = 99, q = 1, quant = c(0.025, 0.975), sims = FALSE) {
  if (is.null(colnames(group)) == TRUE) {
    colnames(group) <- paste("lev", 1:dim(group)[2], 
      sep = "")
  }
  cluster <- rev(colnames(group))
  if (is.null(replace) == TRUE) {
    replace <- rep(TRUE, length(cluster))
  }
  n.levels <- length(cluster)
  out <- vector("list", length = n.levels)
  names(out) <- cluster
  gD <- data.frame(group, stringsAsFactors = FALSE)
  for (i in 1:n.levels) {
    pb <- txtProgressBar(min=0, max=n.levels, style = 3)
    setTxtProgressBar(pb, i)
    next.lev <- cluster[i]
    n.cat <- length(unique(gD[[next.lev]]))
    nlist <- vector("list", length = n.cat)
    names(nlist) <- unique(gD[[next.lev]])
    for (j in 1:n.cat) {
      X <- unique(gD[[next.lev]])[j]
      ndat <- matrix(dat[gD[[next.lev]] == X, 
        ], nrow = sum(gD[[next.lev]] == X))
      if (i == 1) {
        ngroup <- group
        nclust <- cluster
      }
      else {
        ngroup <- group[gD[[next.lev]] == X, 
          1:(n.levels - i + 1)]
        nclust <- cluster[-(1:i - 1)]
      }
      nlist[[j]] <- booties(dat = ndat, group = ngroup, 
        cluster = nclust, replace, reps, q, 
        quant = quant, sims)
    }
    out[[i]] <- nlist
  }
  close(pb)
  return(out)
}
