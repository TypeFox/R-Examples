## functions for generating random trees
## A Hipp, 23 sept 2010 (ahipp@mortonarb.org)

rtreePhylo = function(x, N = 1, ...) {
  ## rtree method that just calls rtree (ape) using an input tree to get taxon names and tips
  ## Arguments:
  ##   x = a tree in phylo format
  ##   N = number of trees to return
  ##   ... = parameters to pass along to rtree
  tips = x$tip.label
  out <- lapply(rep(length(tips), N), rtree, tip.label = tips, ...)
  names(out) <- paste('rtree', 1:N, sep = '')
  class(out) <- 'multiPhylo'
  return(out)
  }

