## Support for trees with unresolved clades
make.clade.tree <- function(tree, clades) {
  if ( !identical(class(tree), "phylo") )
    stop("tree must be a plain 'phylo' tree")
  if ( !all(names(clades) %in% tree$tip.label) )
    stop("Unknown clade representatives")
  if ( !all(sapply(clades, is.character)) )
    stop("'clades' must be a list of character vectors")
  if ( any(duplicated(unlist(clades))) )
    stop("Duplicated species names")
  
  tree$clades <- clades
  class(tree) <- c("clade.tree", class(tree))
  tree
}

## This takes a state vector and a tree that has been split into
## unresolved clades and produces an 'unresolved' data structure.
make.unresolved.bisse <- function(clades, states) {
  if ( !all(unlist(clades) %in% names(states)) )
    stop("Species in 'clades' do not have states information")

  count.true  <- function(x) sum(!is.na(x) & x)
  count.false <- function(x) sum(!is.na(x) & !x)
  
  states.clades <- lapply(clades, function(x) states[x])
  unresolved <-
    data.frame(tip.label=names(clades),
               Nc=sapply(clades, length),
               n0=sapply(states.clades, count.false),
               n1=sapply(states.clades, count.true))
  rownames(unresolved) <- NULL
  unresolved
}

make.unresolved.bd <- function(clades)
  sapply(clades, length)

clades.from.polytomies <- function(tree) {
  from <- tree$edge[,1]
  to   <- tree$edge[,2]
  n.taxa <- length(tree$tip.label)
  is.node <- seq_len(max(tree$edge)) %in% from

  edge.counts <- tapply(to, from, length)
  poly.nodes <- as.integer(names(edge.counts[edge.counts > 2]))

  ## Find who the ancestors of the polytomy nodes are:
  ## TODO: I think this might be descendants?
  ans <- lapply(poly.nodes, ancestors2, tree)

  ## Now, some of these are nested within one another:
  ans1 <- mapply(setdiff, ans, poly.nodes, SIMPLIFY=FALSE)
  clades <-  lapply(ans1[!(poly.nodes %in% unlist(ans1))],
                    function(x) x[x <= n.taxa])

  ## One representative species (arbitrarily the first species in the
  ## clade) will be left in the phylogeny.
  clades.repr <- tree$tip.label[sapply(clades, "[",  1)]
  names(clades) <- clades.repr
  clades.spp <- lapply(clades, function(x) tree$tip.label[x])

  ## Drop all but the representative species from each clade.
  clades.drop <- sort(unlist(lapply(clades, "[", -1)))
  tree2 <- drop.tip.fixed(tree, clades.drop)

  make.clade.tree(tree2, clades.spp)
}

## Generate a sample using the algorithm in FitzJohn et al. 2009.
clades.from.sample <- function(phy, p) {
  n.taxa <- length(phy$tip.label)

  desc <- lapply(seq_len(phy$Nnode) + n.taxa, get.descendants, phy, TRUE)
  anc <- ancestors(phy)

  f <- function(i) {
    check <- rev(na.omit(anc[,i]))[-1] - n.taxa
    for ( j in check ) {
      n <- sum(keep[desc[[j]]])
      if ( n == 1 )
        return(desc[[j]][keep[desc[[j]]]])
      else if ( n > 1 )
        return(NA)
    }
  }
  
  keep <- runif(n.taxa) < p
  repeat {
    drop <- which(!keep)
    tmp <- sapply(drop, f)
    if ( any(is.na(tmp)) ) {
      orphan <- drop[is.na(tmp)]
      keep[orphan[runif(length(orphan)) < p]] <- TRUE
    } else {
      clades <- split(phy$tip.label[drop], phy$tip.label[tmp])
      break
    }
  }

  make.clade.tree(drop.tip.fixed(phy, drop), clades)
}


## Renamed poorly because of a clash with util.R:ancestors
## TODO: I believe that this is actually descendents()
## It can be replaced by:
##   function(x, tree) descendants(x, tree$edge)
ancestors2 <- function(x, tree, tips.only=FALSE) {
  from <- tree$edge[,1]
  to   <- tree$edge[,2]

  is.node <- seq_len(max(tree$edge)) %in% from

  anc <- list(x)
  n <- 1
  while ( length(x) > 0 ) {
    kids <- to[from %in% x]
    anc[[n <- n + 1]] <- kids
    x <- kids[is.node[kids]]
  }

  anc <- sort(unlist(anc))
  if ( tips.only )
    anc[anc <= length(tree$tip.label)]
  else
    anc
}

## Apply a classification to a tree.  Here 'class' is a character 
clades.from.classification <- function(tree, class, check=TRUE) {
  n.tip <- length(tree$tip.label)

  if ( class(class) != "character" )
    stop("'class' must be a character vector")
  if ( length(class) != n.tip )
    stop("'class' must be the same length as the tip.labels")

  ## Doing this, I need to check that all the groups are
  ## monophyletic.  To do this, I need a MRCA function that takes a
  ## number of tips.
  spp.cl <- split(tree$tip.label, class)
  
  if ( check ) {
    mrca <- sapply(spp.cl, function(x)
                   mrca.tipset(tree, match(x, tree$tip.label)))
    desc <- lapply(mrca, descendants, tree$edge)
    chk <- lapply(desc, function(x) tree$tip.label[x[x <= n.tip]])
    ok <- sapply(seq_along(chk), function(i)
                 length(setdiff(chk[[i]], spp.cl[[i]]))) == 0

    if ( any(!ok) ) {
      stop("Some groups had problems: ",
           paste(names(desc)[!ok], collapse=", "))
      i <- which(!ok)[1]
      to.drop <- setdiff(tree$tip.label, chk[[i]])
      tmp <- drop.tip.fixed(tree, setdiff(tree$tip.label, chk[[i]]))
      col <- (tmp$tip.label %in% tree$tip.label[class == names(chk)[i]])+1
      plot2.phylo(tmp, type="f", cex=.5, label.offset=1, font=1,
                  tip.color=col, no.margin=TRUE)
      tmp$tip.label[tmp$tip.label %in% tree$tip.label[class == names(chk)[i]]]
    }
  }

  to.drop <- unlist(lapply(spp.cl, "[", -1))
  if ( length(to.drop) > 0 )
    tree.cl <- drop.tip.fixed(tree, unlist(lapply(spp.cl, "[", -1)))
  else
    tree.cl <- tree
  tmp <- sapply(spp.cl, "[[", 1)
  tree.cl$tip.label <- names(tmp)[match(tree.cl$tip.label, tmp)]

  make.clade.tree(tree.cl, spp.cl)
}

## This is all that is required for the plotting now; the rest gets
## done by the plot-alt code.
plot.clade.tree <- function(x, as.clade.tree=TRUE, transform=identity,
                            ...) {
  if ( abs(transform(1) - 1) > 1e-8 )
    stop("transform(1) must equal 1")
  if ( as.clade.tree ) {
    n.taxa <- transform(sapply(x$clades, length)[x$tip.label])
    n.taxa[is.na(n.taxa)] <- 1
    names(n.taxa) <- x$tip.label
  } else
    n.taxa <- NULL

  plot2.phylo(x, n.taxa=n.taxa, ...)
}

polytomies.to.clades <- function(tree) {
  .Deprecated("clades.from.polytomies")
  clades.from.polytomies(tree)
}
