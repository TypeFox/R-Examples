## I need to tidy this up a little bit to allow for different tree
## types.  I also cannot use functions beginning with simulate(), as
## this is a standard R generic function.
##
## It might be useful to have the simulations use somthing like the
## equilibrium distribution for characters at the bottom of the tree.
##
## It is also worth noting that Luke Harmon has a birthdeath.tree
## function that simulates a tree under a birth death process in
## geiger.

## There is also some missing logic with single branch trees, that I
## need to work through.

## Main interface.  In the hope that I will make this generic over a
## 'model' object, I will design the calling structure in a way that
## is similar to S3 generics/methods.
trees <- function(pars,
                  type=c("bisse", "bisseness", "bd", "classe", "geosse",
                         "musse", "quasse", "yule"), n=1,
                  max.taxa=Inf, max.t=Inf, include.extinct=FALSE,
                  ...) {
  if ( is.infinite(max.taxa) && is.infinite(max.t) )
    stop("At least one of max.taxa and max.t must be finite")
  type <- match.arg(type)
  f <- switch(type,
              bisse=tree.bisse,
              bisseness=tree.bisseness,
              bd=tree.bd,
              classe=tree.classe,
              geosse=tree.geosse,
              musse=tree.musse,
              yule=tree.yule)
  trees <- vector("list", n)  
  i <- 1

  while ( i <= n ) {
    trees[[i]] <- phy <-
      f(pars, max.taxa, max.t, include.extinct, ...)
    if ( include.extinct || !is.null(phy) )
      i <- i + 1
  }

  trees
}

## My dodgy tree format to ape's tree format.  The 'bisse' suffix is
## only here for historical reasons.
me.to.ape.bisse <- function(x, root.state) {
  if ( nrow(x) == 0 )
    return(NULL)
  Nnode <- sum(!x$split) - 1
  n.tips <- sum(!x$split)

  x$idx2 <- NA
  x$idx2[!x$split] <- 1:n.tips
  x$idx2[ x$split] <- order(x$idx[x$split]) + n.tips + 1

  i <- match(x$parent, x$idx)
  x$parent2 <- x$idx2[i]
  x$parent2[is.na(x$parent2)] <- n.tips + 1

  tip.label <- ifelse(subset(x, !split)$extinct,
                      sprintf("ex%d", 1:n.tips),
                      sprintf("sp%d", 1:n.tips))
  node.label <- sprintf("nd%d", 1:Nnode)

  x$name <- NA
  x$name[!x$split] <- tip.label
  ## More useful, but I don't want to clobber anything...
  x$name2 <- c(tip.label, node.label)[x$idx2]

  tip.state <- x$state[match(1:n.tips, x$idx2)]
  names(tip.state) <- tip.label

  node.state <- x$state[match(1:Nnode + n.tips, x$idx2)]
  names(node.state) <- node.label
  node.state["nd1"] <- root.state

  hist <- attr(x, "hist")
  if ( !is.null(hist) ) {
    hist$idx2 <- x$idx2[match(hist$idx, x$idx)]
    hist$name2 <- x$name2[match(hist$idx, x$idx)]
    if ( nrow(hist) > 0 )
      hist <- hist[order(hist$idx2),]
  }

  phy <- reorder(structure(list(edge=cbind(x$parent2, x$idx2),
                                Nnode=Nnode,
                                tip.label=tip.label,
                                tip.state=tip.state,
                                node.label=node.label,
                                node.state=node.state,
                                edge.length=x$len,
                                orig=x,
                                hist=hist),
                           class="phylo"))
  phy$edge.state <- x$state[match(phy$edge[,2], x$idx2)]
  phy
}

## Adapted from prune.extinct.taxa() in geiger
prune <- function(phy, to.drop=NULL) {
  if ( is.null(to.drop) )
    to.drop <- subset(phy$orig, !split)$extinct
  if ( sum(!to.drop) < 2 ) {
    NULL
  } else if ( any(to.drop) ) {
    phy2 <- drop.tip.fixed(phy, phy$tip.label[to.drop])
    ## phy2$orig <- subset(phy2$orig, !extinct) # Check NOTE
    phy2$orig <- phy2$orig[!phy2$orig$extinct,]
    phy2$tip.state <- phy2$tip.state[!to.drop]
    phy2$node.state <- phy2$node.state[phy2$node.label]
    phy2$hist <- prune.hist(phy, phy2)
    phy2
  } else {
    phy
  }
}

## This function aims to covert the "hist" object.  This is fairly
## complicated and possibly can be streamlined a bit.  The big issue
## here is that when extinct species are removed from the tree, it
## leaves unbranched nodes - the history along a branch with such a
## node needs to be joined.
prune.hist <- function(phy, phy2) {
  hist <- phy$hist
  if ( is.null(hist) || nrow(hist) == 0 )
    return(hist)

  ## More interesting is to collect up all of the names and look at the
  ## branches that terminate
  phy.names <- c(phy$tip.label, phy$node.label)
  phy2.names <- c(phy2$tip.label, phy2$node.label)

  ## Next, check what the parent of the nodes is in the new tree, using
  ## the standard names (parent-offspring)
  ## First, for phy2
  po.phy <- cbind(from=phy.names[phy$edge[,1]],
                  to=phy.names[phy$edge[,2]])
  po.phy2 <- cbind(from=phy2.names[phy2$edge[,1]],
                   to=phy2.names[phy2$edge[,2]])

  ## Then find out where the parent/offspring relationship changed:
  ## i <- match(po.phy2[,2], po.phy[,2])
  j <- which(po.phy[match(po.phy2[,2], po.phy[,2]),1] != po.phy2[,1])

  for ( idx in j ) {
    to <- po.phy2[idx,2]
    from <- po.phy2[idx,1]
    ans <- to
    offset <- 0
    repeat {
      to <- po.phy[po.phy[,2] == to,1]
      ans <- c(to, ans)
      if ( is.na(to) )
        stop("Horrible error")
      if ( to == from )
        break
    }
    
    if ( any(ans[-1] %in% hist$name2) ) {
      k <- hist$name2 %in% ans[-1]
      offset <- cumsum(phy$edge.length[match(ans[-1], po.phy[,2])])
      offset <- c(0, offset[-length(offset)])
      hist$x0[k] <- hist$x0[k] - offset[match(hist$name2[k], ans[-1])]
      hist$tc[k] <- hist$t[k] - hist$x0[k]
      hist$name2[k] <- ans[length(ans)]
    }
  }

  ## Prune out the extinct species and nodes that lead to them.  Note
  ## that the root must be excluded as history objects that lead to
  ## the new root (if it has changed) should not be allowed.
  phy2.names.noroot <- phy2.names[phy2.names != phy2$node.label[1]]
  hist <- hist[hist$name2 %in% phy2.names.noroot,]

  ## Remake idx2 to point at the new tree.
  hist$idx2 <- match(hist$name2, phy2.names)

  hist[order(hist$idx2, hist$t),]
}
