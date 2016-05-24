make.tree.bd <- function(pars, max.taxa=Inf, max.t=Inf) {
  extinct <- FALSE
  split   <- FALSE
  parent <- 0

  lambda <- pars[1]
  mu <- pars[2]
  r <- lambda + mu
  len <- 0
  t <- 0
  n.taxa <- 1
  lineages <- 1

  pr.speciation <- lambda/(lambda + mu)

  while ( n.taxa <= max.taxa && n.taxa > 0) {
    ## When does an event happen?
    r.n <- r * n.taxa
    r.tot <- sum(r.n)
    dt <- rexp(1, r.tot)
    t <- t + dt

    if ( t > max.t ) {
      dt <- dt - (t - max.t)
      len[lineages] <- len[lineages] + dt
      t <- max.t
      break
    }

    len[lineages] <- len[lineages] + dt

    ## Pick a lineage for the event to happen to:
    lineage.i <- sample(n.taxa, 1)
    lineage <- lineages[lineage.i]

    if ( runif(1) < pr.speciation ) {
      ## Speciating:
      if ( n.taxa == max.taxa )
        ## Don't add this one
        break
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE
      parent[new.i] <- lineage
      len[new.i] <- 0

      n.taxa <- n.taxa + 1

      ## lineages <- which(!split & !extinct)
      lineages <- c(lineages[-lineage.i], new.i)
    } else {
      ## Extinct
      extinct[lineage] <- TRUE
      ## lineages <- which(!split & !extinct)
      lineages <- lineages[-lineage.i]
      n.taxa <- n.taxa - 1
    }

  }

  info <- data.frame(idx=seq_along(extinct), len=len, parent=parent,
                     extinct=extinct, split=split)
  attr(info, "t") <- t
  info
}

me.to.ape.bd <- function(x) {
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

  phy <- reorder(structure(list(edge=cbind(x$parent2, x$idx2),
                                Nnode=Nnode,
                                tip.label=tip.label,
                                node.label=node.label,
                                edge.length=x$len,
                                orig=x),
                           class="phylo"))
  phy
}


tree.bd <- function(pars, max.taxa=Inf, max.t=Inf,
                    include.extinct=FALSE) {
  info <- make.tree.bd(pars, max.taxa, max.t)
  phy <- me.to.ape.bd(info[-1,])
  if ( include.extinct || is.null(phy) )
    phy
  else
    prune(phy)
}

tree.yule <- function(pars, max.taxa=Inf, max.t=Inf,
                      include.extinct=FALSE) {
  if ( length(pars) != 1 )
    stop("pars must be of length 1")
  info <- make.tree.bd(c(pars, 0), max.taxa, max.t)
  phy <- me.to.ape.bd(info[-1,])
  if ( include.extinct || is.null(phy) )
    phy
  else
    prune(phy)
}
