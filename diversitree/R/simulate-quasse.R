## lambda, mu are speciation and extinction functions and return the
## speciation and extinction rates, given a vector of character states
## 'x'.
## char is a function describing character evolution, randomly chosing
## a new character given a vector of current characters 'x' and a
## period of time 'dt' over which evolution occurs.
tree.quasse <- function(pars, max.taxa=Inf, max.t=Inf,
                        include.extinct=FALSE, x0=NA,
                        single.lineage=TRUE, verbose=FALSE) {
  if ( is.na(x0) )
    stop("x0 must be specified")
  else if ( length(x0) != 1 )
    stop("x0 must be of length 1")
  stopifnot(is.list(pars), all(sapply(pars, is.function)))

  info <- make.tree.quasse(pars, max.taxa, max.t, x0, single.lineage,
                           verbose)
  if ( single.lineage )
    info <- info[-1,]
  phy <- me.to.ape.quasse(info)
  if ( include.extinct || is.null(phy) )
    phy
  else
    prune(phy)
}

make.tree.quasse <- function(pars, max.taxa=Inf, max.t=Inf, x0,
                             single.lineage=TRUE,
                             verbose=FALSE, k=500, ...) {
  lambda <- pars[[1]]
  mu     <- pars[[2]]
  char   <- pars[[3]]
  
  ## Always assuming *two* lineages for now, both in state x0.  I am
  ## adding a very small amount of time here so that the root does not
  ## have a real zero branch length (this happens 1/k of the time, and
  ## breaks quassecache).
  if ( single.lineage ) {
    info <- data.frame(idx=1, len=1e-8, parent=0, state=x0,
                       extinct=FALSE, split=FALSE)
  } else {
    info <- data.frame(idx=1:2, len=1e-8, parent=0, state=x0,
                       extinct=FALSE, split=FALSE)
  }

  lineages <- which(!info$extinct & !info$split)
  n.taxa <- length(lineages)
  t <- 0
  t.left <- max.t

  while ( n.taxa <= max.taxa && n.taxa > 0 && t.left > 0 ) {
    x <- run.until.change(lineages, info, k, lambda, mu, char, t.left)
    lineages <- x[[1]]
    info <- x[[2]]
    n.taxa <- length(lineages)
    t <- t + x[[4]]
    t.left <- t.left - x[[4]]
    if ( verbose )
      cat(sprintf("%s: %d [%2.3f]\n",
                  c("-", " ", "+")[sign(x[[3]])+2], n.taxa, t))
  }

  if ( n.taxa > max.taxa ) {
    ## Drop final speciation event.
    drop <- info[nrow(info)-1,]
    info$split[drop$parent] <- FALSE
    info$state[drop$parent] <- drop$state
    info$len[drop$parent] <- info$len[drop$parent] + drop$len
    info <- info[seq_len(nrow(info)-2),]
  }

  attr(info, "t") <- t
  info
}

## This function extends lineages until either a speciation or
## extinction event happens.  The time interval is rescaled so that
## 'k' character change events are expected between each speciation
## and extinction event.
## TODO: I am not truncating exactly the time interval when we might
## run into the maximum time.  If 'k' is sufficiently high, this
## should not matter much.
run.until.change <- function(lineages, info, k, lambda, mu, char,
                             max.t) {
  i <- 1
  time <- 0
  n.extant <- length(lineages)
  p.change <- 1/k
  niter <- 1
  repeat {
    state <- info$state[lineages]
    lx <- lambda(state)
    mx <- mu(state)
    r <- sum(lx + mx)
    dt <- 1/(r*k)
    
    if ( runif(1) < p.change ) {
      if ( runif(1) < sum(lx)/r ) { # speciation
        i <- sample(n.extant, 1, prob=lx)
        info <- speciate(info, lineages[i])
        lineages <- c(lineages[-i], c(-1,0) + nrow(info))
      } else {
        i <- sample(n.extant, 1, prob=mx)
        info$extinct[lineages[i]] <- TRUE
        lineages <- lineages[-i]
      }
      info$len[lineages]   <- info$len[lineages] + dt
      info$state[lineages] <- char(info$state[lineages], dt)
      time <- time + dt
      break
    }

    info$len[lineages]   <- info$len[lineages] + dt
    info$state[lineages] <- char(info$state[lineages], dt)
    niter <- niter + 1
    time <- time + dt

    if ( time > max.t )
      break
  }

  list(lineages, info, length(lineages) - n.extant, time)
}

## This function rejigs the 'lineages' structure when speciation
## happens, creating new species.
speciate <- function(info, i) {
  j <- 1:2 + nrow(info)
  info[j,"idx"] <- j
  info[j,-1] <- list(len=0, parent=i,
                     state=info$state[i],
                     extinct=FALSE, split=FALSE)
  info$split[i] <- TRUE
  info
}

## Transform the lineages structure produced by sim.tree into an ape
## phylogeny.
me.to.ape.quasse <- function(info) {
  if ( nrow(info) == 0 )
    return(NULL)
  Nnode <- sum(!info$split) - 1
  n.tips <- sum(!info$split)

  info$idx2 <- NA
  info$idx2[!info$split] <- 1:n.tips
  info$idx2[ info$split] <- order(info$idx[info$split]) + n.tips + 1

  i <- match(info$parent, info$idx)
  info$parent2 <- info$idx2[i]
  info$parent2[is.na(info$parent2)] <- n.tips + 1

  tip.label <- ifelse(subset(info, !split)$extinct,
                      sprintf("ex%d", 1:n.tips),
                      sprintf("sp%d", 1:n.tips))
  node.label <- sprintf("nd%d", 1:Nnode)

  info$name <- NA
  info$name[!info$split] <- tip.label

  tip.state <- info$state[match(1:n.tips, info$idx2)]
  names(tip.state) <- tip.label
  phy <- reorder(structure(list(edge=cbind(info$parent2, info$idx2),
                               Nnode=Nnode,
                               tip.label=tip.label,
                               tip.state=tip.state,
                               node.label=node.label,
                               edge.length=info$len,
                               orig=info),
                          class="phylo"))

  phy$edge.state <- info$state[match(phy$edge[,2], info$idx2)]
  phy
}

## Pretty colour level splitting, similar to image()
## col.split <- function(x, cols=heat.colors(10))
##   cols[as.integer(cut(x, length(cols)))]

## ## Use:
## lambda.x <- make.sigmoid(0.07, 0.13, 0, 1)
## mu.x <- function(x) rep(0.02, length(x))
## char.evol <- make.brownian.with.drift(0, 0.01)
## cols <- heat.colors(10)

## ## Here's a simple tree
## set.seed(2)
## ans <- sim.tree(20, 500, lambda.x, mu.x, char.evol, x0=-2)

## ## This shows the entire evolutionary process.
## phy <- me.to.ape(ans)
## plot(phy, label.offset=1.1, cex=.75, show.node.label=TRUE,
##      edge.col=col.split(phy$edge.state))
## tip.col <- col.split(c(range(phy$edge.state),
##                        subset(phy$orig, !split)$state))[-(1:2)]
## tiplabels(pch=19, cex=0.85, adj=1, col=tip.col)

## ## I still haven't worked out how to get the internel nodes correct
## ## after pruning extinct taxa, though.
## phy2 <- prune(phy)
## plot(phy2, label.offset=1)
## tiplabels(pch=19, cex=1, adj=1,
##           col=col.split(subset(phy2$orig, !split)$state))

make.brownian.with.drift <- function(drift, diffusion)
  function(x, dt) x + rnorm(length(x), drift*dt, sqrt(dt*diffusion))
