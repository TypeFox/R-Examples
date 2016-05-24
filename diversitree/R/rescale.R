## I'm going to take code here from arbutus, rather than from geiger,
## because we're tracking it more carefully and it's already tested.
## Once arbutus is released, perhaps we can either depend on arbutus
## or we can pull some of these utilities out.

## The difference here is that SE rescaling is handled separately.  We
## avoid scaling by s^2 (that's for the BM calculations to deal
## with).  Parameter vectors are handled differently (and assumed
## correct).
branching.times.edge <- function(phy) {
  ht <- unname(branching.heights(phy)) # same as node.depth.edgelength()
  list(start=ht[phy$edge[,1]],
       end=ht[phy$edge[,2]],
       max=max(ht))
}

## The direct translation from geiger via arbutus is a rescaling of:
##     tau.max   <- 2 * alpha * t.max
##     tau.start <- 2 * alpha * t.start
##     tau.end   <- 2 * alpha * t.end
##     exp(-tau.max) / (2 * alpha) * (
##       exp(tau.start) * expm1(-tau.start) -
##       exp(tau.end)   * expm1(-tau.end))
## but there is a useful simplification here:
##     exp(x)(1 - exp(-x)) ---> exp(x) - 1
## so we have
##     exp(-tau.max) * (expm1(tau.end) - expm1(tau.start)) / (2 * alpha)
## or
##     exp(-tau.max) * (exp(tau.end) - exp(tau.start)) / (2 * alpha)
## The first of which seems to be more numerically stable as alpha
## tends towards zero.
##
## So, given depth in the tree (t.max - t.start) and len, this
## is equivalent:
##     len <- t.end - t.start
##     exp(-2 * alpha * (t.max - t.start)) * expm1(2 * alpha * len)/(2*alpha)
## where t.max - t.start is of course the depth of the bottom of the
## branch.  That will come in useful elsewhere.  See model-ou-pruning.R
make.rescale.phylo.ou <- function(phy) {
  t <- branching.times.edge(phy)
  t.start <- t$start
  t.end   <- t$end
  t.max   <- t$max

  function(alpha) {
    if (alpha > 0) {
      phy$edge.length <-
        exp(-2 * alpha * t.max) * (
          # exp  (2 * alpha * t.end) - exp  (2 * alpha * t.start)) /
          expm1(2 * alpha * t.end) - expm1(2 * alpha * t.start)) /
          (2 * alpha)
    }
    phy
  }
}
make.rescale.phylo.eb <- function(phy, pars) {
  t <- branching.times.edge(phy)
  t.start <- t$start
  t.end   <- t$end
  t.max   <- t$max

  function(a) {
    if (a != 0) {
      phy$edge.length <-
        (exp(a * t.end) - exp(a * t.start)) / a
    }
    phy
  }
}

make.rescale.phylo.lambda <- function(phy) {
  n.tip <- length(phy$tip.label)
  interns <- which(phy$edge[, 2] >  n.tip)
  externs <- which(phy$edge[, 2] <= n.tip)

  t <- branching.times.edge(phy)
  height <- t$end[externs]

  function(lambda) {
    phy$edge.length[interns] <-
      phy$edge.length[interns] * lambda
    phy$edge.length[externs] <-
      phy$edge.length[externs] * lambda + (1 - lambda) * height
    phy
  }
}

rescale.phylo.se <- function(phy, se) {
  tips <- phy$edge[,2] <= length(phy$tip.label)
  phy$edge.length[tips] <- phy$edge.length[tips] + se
  phy
}

## A little wrapping helper function.
make.rescale.phylo <- function(phy, model) {
  switch(model,
         ou=make.rescale.phylo.ou,
         eb=make.rescale.phylo.eb,
         lambda=make.rescale.phylo.lambda,
         stop("Unknown model ", dQuote(model)))(phy)
}
