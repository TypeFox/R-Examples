source("helper-diversitree.R")
context("Tree rescaling")

## Here is some rescaling functions I trust.
if (suppressWarnings(require("arbutus", quietly=TRUE))) {
  set.seed(1)
  phy <- tree.yule(1, max.taxa=30)

  ## This is a regression test.
  test_that("Time extraction works with different tree orders", {
    t <- diversitree:::branching.times.edge(phy)
    expect_that(t$end - t$start, equals(phy$edge.length))
    tree <- reorder(phy, "pruningwise")
    t <- diversitree:::branching.times.edge(tree)
    expect_that(t$end - t$start, equals(tree$edge.length))
    tree <- reorder(phy, "cladewise")
    t <- diversitree:::branching.times.edge(tree)
    expect_that(t$end - t$start, equals(tree$edge.length))
  })

  test_that("OU rescaling agrees with arbutus", {
    r.ou <- diversitree:::make.rescale.phylo.ou(phy)

    ## Nothing happens with alpha 0
    expect_that(r.ou(0), is_identical_to(phy))
    p <- function(alpha) list(alpha=alpha, sigsq=1, SE=0)

    expect_that(r.ou(.1),
                equals(arbutus:::model.phylo.ou(phy, p(.1))))
    ## Calculations can get unstable with very small alpha - can cause
    ## complete failure.  Unfortunately all.equal.phylo does not honour
    ## additional arguments, so unclassing is required here.
    expect_that(unclass(r.ou(1e-10)),
                equals(unclass(arbutus:::model.phylo.ou(phy, p(1e-10))),
                       tolerance=1e-6))
  })

  test_that("EB rescaling agrees with arbutus", {
    r.eb <- diversitree:::make.rescale.phylo.eb(phy)

    ## Nothing happens with a 0
    expect_that(r.eb(0), is_identical_to(phy))

    p <- function(a) list(a=a, sigsq=1, SE=0)
    expect_that(r.eb(.1),
                equals(arbutus:::model.phylo.eb(phy, p(.1))))
    expect_that(r.eb(-.1),
                equals(arbutus:::model.phylo.eb(phy, p(-.1))))
  })

  test_that("Lambda rescaling agrees with arbutus", {
    r.lambda <- diversitree:::make.rescale.phylo.lambda(phy)

    ## Nothing happens with lambda 1
    expect_that(r.lambda(1), is_identical_to(phy))

    p <- function(lambda) list(lambda=lambda, sigsq=1, SE=0)
    expect_that(r.lambda(.3),
                equals(arbutus:::model.phylo.lambda(phy, p(.3))))
  })

  test_that("SE rescaling agrees with arbutus", {
    expect_that(diversitree:::rescale.phylo.se(phy, 0),
                is_identical_to(phy))
    expect_that(diversitree:::rescale.phylo.se(phy, .1),
                equals(arbutus:::model.phylo.se(phy, list(SE=.1))))
  })

  ## Now, try on an non-ultrametric tree:
  set.seed(1)
  phy <- rtree(30)

  test_that("OU rescaling agrees with arbutus", {
    r.ou <- diversitree:::make.rescale.phylo.ou(phy)

    ## Nothing happens with alpha 0
    expect_that(r.ou(0), is_identical_to(phy))
    p <- function(alpha) list(alpha=alpha, sigsq=1, SE=0)

    expect_that(r.ou(.1),
                equals(arbutus:::model.phylo.ou(phy, p(.1))))
  })

  test_that("EB rescaling agrees with arbutus", {
    r.eb <- diversitree:::make.rescale.phylo.eb(phy)

    ## Nothing happens with a 0
    expect_that(r.eb(0), is_identical_to(phy))

    p <- function(a) list(a=a, sigsq=1, SE=0)
    expect_that(r.eb(.1),
                equals(arbutus:::model.phylo.eb(phy, p(.1))))
    expect_that(r.eb(-.1),
                equals(arbutus:::model.phylo.eb(phy, p(-.1))))
  })

  test_that("Lambda rescaling agrees with arbutus", {
    r.lambda <- diversitree:::make.rescale.phylo.lambda(phy)

    ## Nothing happens with lambda 1
    expect_that(r.lambda(1), is_identical_to(phy))

    p <- function(lambda) list(lambda=lambda, sigsq=1, SE=0)
    expect_that(r.lambda(.3),
                equals(arbutus:::model.phylo.lambda(phy, p(.3))))
  })

  test_that("SE rescaling agrees with arbutus", {
    expect_that(diversitree:::rescale.phylo.se(phy, 0),
                is_identical_to(phy))
    expect_that(diversitree:::rescale.phylo.se(phy, .1),
                equals(arbutus:::model.phylo.se(phy, list(SE=.1))))
  })
}
