##' Number of possible distinct trees that can be generated from a
##' phyndr set.
##' @title Number of distinct trees
##' @param phy A phyndr phylogeny
##' @export
phyndr_n_distinct <- function(phy) {
  prod(viapply(phy$clades, length))
}

##' Sample distinct set of tip labels from a phyndr tree.
##' \code{phyndr_sample} returns a tree (phylo object) while
##' \code{phyndr_sample_n} generates a list of trees.
##' @title Sample one distinct set of tip labels from a phyndr tree
##' @param phy A phyndr tree
##' @export
phyndr_sample <- function(phy) {
  if (!inherits(phy, "phyndr")) {
    stop("Expected a phyndr tree")
  }
  repl <- vcapply(phy$clades, sample, 1L)
  phy$tip.label[match(names(repl), phy$tip.label)] <- unname(repl)
  phy$clades <- NULL
  class(phy) <- setdiff(class(phy), "phyndr")
  phy
}

##' @export
##' @rdname phyndr_sample
##' @param n Number of trees to generate
phyndr_sample_n <- function(phy, n) {
  replicate(n, phyndr_sample(phy), FALSE)
}

##' Generate all distinct sets of tip labels from a phyndr tree
##' @title Generate all distinct sets of tip labels from a phyndr tree
##' @param phy A phyndr tree
##' @param max_n Maximum number of trees to generate; if the number of
##' distinct trees is larger than this, we'll throw an error instead.
##' @export
phyndr_combn <- function(phy, max_n=1000) {
  if (!inherits(phy, "phyndr")) {
    stop("Expected a phyndr tree")
  }
  if (phyndr_n_distinct(phy) > max_n) {
    stop("Too many trees would be generated")
  }

  ## Do all the processing once:
  idx <- match(names(phy$clades), phy$tip.label)
  phy2 <- phy
  phy2$clades <- NULL
  class(phy2) <- setdiff(class(phy2), "phyndr")

  f <- function(x) {
    phy2$tip.label[idx] <- unname(x)
    phy2
  }

  repl <- as.matrix(do.call("expand.grid", phy$clades))
  apply(repl, 1, f)
}
