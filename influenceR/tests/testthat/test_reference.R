library(testthat)
load("flo_results.RData")

context("ref-tests")

ens_test <- function(g) {
  ens <- vector("numeric", length(igraph::V(g)))
  for (i in igraph::V(g)) {
    s <- 0
    ni <- igraph::neighbors(g, i)
    di <- igraph::degree(g, i)
    for (j in ni) {
      t <- 0
      Q <- igraph::intersection(ni, igraph::neighbors(g, j))
      for (q in Q) {
        t <- t + 1/di
      }
      s <- s + 1 - t
    }
    ens[i] <- s
  }
  ens
}

bridge_test <- function(g) {
  n <- length(igraph::V(g))
  
  cohesion <- function(g) {
    D <- igraph::distances(g)
    Dr <- 1/D
    Dr[Dr==Inf] <- 0
    s <- sum(Dr)
    s/(n*(n-1))
  }
  
  C <- cohesion(g)
  
  bridge <- NULL
  
  for (v in igraph::V(g)) {
    edges <- igraph::incident(g, v)
    s <- 0
    for (e in edges) {
      eset <- igraph::E(g)[igraph::E(g) != e]
      g_ <- igraph::subgraph.edges(g, eset)
      C_ <- cohesion(g_)
      s <- s + (C - C_)
    }
    bridge <- c(bridge, s)
  }
  
  bridge / (2*igraph::degree(g))
}


test_that("ens matches reference function", {
  expect_equal(influenceR::ens(flo_graph), ens_test(flo_graph))
})


test_that("influenceR::constraint matches igraph", {
  ig_constraint <- igraph::constraint(flo_graph)
  ig_constraint[is.nan(ig_constraint)] <- 0
  expect_less_than(sum(abs(influenceR::constraint(flo_graph) - ig_constraint)), 0.1)
})

test_that("influenceR::betweenness is 2x igraph version", {
  expect_less_than(sum(abs(influenceR::betweenness(flo_graph) - 2*igraph::betweenness(flo_graph))), 0.0001)
})

test_that("bridging matches reference function", {
  flo_bridge <- bridge_test(flo_graph)
  flo_bridge[is.nan(flo_bridge)] <- 0
  expect_less_than(sum(abs(influenceR::bridging(flo_graph) - flo_bridge)), 0.0001)
})