require(testthat)
require(pez)
context("Metapopulation and phylogeny simulation")

#Tests
# - no explicit test of edge2phylo because these test it
test_that("sim.meta.comm", expect_that(names(sim.meta.comm()), equals(c("comm","environment"))))
test_that("sim.meta.phy.comm", expect_that(class(sim.meta.phy.comm()), equals(c("comparative.comm", "comparative.data"))))
