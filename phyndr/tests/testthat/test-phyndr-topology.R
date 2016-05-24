context("phyndr_topology")

test_that("phyndr_topology", {
  phy <- read.tree("Conifer-timetree.tre")
  topology <- read.tree("pinales_topotree.tre")

  set.seed(1)
  keep <- sample(phy$tip.label, 100)
  phy <- drop_tip(phy, setdiff(phy$tip.label, keep))

  set.seed(1)
  data_species <- sample(union(phy$tip.label, topology$tip.label), 200)

  res <- phyndr_topology(phy, data_species, topology)
  expect_that(length(res$tip.label), equals(44))

  expect_that(res$clades, is_a("list"))
  expect_that(phyndr_n_distinct(res), equals(208))

  expect_that(phyndr_sample(phy), throws_error("Expected a phyndr"))
  set.seed(1)
  phy2 <- phyndr_sample(res)
  expect_that(phy2, not(is_a("phyndr")))
  expect_that(phy2$clades, is_null())
  expect_that(phy2$tip.label, not(equals(res$tip.label)))

  set.seed(1)
  tmp <- phyndr_sample_n(res, 5)
  expect_that(tmp[[1]], equals(phy2))
  expect_that(length(tmp), equals(5))

  expect_that(phyndr_sample_n(res, 0), is_identical_to(list()))

  set.seed(1)
  expect_that(phyndr_sample_n(res, 1), is_identical_to(list(phy2)))

  expect_that(phyndr_combn(phy), throws_error("Expected a phyndr"))
  phyl <- phyndr_combn(res)
  expect_that(phyl, is_a("list"))
  expect_that(length(phyl), equals(phyndr_n_distinct(res)))
  expect_that(phyl[[1]], is_a("phylo"))
  expect_that(phyl[[1]], not(is_a("phyndr")))
  expect_that(phyl[[1]]$tip.label, not(equals(res$tip.label)))

  ok <- vlapply(phyl, function(x) identical(x$tip.label, res$tip.label))
  expect_that(sum(ok), equals(1L))

  expect_that(phyndr_combn(res, 10),
              throws_error("Too many trees would be generated"))
})
