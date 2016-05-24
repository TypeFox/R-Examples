context("SumStat MCMF")

ss <- matrix(c(1, 0, 0, 1,
               0, 1, 0, 0,
               1, 0, 1, 0,
               1, 0, 0, 0), 4, 4, byrow = TRUE)
seg_sites <- list(create_segsites(ss, c(0.1, 0.2, 0.5, 0.7)))


test_that("mcmf is correctly calculation for normal loci", {
  expect_equal(calc_mcmf(seg_sites, 1:4, FALSE), .5)
  expect_equal(calc_mcmf(seg_sites, c(1, 3, 4), FALSE), .5)
  expect_equal(calc_mcmf(seg_sites, 2:4, FALSE), 2 / 3)
  expect_equal(calc_mcmf(seg_sites, 3:4, FALSE), 1)
})


test_that("mcmf is correctly calculation for locus trios", {
  seg_sites <- list(create_segsites(cbind(ss, ss, ss),
                                    rep(c(0.1, 0.2, 0.5, 0.7), 3),
                                    rep(c(-1, 0, 1), each = 4)))

  expect_equal(calc_mcmf(seg_sites, 1:4), c(4 / 12))
  expect_equal(calc_mcmf(seg_sites, 2:4), c(4 / 9))
  expect_equal(calc_mcmf(seg_sites, 3:4), c(2 / 3))

  ss <- matrix(c(0, 0, 0, 1,
                 0, 0, 1, 0,
                 0, 0, 1, 0,
                 0, 0, 1, 0), 4, 4, byrow = TRUE)
  seg_sites[[2]] <- create_segsites(cbind(ss, ss, ss),
                                    rep(c(0.1, 0.2, 0.5, 0.7), 3),
                                    rep(c(-1, 0, 1), each = 4))
  expect_equal(calc_mcmf(seg_sites, 1:4), c(c(4 / 12), c(4 / 6)))
  expect_equal(calc_mcmf(seg_sites, 2:4), c(c(4 / 9), NA))
  expect_error(calc_mcmf(seg_sites, 1:5))

  seg_sites <- list(create_segsites(matrix(0, 4, 0), numeric()))
  expect_true(is.na(calc_mcmf(seg_sites, 1:4)))
})


test_that("mcmf is correctly calculation for diplod loci", {
  expect_equal(calc_mcmf(seg_sites, 1:2, FALSE, ploidy = 2), .75)
  expect_equal(calc_mcmf(seg_sites, 1, FALSE, ploidy = 2), 1)
  expect_equal(calc_mcmf(seg_sites, 2, FALSE, ploidy = 2), 1)

  seg_sites[[1]] <- create_segsites(rbind(ss, ss), c(0.1, 0.2, 0.5, 0.7))
  expect_equal(calc_mcmf(seg_sites, 1:3, FALSE, ploidy = 2), .75)
  expect_equal(calc_mcmf(seg_sites, 1:4, FALSE, ploidy = 2), .75)
  expect_error(calc_mcmf(seg_sites, 1:5, FALSE, ploidy = 2))
  expect_error(calc_mcmf(seg_sites, 1:4, FALSE, ploidy = 3))
})


test_that("mcmf is correctly calculation for trioplod loci", {
  seg_sites[[1]] <- create_segsites(rbind(ss, ss), c(0.1, 0.2, 0.5, 0.7))
  expect_equal(calc_mcmf(seg_sites, 1:2, FALSE, ploidy = 3), .75)
})


test_that("initialzation of statistic works", {
  stat <- sumstat_mcmf(population = 1)
  expect_equal(stat$calculate(seg_sites, NULL, NULL, coal_model(4)), .5)

  seg_sites <- list(create_segsites(cbind(ss, ss, ss),
                                    rep(c(0.1, 0.2, 0.5, 0.7), 3),
                                    rep(c(-1, 0, 1), each = 4)))
  expect_equal(stat$calculate(seg_sites, NULL, NULL,
                              coal_model(4) + locus_trio()), 1 / 3)
})


test_that("mcmf statistics is correct for diploid models", {
  stat <- sumstat_mcmf(population = 1)
  model <- coal_model(2, ploidy = 2)
  expect_equal(stat$calculate(seg_sites, NULL, NULL, model), .5)
  expect_equal(stat$calculate(seg_sites, NULL, NULL,
                              model + feat_unphased(2)), .75)
})
