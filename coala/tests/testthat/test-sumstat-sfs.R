context("SumStat SFS")

test_that("calculation of the SFS is correct", {
  seg.sites <- list(create_segsites(matrix(c(1, 0, 0, 0,
                                             1, 1, 0, 1,
                                             0, 0, 1, 1,
                                             1, 1, 0, 1), 4, 4, byrow = TRUE),
                                    c(0.1, 0.2, 0.5, 0.7)))

  model <- coal_model(4, 1)
  stat_sfs_all <- sumstat_sfs(population = "all")
  stat_sfs_1 <- sumstat_sfs(population = 1)
  stat_sfs_2 <- sumstat_sfs(population = 2)

  expect_equal(stat_sfs_1$calculate(seg.sites, NULL, NULL, model), c(1, 1, 2))
  expect_equal(stat_sfs_all$calculate(seg.sites, NULL, NULL, model), c(1, 1, 2))

  model <- coal_model(c(2, 2), 1)
  expect_equal(stat_sfs_all$calculate(seg.sites, NULL, NULL, model), c(1, 1, 2))
  expect_equal(stat_sfs_1$calculate(seg.sites, NULL, NULL, model), 2)
  expect_equal(stat_sfs_2$calculate(seg.sites, NULL, NULL, model), 3)

  seg.sites[[2]] <- create_segsites(matrix(c(1, 1, 1,
                                             1, 1, 1,
                                             1, 1, 1,
                                             1, 1, 1),  4, 3),
                                    c(0.1, 0.5, 0.7))
  expect_equal(stat_sfs_all$calculate(seg.sites, NULL, NULL, model), c(1, 1, 2))
  expect_equal(stat_sfs_1$calculate(seg.sites, NULL, NULL, model), c(2))
  expect_equal(stat_sfs_2$calculate(seg.sites, NULL, NULL, model), c(3))

  seg.sites[[3]] <- create_empty_segsites(4)
  expect_equal(stat_sfs_all$calculate(seg.sites, NULL, NULL, model), c(1, 1, 2))
  expect_equal(stat_sfs_1$calculate(seg.sites, NULL, NULL, model), 2)
  expect_equal(stat_sfs_2$calculate(seg.sites, NULL, NULL, model), 3)
})


test_that("calculation of sfs works with trios", {
  ss <- matrix(c(1, 0, 0, 0,
                 1, 1, 0, 1,
                 1, 0, 0, 1,
                 1, 0, 0, 1), 4, 4, byrow = TRUE)

  seg.sites <- list(create_segsites(cbind(ss, ss, ss),
                                    rep(c(0.1, 0.2, 0.5, 0.7), 3),
                                    rep(c(-1, 0, 1), each = 4)))

  model <- coal_model(4, 1)
  stat_sfs_all <- sumstat_sfs(population = "all")
  expect_equal(stat_sfs_all$calculate(seg.sites, NULL, NULL, model), c(1, 0, 1))
})


test_that("SFS is calculate with an outgroup present", {
  model <- coal_model(c(2, 1, 1), 1) + feat_outgroup(3)
  stat <- sumstat_sfs("sfs", "all")

  seg_sites <- list(create_segsites(matrix(c(1, 0, 0, 0,
                                             1, 1, 0, 1,
                                             1, 0, 0, 1), 3, 4, byrow = TRUE),
                                    1:4 / 4))
  expect_equal(stat$calculate(seg_sites, NULL, NULL, model), c(1, 1))

  stat <- sumstat_sfs("jsfs", 3)
  expect_error(stat$calculate(seg_sites, NULL, NULL, model))

  model <- coal_model(c(2, 1, 1), 1) + feat_outgroup(2)
  expect_equal(sumstat_sfs("jsfs", 1)$calculate(seg_sites, NULL, NULL, model),
               2)
  expect_equal(sumstat_sfs("jsfs", 3)$calculate(seg_sites, NULL, NULL, model),
               numeric(0))
})
