context("SumStat JSFS")

test_that("calculation of the JSFS is correct", {
  seg_sites_3_pop <- list(create_segsites(matrix(c(1, 0, 1, 0,
                                                   1, 1, 0, 1,
                                                   0, 0, 0, 1,
                                                   1, 0, 0, 1,
                                                   1, 1, 1, 1,
                                                   1, 0, 1, 0,
                                                   1, 0, 1 ,0 ),
                                                 nrow=7, ncol=4, byrow = TRUE),
                                          c(0.1, 0.2, 0.5, 0.7)))
  jsfs_3_pop <- calc_jsfs(seg_sites_3_pop, list(1:2, 3:4, 5:7))
  expect_is(jsfs_3_pop, "array")
  expect_equivalent(jsfs_3_pop, array(c(0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 1, 0, 0, 0, 0, 0, 1, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 1, 0, 0, 0, 1, 0, 0, 0),
                                      dim=c(2, 2, 3) + 1))


  seg_sites_3_pop <- list(create_segsites(matrix(c(1, 0, 1, 0,
                                                   1, 1, 0, 1,
                                                   0, 0, 1, 1,
                                                   0, 0, 0, 1,
                                                   1, 1, 1, 1,
                                                   1, 0, 1, 0,
                                                   1, 0, 1 ,0 ),
                                                 nrow=7, ncol=4, byrow = TRUE),
                                          c(0.1, 0.2, 0.5, 0.7)))
  jsfs_3_pop <- calc_jsfs(seg_sites_3_pop, list(1, 2:6, 7))
  expect_equivalent(jsfs_3_pop, array(c(0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0),
                                      dim=c(1, 5, 1) + 1))

  seg_sites <- list(create_segsites(matrix(c(1, 0, 1, 0,
                                             1, 1, 0, 1,
                                             0, 0, 0, 1,
                                             1, 0, 0, 1), 4, 4, byrow = TRUE),
                                    c(0.1, 0.2, 0.5, 0.7)))

  jsfs <- calc_jsfs(seg_sites, list(1:2, 3:4))
  expect_is(jsfs, "matrix")
  expect_equivalent(jsfs, matrix(c(0, 0, 0,
                                   2, 0, 1,
                                   0, 1, 0), 3, 3, byrow = TRUE))

  jsfs <- calc_jsfs(seg_sites, list(1, 2:4))
  expect_equivalent(jsfs, matrix(c(0, 1, 0, 1,
                                   1, 0, 1, 0), 2, 4, byrow = TRUE))

  jsfs2 <- calc_jsfs(seg_sites, list(2:4, 1))
  expect_equivalent(jsfs2, t(jsfs))

  jsfs <- calc_jsfs(seg_sites, list(c(1,3), c(2,4)))
  expect_equivalent(jsfs, matrix(c(0, 1, 0,
                                   1, 0, 2,
                                   0, 0, 0), 3, 3, byrow = TRUE)) #nolint

  expect_error(calc_jsfs(seg_sites, c(1,3), c(2,5)))
  expect_error(calc_jsfs(seg_sites, c(1,7), c(2,3)))
  calc_jsfs(seg_sites, list(numeric(), 1:4))


  seg_sites <- list(create_segsites(matrix(c(1, 1, 1, 1,
                                             0, 0, 1, 1,
                                             1, 1, 0, 0), 4, 3),
                                    c(0.1, 0.5, 0.7)))

  jsfs <- calc_jsfs(seg_sites, list(1:2, 3:4))
  expect_true(is.matrix(jsfs))
  expect_equal(dim(jsfs), c(3, 3))
  expect_equal(sum(jsfs), 2)
  expect_equal(jsfs[1, 3], 1)
  expect_equal(jsfs[3, 1], 1)


  seg_sites[[2]] <- create_segsites(matrix(c(0, 1, 1, 1,
                                             0, 1, 1, 1,
                                             0, 1, 1, 1),  4, 3),
                                    c(0.1, 0.5, 0.7))
  jsfs <- calc_jsfs(seg_sites, list(1:2, 3:4))
  expect_true(is.matrix(jsfs))
  expect_equal(dim(jsfs), c(3, 3))
  expect_equal(sum(jsfs), 5)
  expect_equal(jsfs[2, 3], 3)
  expect_equal(jsfs[1, 3], 1)
  expect_equal(jsfs[3, 1], 1)


  seg_sites[[3]] <- create_empty_segsites(4)
  jsfs <- calc_jsfs(seg_sites, list(1:2, 3:4))
  expect_true(is.matrix(jsfs))
  expect_equal(dim(jsfs), c(3, 3))
  expect_equal(sum(jsfs), 5)
  expect_equal(jsfs[2, 3], 3)
  expect_equal(jsfs[1, 3], 1)
  expect_equal(jsfs[3, 1], 1)
})


test_that("calc_jsfs works with trios", {
  ss <- matrix(c(1, 0, 0, 0,
                 1, 1, 0, 1,
                 1, 0, 0, 1,
                 1, 0, 0, 1), 4, 4, byrow = TRUE)

  seg_sites <- list(create_segsites(cbind(ss, ss, ss),
                                    rep(c(0.1, 0.2, 0.5, 0.7), 3),
                                    rep(c(-1, 0, 1), each = 4)))

  jsfs <- calc_jsfs(seg_sites, list(1:2, 3:4))
  expect_equivalent(jsfs, matrix(c(0, 0, 0,
                                   1, 0, 1,
                                   0, 0, 0), 3, 3, byrow = TRUE))
})


test_that("JSFS sumstat works", {
  stat <- sumstat_jsfs("jsfs_test", c(1, 2))
  model <- coal_model(c(2, 2), 1)

  seg_sites <- list(create_segsites(matrix(c(1, 0, 0, 0,
                                             1, 1, 0, 1,
                                             1, 0, 0, 1,
                                             1, 0, 0, 1), 4, 4, byrow = TRUE),
                                    c(0.1, 0.2, 0.5, 0.7)))

  expect_equal(stat$get_name(), "jsfs_test")
  expect_equivalent(stat$calculate(seg_sites, NULL, NULL, model),
                    matrix(c(0, 0, 0,
                             1, 0, 1,
                             0, 0, 0), 3, 3, byrow = TRUE))
})


test_that("JSFS is caluculated with an outgroup present", {
  stat <- sumstat_jsfs("jsfs", c(1, 2))
  model <- coal_model(c(2, 1, 1), 1) + feat_outgroup(3)
  seg_sites <- list(create_segsites(matrix(c(1, 0, 0, 0,
                                             1, 1, 0, 1,
                                             1, 0, 0, 1), 3, 4, byrow = TRUE),
                                    1:4 / 4))

  expect_equivalent(stat$calculate(seg_sites, NULL, NULL, model),
                    matrix(c(0, 0,
                             1, 1,
                             0, 0), 3, 2, byrow = TRUE))

  model <- coal_model(c(2, 1, 1), 1) + feat_outgroup(2)
  expect_error(stat$calculate(seg_sites, NULL, NULL, model))

  stat <- sumstat_jsfs("jsfs", c(1, 3))
  expect_equivalent(stat$calculate(seg_sites, NULL, NULL, model),
                    matrix(c(0, 0,
                             1, 1,
                             0, 0), 3, 2, byrow = TRUE))
})


test_that("Per locus calculation of JSFS works", {
  seg_sites_list <- list(create_test_segsites(),
                         create_test_segsites(),
                         create_test_segsites())

  model <- coal_model(c(2, 1), 3)
  jsfs <- sumstat_jsfs("jsfs", c(1, 2), per_locus = TRUE)
  stat <- jsfs$calculate(seg_sites_list, NULL, NULL, model)
  expect_that(stat, is_a("list"))
  expect_equal(stat[[1]], stat[[2]])
  expect_equal(stat[[1]], stat[[3]])

  jsfs_full <- sumstat_jsfs("jsfs", c(1, 2), per_locus = FALSE)
  expect_equal(3 * stat[[1]],
               jsfs_full$calculate(seg_sites_list, NULL, NULL, model))
})
