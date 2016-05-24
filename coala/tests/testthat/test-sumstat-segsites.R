context("SumStat SegSites")

test_that("SegSites statistic works", {
  stat <- sumstat_seg_sites("segsites_test")

  seg_sites <- list(create_segsites(matrix(c(1, 0, 0, 0,
                                             1, 1, 0, 1,
                                             1, 0, 0, 1,
                                             1, 0, 0, 1), 4, 4, byrow = TRUE),
                                    c(0.1, 0.2, 0.5, 0.7)))

  expect_equal(stat$get_name(), "segsites_test")
  expect_equal(stat$calculate(seg_sites, NULL, NULL), seg_sites)
})


test_that("SegSites are convert for trios", {
  seg_sites <- list(create_segsites(matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 0,
                                             1, 1, 0, 1, 1, 1, 0, 1, 1,
                                             0, 0, 1, 1, 0, 0, 1, 1, 1,
                                             1, 0, 0, 1, 1, 0, 0, 1, 0),
                                           4, 9, byrow = TRUE), 1:9 / 10))
  expect_equal(conv_for_trios(seg_sites, coal_model(4) + locus_single(100)),
               seg_sites)

  model <- coal_model(4) + locus_trio(c(25, 30, 25), c(10, 10))
  seg_sites_trio <- seg_sites
  seg_sites_trio[[1]] <- create_segsites(get_snps(seg_sites_trio[[1]][ , c(1:2, 4:6, 8:9)]), #nolint
                                         c(.4, 0.8, c(5, 15, 25) / 30, .2, .6),
                                         c(-1, -1, 0, 0, 0, 1, 1))
  expect_equal(conv_for_trios(seg_sites, model), seg_sites_trio)


  seg_sites <- list(seg_sites[[1]], seg_sites[[1]])
  seg_sites_trio <- list(seg_sites_trio[[1]], seg_sites_trio[[1]])

  # Two loci groups
  model <- model + locus_trio(c(25, 30, 25), c(10, 10))
  expect_equal(conv_for_trios(seg_sites, model), seg_sites_trio)

  # Multiple loci in a group
  model <- coal_model(4) + locus_trio(c(25, 30, 25), c(10, 10), number = 2)
  expect_equal(conv_for_trios(seg_sites, model), seg_sites_trio)
})
