context("SumStat iHS")

seg_sites <- create_segsites(matrix(c(1, 0, 0, 0, 1,
                                      1, 1, 0, 1, 0,
                                      1, 0, 1, 1, 1,
                                      0, 0, 0, 1, 0), 4, 5, byrow = TRUE),
                             c(0.1, 0.2, 0.5, 0.7, 0.9))
model <- coal_model(4, 1, 337)
pos <- get_snp_positions(list(seg_sites), model, relative = FALSE)[[1]]


test_that("snp masks are generated corretly", {
  skip_if_not_installed("rehh")

  stat_ihh <- sumstat_ihh(population = 1)
  expect_equal(stat_ihh$create_snp_mask(seg_sites), rep(TRUE, 5))

  stat_ihh <- sumstat_ihh(population = 1, max_snps = 2)
  for (i in 1:10) {
    snp_mask <- stat_ihh$create_snp_mask(seg_sites)
    expect_true(all(snp_mask %in% 1:5))
    expect_equal(length(snp_mask), 2)
    expect_equal(sort(snp_mask), snp_mask)
  }
})


test_that("generation of rehh data works", {
  skip_if_not_installed("rehh")
  stat_ihh <- sumstat_ihh(population = 1)
  rehh_data <- stat_ihh$create_rehh_data(seg_sites, 1:4, model)
  expect_equivalent(rehh_data@haplo, as.matrix(seg_sites) + 1)
  expect_equal(rehh_data@position, pos)
  expect_equal(rehh_data@snp.name, as.character(1:5))
  expect_equal(rehh_data@nhap, 4)
  expect_equal(rehh_data@nsnp, 5)

  rehh_data <- stat_ihh$create_rehh_data(create_empty_segsites(5), 1:5, model)
  expect_equal(rehh_data@haplo, matrix(0, 5, 0))

  rehh_data <- stat_ihh$create_rehh_data(seg_sites, numeric(), model)
  expect_equal(rehh_data@haplo, matrix(0, 0, 0))
})


test_that("SNPs not segregating in individuals are removed from rehh_data", {
  skip_if_not_installed("rehh")
  stat_ihh <- sumstat_ihh(population = 1)
  rehh_data <- stat_ihh$create_rehh_data(seg_sites, 1:2, model)
  expect_equivalent(rehh_data@haplo, as.matrix(seg_sites[1:2, c(2, 4, 5)]) + 1)
  expect_equal(rehh_data@position, pos[c(2, 4, 5)])
  expect_equal(rehh_data@snp.name, as.character(1:3))
  expect_equal(rehh_data@nhap, 2)
  expect_equal(rehh_data@nsnp, 3)
})


test_that("selection of snps works", {
  stat_ihh <- sumstat_ihh(population = 1, max_snps = 2)
  rehh_data <- stat_ihh$create_rehh_data(seg_sites, 1:4, model)
  expect_equal(dim(rehh_data@haplo), c(4, 2))
  expect_equal(rehh_data@nsnp, 2)
  expect_equal(rehh_data@nhap, 4)

  stat_ihh <- sumstat_ihh(population = 1, max_snps = 3)
  rehh_data <- stat_ihh$create_rehh_data(seg_sites, 1:4, model)
  expect_equal(dim(rehh_data@haplo), c(4, 3))
  expect_equal(rehh_data@nsnp, 3)

  stat_ihh <- sumstat_ihh(population = 1, max_snps = 2)
  rehh_data <- stat_ihh$create_rehh_data(seg_sites, 1:2, model)
  expect_equal(dim(rehh_data@haplo), c(2, 2))
  expect_equal(rehh_data@nsnp, 2)
  expect_equal(rehh_data@nhap, 2)
})


test_that("calculation of ihh works", {
  skip_if_not_installed("rehh")
  stat_ihh <- sumstat_ihh()
  ihh <- stat_ihh$calculate(list(seg_sites), NULL, NULL, model)
  expect_that(ihh, is_a("data.frame"))
  expect_equal(dim(ihh), c(5, 6))
  expect_equal(ihh$CHR, rep(1, 5))
  expect_equal(ihh$POSITION, pos)

  ihh <- stat_ihh$calculate(list(seg_sites, seg_sites), NULL, NULL,
                            coal_model(4, 2, 337))
  expect_that(ihh, is_a("data.frame"))
  expect_equal(dim(ihh), c(10, 6))
  expect_equal(ihh$CHR, rep(1:2, each = 5))
  expect_equal(ihh$POSITION, c(pos, pos))
})


test_that("calculation of ihs works", {
  skip_if_not_installed("rehh")

  stat_ihh <- sumstat_ihh(calc_ihs = TRUE)
  ihh <- stat_ihh$calculate(list(seg_sites), NULL, NULL, model)
  expect_that(ihh, is_a("list"))
  expect_that(ihh[[1]], is_a("data.frame"))
  expect_that(ihh[[2]], is_a("data.frame"))
  expect_equal(dim(ihh[[1]]), c(5, 6))
  expect_equal(ncol(ihh[[2]]), 3)

  model2 <- coal_model(50, 2, 1000) + feat_mutation(10) + sumstat_seg_sites()
  seg_sites2 <- simulate(model2)$seg_sites

  stat_ihh <- sumstat_ihh(calc_ihs = TRUE)
  ihh <- stat_ihh$calculate(seg_sites2, NULL, NULL, model2)
  expect_that(ihh, is_a("list"))
  expect_that(ihh[[1]], is_a("data.frame"))
  expect_that(ihh[[2]], is_a("data.frame"))
  expect_true(all(ihh$CHR %in% 1:2))
})


test_that("calculation of iHS works with few SNPS", {
  skip_if_not_installed("rehh")

  seg_sites2 <- create_segsites(matrix(c(1, 0, 1, 0), 4, 1), 0.5)
  model2 <- coal_model(4, 1, 337)

  stat_ihh <- sumstat_ihh(calc_ihs = TRUE)
  ihh <- stat_ihh$calculate(list(seg_sites2), NULL, NULL, model2)
  expect_that(ihh, is_a("list"))
  expect_that(ihh[[1]], is_a("data.frame"))
  expect_that(ihh[[2]], is_a("data.frame"))
  expect_equal(dim(ihh[[2]]), c(1, 3))
  expect_true(is.na(ihh[[2]][1, 3]))
})


test_that("ihh works with trios", {
  skip_if_not_installed("rehh")
  model <- model_trios()
  stats <- simulate(model)
  ihh <- sumstat_ihh(population = 1)
  stat <- ihh$calculate(stats$seg_sites, NULL, NULL, model)
  expect_that(stat, is_a("data.frame"))
})


test_that("ihh works with empty segsites", {
  skip_if_not_installed("rehh")
  model <- model_trios()
  seg_sites <- list(seg_sites, create_segsites(matrix(0, 5, 0), numeric(0)))
  ihh <- sumstat_ihh(population = 1)

  stat <- ihh$calculate(seg_sites, NULL, NULL, coal_model(4, 2, 337))
  expect_that(stat, is_a("data.frame"))
  expect_equal(dim(stat), c(5, 6))

  seg_sites <- list(create_segsites(matrix(0, 0, 0), numeric(0)))
  stat <- ihh$calculate(seg_sites, NULL, NULL, model)
  expect_that(stat, is_a("data.frame"))
  expect_equal(dim(stat), c(0, 6))

  seg_sites <- list(create_segsites(matrix(1, 5, 10), 1:10 / 11))
  stat <- ihh$calculate(seg_sites, NULL, NULL, model)
  expect_that(stat, is_a("data.frame"))
  expect_equal(dim(stat), c(0, 6))

  ihh <- sumstat_ihh(population = 1, calc_ihs = TRUE)
  stat <- ihh$calculate(seg_sites, NULL, NULL, model)
  expect_that(stat, is_a("list"))
  expect_that(stat[[1]], is_a("data.frame"))
  expect_that(stat[[2]], is_a("data.frame"))
})
