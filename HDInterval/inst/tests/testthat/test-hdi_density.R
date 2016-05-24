
# test_that code for hdi.density


context("hdi.density")

test_that("hdi.density gives correct output",  {
  set.seed(123)
  tst <- density(rnorm(1e4))
  # with allowSplit = FALSE, the default
  expect_that(names(hdi(tst)), equals(c("lower", "upper")))
  expect_that(round(hdi(tst), 6), is_equivalent_to(c(-1.997715, 1.949980)))
  expect_that(round(hdi(tst, 0.6), 6), is_equivalent_to(c(-0.826789, 0.862691)))
  expect_that(round(hdi(tst, 0.9999999), 6), is_equivalent_to(c(-4.255930, 4.275105)))
  expect_that(hdi(tst, 0), throws_error("credMass must be between 0 and 1"))
  # with allowSplit = TRUE, but unimodal density (values should be same as above,
  #   but it's a matrix)
  expect_that(colnames(hdi(tst, allowSplit=TRUE)), equals(c("begin", "end")))
  expect_that(as.vector(round(hdi(tst, allowSplit=TRUE), 6)), is_equivalent_to(c(-1.997715, 1.949980)))
  expect_that(as.vector(round(hdi(tst, 0.6, allowSplit=TRUE), 6)), is_equivalent_to(c(-0.826789, 0.862691)))
  expect_that(as.vector(round(hdi(tst, 0.9999999, allowSplit=TRUE), 6)), is_equivalent_to(c(-4.255930, 4.275105)))
  # With bimodal density and allowSplit = FALSE (default)
  tst2 <- density(c(rnorm(1e5), rnorm(5e4, 7)))
  expect_warning(hdi(tst2), "The HDI is discontinuous but allowSplit = FALSE")
  expect_warning(expect_that(names(hdi(tst2)), equals(c("lower", "upper"))))
  expect_warning(expect_that(round(hdi(tst2), 6), is_equivalent_to(c(-1.924750, 8.458657))))
  expect_warning(expect_that(round(hdi(tst2, 0.6), 6), is_equivalent_to(c(-1.740528, 1.709442))))
  expect_that(round(hdi(tst2, 0.9999999), 6), is_equivalent_to(c(-5.073267, 11.707659)))
  # With bimodal density and allowSplit = TRUE
  expect_that(dim(hdi(tst2, allowSplit=TRUE)), equals(c(2, 2)))
  expect_that(colnames(hdi(tst2, allowSplit=TRUE)), equals(c("begin", "end")))
  expect_that(as.vector(round(hdi(tst2, allowSplit=TRUE), 6)),
    is_equivalent_to(c(-2.192709,  5.209656,  2.195118,  8.827101)))
  expect_that(as.vector(round(hdi(tst2, 0.6, allowSplit=TRUE), 6)), 
    is_equivalent_to(c(-1.254853,  6.649935,  1.257262,  7.286337)))
  expect_that(as.vector(round(hdi(tst2, 0.9999999, allowSplit=TRUE), 6)),
    is_equivalent_to(c(-5.073267, 11.707659)))
}  )

