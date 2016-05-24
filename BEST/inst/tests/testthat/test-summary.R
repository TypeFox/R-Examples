
# Tests for summary.BEST

# library(BEST)
# library(testthat)
# test_file("test-summary.R")

context("summary.BEST")

# Fake BEST object (2 groups)
# JAGS returns different values for each run, so can't check for exact values of
#   output for summary, etc., so create a reproducable fake BEST object.
# NB this is actually a matrix, not an mcmc.list object.
set.seed(123)
len <- 1e5
fake <- data.frame(mu1 = rnorm(len, 4.7, 0.47),
              mu2 = rnorm(len, 3.2, 0.39),
              nu = exp(rnorm(len, 3.1, 0.9)),
              sigma1 = exp(rnorm(len, -0.08, 0.42)),
              sigma2 = exp(rnorm(len, -0.28, 0.42)))
class(fake) <- c("BEST", "data.frame")

test_that("summary.BEST with 2 groups and default values gives correct output",  {
  tst <- summary(fake)
  expect_that(class(tst), equals(c("summary.BEST", "matrix")))
  expect_that(colnames(tst),
    equals(c("mean", "median", "mode", "HDI%", "HDIlo", "HDIup","compVal",
             "%>compVal", "ROPElow", "ROPEhigh", "%InROPE")))
  expect_that(rownames(tst),
    equals(c("mu1", "mu2", "muDiff", "sigma1", "sigma2", "sigmaDiff",
             "nu", "log10nu", "effSz")))
  expect_that(round(tst[, "mean"], 5), 
    is_equivalent_to(c(4.70046, 3.20203, 1.49843, 1.00882, 0.82330, 0.18552,
                      33.33190, 1.34608, 1.73125)))
  expect_that(round(tst[, "median"], 5),
    is_equivalent_to(c(4.70044, 3.20177, 1.49882, 0.92136, 0.75354, 0.16123,
                      22.23564, 1.34705, 1.59928)))
  expect_that(round(tst[, "mode"], 5),
    is_equivalent_to(c(4.69105, 3.21097, 1.59130, 0.76825, 0.64755, 0.12847,
                      10.17372, 1.36961, 1.36379)))
  expect_that(tst[, "HDI%"],
    is_equivalent_to(rep(95, 9)))
  expect_that(round(tst[, "HDIlo"], 5),
    is_equivalent_to(c(3.78332, 2.44174, 0.30195, 0.31733, 0.25789, -0.98030,
                       1.04944, 0.57959, 0.10968)))
  expect_that(round(tst[, "HDIup"], 5),
    is_equivalent_to(c(5.62369, 3.97225, 2.70333, 1.88649, 1.53835, 1.34601,
                      97.71139, 2.10913, 3.61654)))
})

test_that("summary.BEST with 2 groups and non-default values gives correct output",  {
  tst <- summary(fake, credMass = 0.8,
    ROPEm = c(-0.1, 0.1), ROPEsd = c(-0.1, 0.1), ROPEeff = c(-0.1, 0.1),
    compValm = 1.5, compValsd = 0, compValeff = 0)
  expect_that(class(tst), equals(c("summary.BEST", "matrix")))
  expect_that(colnames(tst),
    equals(c("mean", "median", "mode", "HDI%", "HDIlo", "HDIup","compVal",
             "%>compVal", "ROPElow", "ROPEhigh", "%InROPE")))
  expect_that(rownames(tst),
    equals(c("mu1", "mu2", "muDiff", "sigma1", "sigma2", "sigmaDiff",
             "nu", "log10nu", "effSz")))
  expect_that(round(tst[, "mean"], 5), 
    is_equivalent_to(c(4.70046, 3.20203, 1.49843, 1.00882, 0.82330, 0.18552,
                      33.33190, 1.34608, 1.73125)))
  expect_that(round(tst[, "median"], 5),
    is_equivalent_to(c(4.70044, 3.20177, 1.49882, 0.92136, 0.75354, 0.16123,
                      22.23564, 1.34705, 1.59928)))
  expect_that(round(tst[, "mode"], 5),
    is_equivalent_to(c(4.69105, 3.21097, 1.59130, 0.76825, 0.64755, 0.12847,
                      10.17372, 1.36961, 1.36379)))
  expect_that(tst[, "HDI%"],
    is_equivalent_to(rep(80, 9)))
  expect_that(round(tst[, "HDIlo"], 5),
    is_equivalent_to(c(4.08736, 2.70201, 0.72404, 0.42834, 0.36181, -0.51840,
                       2.13345, 0.83916, 0.50194)))
  expect_that(round(tst[, "HDIup"], 5),
    is_equivalent_to(c(5.28987, 3.70624, 2.29257, 1.38703, 1.14570, 0.84175,
                      48.01003, 1.84015, 2.65435)))
  notna <- c(3,6,9)
  expect_that(tst[notna, "compVal"],
    is_equivalent_to(c(1.5, 0.0, 0.0)))
  expect_that(tst[-notna, "compVal"],
    is_equivalent_to(rep(NA_real_, 6)))
  expect_that(round(tst[notna, "%>compVal"], 5),
    is_equivalent_to(c(49.922, 63.235, 99.289)))
  expect_that(tst[-notna, "%>compVal"],
    is_equivalent_to(rep(NA_real_, 6)))
  expect_that(tst[notna, "ROPElow"],
    is_equivalent_to(c(-0.1, -0.1, -0.1)))
  expect_that(tst[-notna, "ROPElow"],
    is_equivalent_to(rep(NA_real_, 6)))
  expect_that(tst[notna, "ROPEhigh"],
    is_equivalent_to(c(0.1, 0.1, 0.1)))
  expect_that(tst[-notna, "ROPEhigh"],
    is_equivalent_to(rep(NA_real_, 6)))
  expect_that(round(tst[notna, "%InROPE"], 5),
    is_equivalent_to(c(0.662, 15.724, 0.645)))
  expect_that(tst[-notna, "%InROPE"],
    is_equivalent_to(rep(NA_real_, 6)))
})

# Fake BEST object (1 group)
set.seed(123)
len <- 1e5
fake <- data.frame('mu' = rnorm(len, 1.5, 0.25),
              'nu' = exp(rnorm(len, 3, 1)),
              'sigma' = exp(rnorm(len, -0.75, 0.42)))
class(fake) <- c("BEST", "data.frame")

test_that("summary.BEST with 1 group and default values gives correct output",  {
  tst <- summary(fake)
  expect_that(class(tst), equals(c("summary.BEST", "matrix")))
  expect_that(colnames(tst),
    equals(c("mean", "median", "mode", "HDI%", "HDIlo", "HDIup","compVal",
             "%>compVal", "ROPElow", "ROPEhigh", "%InROPE")))
  expect_that(rownames(tst),
    equals(c("mu", "sigma", "nu", "log10nu", "effSz")))
  expect_that(round(tst[, "mean"], 5), 
    is_equivalent_to(c(1.50024, 0.51587, 33.42182, 1.30515, 3.46960)))
  expect_that(round(tst[, "median"], 5),
    is_equivalent_to(c(1.50024, 0.47274, 20.17708, 1.30486, 3.14097)))
  expect_that(round(tst[, "mode"], 5),
    is_equivalent_to(c(1.49524, 0.39606, 7.64812, 1.31510, 2.56656)))
  expect_that(tst[, "HDI%"],
    is_equivalent_to(rep(95, 5)))
  expect_that(round(tst[, "HDIlo"], 5),
    is_equivalent_to(c(1.01240, 0.16139, 0.48272, 0.45850, 0.95083)))
  expect_that(round(tst[, "HDIup"], 5),
    is_equivalent_to(c(1.99132, 0.96267, 105.88948, 2.16285, 6.73603)))
})

test_that("summary.BEST with 1 group and non-default values gives correct output",  {
  tst <- summary(fake, credMass = 0.8,
    ROPEm = c(-0.1, 0.1), ROPEsd = c(0, 1), ROPEeff = c(-0.1, 0.1),
    compValm = 0, compValsd = 2, compValeff = 0)
  expect_that(class(tst), equals(c("summary.BEST", "matrix")))
  expect_that(colnames(tst),
    equals(c("mean", "median", "mode", "HDI%", "HDIlo", "HDIup","compVal",
             "%>compVal", "ROPElow", "ROPEhigh", "%InROPE")))
  expect_that(rownames(tst),
    equals(c("mu", "sigma", "nu", "log10nu", "effSz")))
  expect_that(round(tst[, "mean"], 5), 
    is_equivalent_to(c(1.50024, 0.51587, 33.42182, 1.30515, 3.46960)))
  expect_that(round(tst[, "median"], 5),
    is_equivalent_to(c(1.50024, 0.47274, 20.17708, 1.30486, 3.14097)))
  expect_that(round(tst[, "mode"], 5),
    is_equivalent_to(c(1.49524, 0.39606, 7.64812, 1.31510, 2.56656)))
  expect_that(tst[, "HDI%"],
    is_equivalent_to(rep(80, 5)))
  expect_that(round(tst[, "HDIlo"], 5),
    is_equivalent_to(c(1.17413, 0.21836, 1.08642, 0.74834, 1.32132)))
  expect_that(round(tst[, "HDIup"], 5),
    is_equivalent_to(c(1.81376, 0.70775, 47.20796, 1.86661, 4.82158)))
  notna <- c(1,2,5)
  expect_that(tst[notna, "compVal"],
    is_equivalent_to(c(0, 2, 0)))
  expect_that(tst[-notna, "compVal"],
    is_equivalent_to(rep(NA_real_, 2)))
  expect_that(round(tst[notna, "%>compVal"], 5),
    is_equivalent_to(c(100.000, 0.037, 100.000)))
  expect_that(tst[-notna, "%>compVal"],
    is_equivalent_to(rep(NA_real_, 2)))
  expect_that(tst[notna, "ROPElow"],
    is_equivalent_to(c(-0.1, 0, -0.1)))
  expect_that(tst[-notna, "ROPElow"],
    is_equivalent_to(rep(NA_real_, 2)))
  expect_that(tst[notna, "ROPEhigh"],
    is_equivalent_to(c(0.1, 1, 0.1)))
  expect_that(tst[-notna, "ROPEhigh"],
    is_equivalent_to(rep(NA_real_, 2)))
  expect_that(round(tst[notna, "%InROPE"], 5),
    is_equivalent_to(c(0.000, 96.358, 0.000)))
  expect_that(tst[-notna, "%InROPE"],
    is_equivalent_to(rep(NA_real_, 2)))
})





