
# test_that code for hdi generic


context("Highest Density Interval (HDI)")

test_that("hdi.default gives correct output",  {
  set.seed(123)
  tst <- rnorm(1e4)
  expect_that(names(hdi(tst)), equals(c("lower", "upper")))
  expect_that(round(hdi(tst), 6), is_equivalent_to(c(-2.02611, 1.88284)))
  expect_that(hdi(5), is_equivalent_to(c(5, 5)))
  expect_that(round(hdi(tst, 0.6), 6), is_equivalent_to(c(-0.826547, 0.851867)))
  expect_that(round(hdi(tst, 0.9999999), 6), is_equivalent_to(c(-3.845320, 3.847768)))
  tst[5:105] <- NA
  expect_that(round(hdi(tst), 6), is_equivalent_to(c(-2.029857, 1.881446)))
  expect_that(hdi(rep(NA_real_, 1e4)), is_equivalent_to(rep(NA_real_, 2)))
  expect_that(hdi(numeric(0)), is_equivalent_to(rep(NA_real_, 2)))
  tstInf <- c(rep(-Inf, 1e3), rep(Inf, 1e3))
  expect_that(hdi(tstInf), is_equivalent_to(c(-Inf, Inf)))
  expect_that(hdi(letters), throws_error("No applicable method for class"))
}  )


test_that("hdi.matrix gives correct results", {
  set.seed(123)
  len <- 1e5
  fake <- cbind('mu' = rnorm(len, 4.7, 0.47),
                #'nu' = exp(rnorm(len, 3.1, 0.9)),
                'nu' = NA,
                'sigma' = exp(rnorm(len, -0.28, 0.42)))
  fakeres <- hdi(fake)
  expect_that(dim(fakeres), equals(c(2, 3)))
  expect_that(rownames(fakeres), equals(c("lower", "upper")))
  expect_that(colnames(fakeres), equals(c("mu", "nu", "sigma")))
  expect_that(round(fakeres[1, ], 6), is_equivalent_to(c(3.783317, NA_real_, 0.259514)))
  expect_that(round(fakeres[2, ], 6), is_equivalent_to(c(5.623685, NA_real_, 1.554167)))
  expect_that(round(hdi(fake, 0.6)[1, ], 6), is_equivalent_to(c(4.313969, NA_real_, 0.432604)))
  expect_that(hdi(fake, -1), throws_error("credMass must be between 0 and 1"))
  na_mat <- matrix(NA_real_, 1e4, 3)
  expect_that(hdi(na_mat), is_equivalent_to(matrix(NA_real_, 2, 3)))
} )

test_that("hdi.data.frame gives correct results", {
  set.seed(123)
  len <- 100
  fake <- data.frame('mu' = rnorm(len, 4.7, 0.47),
                'sigma' = exp(rnorm(len, -0.28, 0.42)),
                'nu' = NA,
                'letters' = letters[1:25],
                'factor' = as.factor(1:5), stringsAsFactors=FALSE)
  # str(fake)
  fakeres <- hdi(fake)
  expect_that(dim(fakeres), equals(c(2, 5)))
  expect_that(rownames(fakeres), equals(c("lower", "upper")))
  expect_that(colnames(fakeres), equals(c("mu", "sigma", "nu", "letters", "factor")))
  expect_that(round(fakeres[1, 1:2], 6), is_equivalent_to(c(4.105264, 0.319064)))
  expect_that(round(fakeres[2, 1:2], 6), is_equivalent_to(c(5.728047, 1.685098)))
  expect_that(as.vector(fakeres[, 3:5]), is_equivalent_to(rep(NA_real_, 6)))
  expect_that(as.vector(round(hdi(fake, 0.6)[, 1:2], 6)),
      is_equivalent_to(c(4.366676, 5.029637, 0.383080, 0.786412)))
  expect_that(hdi(fake, -1), throws_error("credMass must be between 0 and 1"))
} )


# How to test hdi.mcmc.list?

test_that("hdi.function gives correct results", {
  funcres <- hdi(qgamma, shape=2.5, rate=2)
  expect_that(names(funcres), equals(c("lower", "upper")))
  expect_that(round(funcres, 6),
    is_equivalent_to(c(0.074060, 2.797866)))
  expect_that(round(hdi(qgamma, 0.6, shape=2.5, rate=2), 6),
    is_equivalent_to(c(0.328839, 1.432596)))
  expect_that(hdi(qbeta, shape=2.5, rate=2),
    throws_error("Incorrect arguments for the inverse cumulative density function qbeta"))
} )

test_that("hdi gives correct error messages",  {
  set.seed(123)
  tst <- rnorm(1e4)
  expect_that(hdi(tst, 0), throws_error("credMass must be between 0 and 1"))
  expect_that(hdi(tst, 1), throws_error("credMass must be between 0 and 1"))
  expect_that(hdi(tst, NA), throws_error("credMass must be between 0 and 1"))
  expect_that(hdi(tst, (1:3)/5), throws_error("credMass must be between 0 and 1"))
  expect_that(hdi(tst, -1), throws_error("credMass must be between 0 and 1"))
  expect_that(hdi(tst, 2), throws_error("credMass must be between 0 and 1"))
  tstInf <- c(rep(-Inf, 1e3), rep(Inf, 1e3))
  expect_that(hdi(tstInf), is_equivalent_to(c(-Inf, Inf)))
  expect_that(hdi(letters), throws_error("No applicable method for class"))
}  )


