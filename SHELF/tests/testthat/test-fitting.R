
test_that("normal distribution fitting works",{
  skip_on_cran()
  m <- 10
  s <- 20
  vals <- c(m - s, m , m + 2 * s)
  myfit <- fitdist(vals, pnorm(vals, m, s ))
  norm.parameters <- unlist(myfit$Normal)
  best.name <- unlist(myfit$best.fitting)
  attributes(best.name) <- attributes(norm.parameters) <- NULL
  expect_equal(norm.parameters, c(m, s))
  expect_equal(best.name, 1)
})

test_that("precision fitting works - normal",{
  skip_on_cran()
  med <- 10
  k <- 1
  # sigma^-2 ~ gamma(a, b)
  a <- 3
  b <- 4
  sigmasq <- 1 / qgamma(c(0.05, 0.95), a, b)
  probs <- pnorm(rep(med + k, 2), med, sigmasq^0.5) - 0.5
  pfit <- fitprecision(c(med, med + k), probs, pplot = F)
  gamma.parameters <- unlist(pfit$Gamma)
  attributes(gamma.parameters) <- NULL
  expect_equal(gamma.parameters, c(a, b), tolerance = 1e-4)
})

test_that("precision fitting works - lognormal",{
  skip_on_cran()
  med <- 10
  k <- 5
  # sigma^-2 ~ gamma(a, b)
  a <- 3
  b <- 4
  sigmasq <- 1 / qgamma(c(0.05, 0.95), a, b)
  probs <- plnorm(rep(med + k, 2), log(med), sigmasq^0.5) - 0.5
  pfit <- fitprecision(interval = c(med, med + k), propvals = probs,
                       trans = "log", pplot = F)
  gamma.parameters <- unlist(pfit$Gamma)
  attributes(gamma.parameters) <- NULL
  expect_equal(gamma.parameters, c(a, b), tolerance = 1e-4)
})
