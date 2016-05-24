context("Testing IGMM \n")
set.seed(40)
nobs <- 1e3
yy <- rnorm(n = nobs, mean = 3, sd = 0.2)

test_that("IGMM estimates c(mu, sigma) are approx correct for a Normal distribution", {
  for (tt in c("s", "h", "hh")) {
    cat("Testing IGMM type ", tt, "\n")
    mod <- IGMM(yy, type = tt)
    # mean is approx equal
    expect_gt(mod$tau["mu_x"], 3 - 0.2 * 2 / sqrt(nobs))
    expect_lt(mod$tau["mu_x"], 3 + 0.2 * 2 / sqrt(nobs))
    
    # TODO: replace with actual CI for sigma
    expect_gt(mod$tau["sigma_x"], 0.2 - 2 / sqrt(nobs))
    expect_lt(mod$tau["sigma_x"], 0.2 + 2 / sqrt(nobs))  
    
    other.params <- mod$tau[!grepl("mu_x|sigma_x", names(mod$tau))]
    expect_equal(lp_norm(other.params, 1), 0, tol = 1e-1)
  }
})

yy.neg <- rLambertW(n = 1000, theta = list(beta = c(3, 0.2), gamma = -0.3),
                    distname = "normal")
test_that("IGMM estimate of gamma is negative for negatively skewed", {
  mod <- IGMM(yy.neg, type = "s")
  # mean is approx equal
  expect_gt(mod$tau["mu_x"], 3 - 0.2 * 2 / sqrt(nobs))
  expect_lt(mod$tau["sigma_x"], 3 + 0.2 * 2 / sqrt(nobs))
  
  # TODO: replace with actual CI for sigma
  expect_gt(mod$tau["sigma_x"], 0.2 - 2 / sqrt(nobs))
  expect_lt(mod$tau["sigma_x"], 0.2 + 2 / sqrt(nobs))  
  
  expect_lt(mod$tau["gamma"], -0.2)
})

test_that("IGMM estimate of delta_l > delta_r for negatively skewed", {
  mod <- IGMM(yy.neg, type = "hh")
  # mean is approx equal
  expect_gt(mod$tau["mu_x"], 3 - 0.2 * 3 / sqrt(nobs) - 0.025)
  expect_lt(mod$tau["mu_x"], 3 + 0.2 * 3 / sqrt(nobs) + 0.025)
  
  # TODO: replace with actual CI for sigma
  expect_gt(mod$tau["sigma_x"], 0.2 - 3 / sqrt(nobs))
  expect_lt(mod$tau["sigma_x"], 0.2 + 3 / sqrt(nobs))  
  
  expect_gt(mod$tau["delta_l"], mod$tau["delta_r"])
})


yy.cauchy <- rcauchy(n = nobs)
test_that("IGMM estimate of delta > 1 for Cauchy", {
  mod.cauchy <- IGMM(yy.cauchy, type = "h")
  # mean is approx equal 0
  expect_gt(mod.cauchy$tau["mu_x"], 0 - 2 / sqrt(nobs))
  expect_lt(mod.cauchy$tau["mu_x"], 0 + 2 / sqrt(nobs))

  expect_gt(mod.cauchy$tau["delta"], 0.5)
})

