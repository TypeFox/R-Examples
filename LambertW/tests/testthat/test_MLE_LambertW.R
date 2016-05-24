context("Testing MLE_LambertW \n")
set.seed(20)
nobs <- 2 * 1e3
xx <- rnorm(n = nobs, mean = -2, sd = 0.5)

test_that("MLE estimates c(mu, sigma) correctly for a Normal distribution", {
  for (tt in c("s", "h", "hh")) {
    mod <- MLE_LambertW(xx, type = tt, distname = "normal")
    # mean is approx equal
    expect_gt(mod$theta$beta["mu"], -2 - .5 * 2 / sqrt(nobs))
    expect_lt(mod$theta$beta["mu"], -2 + .5 * 2 / sqrt(nobs))
    
    # TODO: replace with actual CI for sigma
    expect_gt(mod$theta$beta["sigma"], .5 - 2 / sqrt(nobs))
    expect_lt(mod$theta$beta["sigma"], .5  + 2 / sqrt(nobs))
    
    # other parameters are zero for Gaussian data
    other.params <- mod$tau[!grepl("mu_x|sigma_x|alpha", names(mod$tau))]
    expect_equal(lp_norm(other.params, 2), 0, tol = 1e-1)
  }
})


theta.list <- list("s" = list(beta = c(mu = 2, sigma = .1), gamma = -0.1),
                   "h" = list(beta = c(mu = 2, sigma = .1), delta = 0.25),
                   "hh" = list(beta = c(mu = 2, sigma = .1), delta = c(0.1, 0.3)))

yy.list <- lapply(theta.list,
                  function(tt) {
                    rLambertW(n = nobs, distname = "normal",
                              theta = tt)
                  })


test_that("MLE estimates are close to truth", {
  for (tt in names(yy.list)) {
    mod <- MLE_LambertW(yy.list[[tt]], type = tt, distname = "normal")
    # mean is approx equal
    expect_gt(mod$theta$beta["mu"], 
                     theta.list[[tt]]$beta["mu"] - 
                       theta.list[[tt]]$beta["sigma"] * 2 / sqrt(nobs))
    expect_lt(mod$theta$beta["mu"], 
                     theta.list[[tt]]$beta["mu"] + 
                       theta.list[[tt]]$beta["sigma"] * 2 / sqrt(nobs))
    
    # TODO: replace with actual CI for sigma
    expect_gt(mod$theta$beta["sigma"], 
              theta.list[[tt]]$beta["sigma"] - 2 / sqrt(nobs))
    expect_lt(mod$theta$beta["sigma"],
              theta.list[[tt]]$beta["sigma"] + 2 / sqrt(nobs))
    
    # other parameters are zero for Gaussian data
    for (nn in names(theta.list[[tt]])) {
      if (nn == "beta") {
        next
      }
      # relative difference
      expect_equal(lp_norm(theta.list[[tt]][[nn]] - mod$theta[[nn]], 2), 0,
                   tol = 1e-1, info = paste0("type ", tt))
    }
  }
})

