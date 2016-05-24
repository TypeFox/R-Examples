context("Testing bootstrap_LambertW_fit \n")
set.seed(30)
nobs <- 1e2

require(boot)

yy <- rnorm(n = nobs, mean = 10, sd = 5)

est.models <- list("IGMM" = list(fit = IGMM(type = "s", yy),
                                 est.name = "tau"),
                   "MLE" = list(fit = MLE_LambertW(yy, type = "s", 
                                                   distname = "normal"),
                                est.name = "theta"))

for (nn in names(est.models)) {
  cat("Testing", nn, "fit \n")
  
  tmp.model <- est.models[[nn]]$fit
  est.name.tmp <- est.models[[nn]]$est.name
  
  boot.est <- bootstrap.LambertW_fit(tmp.model, R = 10)  
  
  test_that("bootstrap returns correct type of 'boot' output and bootstrap estimates are close to truth", {
    expect_true(inherits(boot.est, "boot"))
    expect_equal(boot.est$t0, 
                 flatten_theta(tmp.model[[est.name.tmp]])[names(boot.est$t0)], 
                 tol = 1e-3)
    
    boot.se <- apply(boot.est$t, 2, sd)
    # bootstrap does have variation (it is indeed subsampling
    # and not using same data every turn)
    expect_true(all(boot.se > 0))
    # avg of bootstrap averages are close to true estimate
    # (use empirical confidence interval)
    expect_true(all(colMeans(boot.est$t) + 2 * boot.se > c(10, 5, 0)))
    expect_true(all(colMeans(boot.est$t) - 2 * boot.se < c(10, 5, 0)))
      
  })

  
  test_that("bootstrap arguments are correct",  {
    expect_error(bootstrap.LambertW_fit("yy", object = tmp.model),
                 info = "first argument must be data")
    expect_error(bootstrap.LambertW_fit(yy, object = tmp.model,
                                        sample.size = -1),
                 info = "sample size can not be negative")
    expect_error(bootstrap.LambertW_fit(yy, object = tmp.model,
                                        sample.size =c(0, 10)))
    expect_error(bootstrap.LambertW_fit(yy, object = tmp.model,
                                        sample.size = length(yy) * 2),
                 info = "sample size can't be larger than original data")
  })
}