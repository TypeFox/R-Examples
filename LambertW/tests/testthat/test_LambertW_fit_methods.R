context("Testing MLE_LambertW \n")
set.seed(20)
nobs <- 2 * 1e3

theta.list <- list("s" = list(beta = c(mu = 2, sigma = .1), gamma = -0.1),
                   "h" = list(beta = c(mu = 2, sigma = .1), delta = 0.25),
                   "hh" = list(beta = c(mu = 2, sigma = .1), delta = c(0.1, 0.3)))

yy.list <- lapply(theta.list,
                  function(tt) {
                    yy <- rLambertW(n = nobs, distname = "normal",
                                    theta = tt)
                    yy[yy < 0] <- 1e-4
                    return(yy)
                  })

test_that("plot works", {
  for (tt in names(yy.list)) {
    for (dd in c("t", "normal", "exp")) {
      mod <- MLE_LambertW(yy.list[[tt]], type = tt, distname = dd)
      expect_true(is.list(summary(mod)))
      tmp <- plot(mod)
    }
  } 
})

