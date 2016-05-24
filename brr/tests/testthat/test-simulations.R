context("Simulations")

test_that("Check PGB2 using simulations",{
  set.seed(666)
  a <- 2 ; c <- 5 ; d <- 30; tau <- 6
  sims <- rPGB2(2e6, a, c, d, tau)
  p1 <- ecdf(sims)(1:3)
  p2 <- pPGB2(1:3, a, c, d, tau)
  expect_equal(p1, p2, tolerance=1e-3)
  s <- summary_PGB2(a, c, d, tau)
  expect_equal(mean(sims), s$mean, tolerance=1e-3)
  expect_equal(sd(sims), s$sd, tolerance=1e-3)
  #
  a <- 2 ; c <- 5 ; d <- 30; tau <- 2
  expect_true(all(dPGB2(0:10, a, c, d, tau=1)==dbeta_nbinom(0:10, a, c, d)))
  tau <- 2
  nsims <- 1e6
  sims <- rpois(nsims, rgamma(nsims, a, rbeta2(nsims, c, d, scale=tau)))
  summ <- summary_PGB2(a,c,d,tau)
  p <- pPGB2(12, a, c, d, tau)
  expect_equal(p, length(sims[sims<=12])/nsims, tolerance=1e-3)
  expect_equal(mean(sims), summ$mean, tolerance=1e-3)
  expect_equal(sd(sims), summ$sd, tolerance=0.02)
})

test_that("Check PGIB using simulations",{
  set.seed(666)
  a <- 2 ; alpha <- 5 ; beta <- 30; rho <- 0.6
  sims <- rPGIB(1e6, a, alpha, beta, rho)
  p1 <- ecdf(sims)(1:3)
  p2 <- pPGIB(1:3, a,  alpha, beta, rho)
  expect_equal(p1, p2, tolerance=1e-3, scale=1)
  s <- summary_PGIB(a, alpha, beta, rho)
  expect_equal(mean(sims), s$mean, tolerance=1e-3)
  expect_equal(sd(sims), s$sd, tolerance=1e-2)
})

test_that("Check GIB using siumations",{
  set.seed(666)
  a <- 50; b <- 1; c <- 5; d <- 4
  sims <- rGIB(1e6,a,c,d,b)
  summ <- summary_GIB(a=a, alpha=c, beta=d, rho=b)
  expect_equal(mean(sims), summ$mean, tolerance=1e-3)
  expect_equal(sd(sims), summ$sd, tolerance=1e-3)
  #
  a <- 3; alpha <- 4; beta <- 2; rho <- 2.5
  sims <- rgamma(1e6, a, rho/rbeta(1e6,beta,alpha))
  summ <- summary_GIB(a, alpha, beta, rho)
  expect_equal(mean(sims), summ$mean, tolerance=1e-3)
  expect_equal(sd(sims), summ$sd, tolerance=1e-3)
})

test_that("Check BNB using simulations", {
  require(magrittr)
  set.seed(666)
  a <- 2 ; c <- 5 ; d <- 30
  nsims <- 1e6
  sims <- rpois(nsims, rgamma(nsims, a, rbeta2(nsims, c, d, scale=1)))
  summ <- summary_beta_nbinom(a, c, d)
  expect_equal(mean(sims), summ$mean, tolerance=1e-3)
  expect_equal(sd(sims), summ$sd, tolerance=0.05)
  expect_equal(pbeta_nbinom(12, a, c, d), ecdf(sims)(12), tolerance=0.001)
})

test_that("Check GB2 using simulations", {
  set.seed(666)
  a <- 2 ; c <- 4 ; d <- 3-1e-5 
  tau <- 20/12
  nsims <- 1e6
  sims <- rGB2(nsims, a, c, d, tau)
  p <-  pGB2(1, a, c, d, tau)
  expect_equal(p, length(sims[sims<=1])/nsims, tolerance=1e-3)
  I <- integrate(function(x) dGB2(x, a, c, d, tau), lower=0, upper=1)
  expect_equal(p, I$value, tolerance=1e-5)
})

test_that("Check dprior_lambda using simulations", {
  set.seed(666)
  a <- 2 ; b <- 2 ; c <- 2.5 ; d <- 3 ; S <- T <- 10
  nsims <- 1e6
  sims <- rprior_lambda(nsims, a, b, c, d, S, T)
  summ <- sprior_lambda(a,b,c,d,S,T)
  expect_equal(mean(sims), summ$mean, tolerance=0.05)
  expect_equal(sd(sims), summ$sd, tolerance=0.05)
  expect_equal(pprior_lambda(1, a,b,c,d,S,T), ecdf(sims)(1), tolerance=0.001)
})
