context("Testing MLE_LambertW \n")
set.seed(20)
nobs <- 2 * 1e3

theta.list <- list("s" = list(beta = c(c = 2, s = .1, df = 4), gamma = -0.1,
                              distname = "t"),
                   "h" = list(beta = c(lambda = 0.5), delta = 0.25,
                              distname = "exp"),
                   "hh" = list(beta = c(mu = 2, sigma = .1), delta = c(0.1, 0.3),
                               distname = "normal"))

yy.list <- lapply(theta.list,
                  function(tt) {
                    rLambertW(n = nobs, distname = tt$distname,
                              theta = tt)
                  })

test_that("total loglikelihood equals input + penalty", {
  for (tt in names(yy.list)) {
    dd.tmp <- theta.list[[tt]]$distname 
    ll <- loglik_LambertW(y = yy.list[[tt]], type = tt, 
                          distname = dd.tmp,
                          theta = theta.list[[tt]])
    # loglikelihood is really sum(log(pdf(x)))
    expect_equal(ll$loglik.LambertW,
                 sum(dLambertW(yy.list[[tt]], distname = dd.tmp,
                               theta = theta.list[[tt]], log = TRUE)))
    if (tt == "s") {
      expect_true(is.na(ll$loglik.input))
      expect_true(is.na(ll$loglik.penalty))
      expect_true(is.numeric(ll$loglik.LambertW))
    }
    if (tt == "h" || tt == "hh") {
      expect_equal(ll$loglik.input + ll$loglik.penalty, ll$loglik.LambertW)
    }
  } 
})


test_that("return.neg returns the negative loglikelihood", {
  for (tt in names(yy.list)) {
    dd.tmp <- theta.list[[tt]]$distname 
    ll <- loglik_LambertW(y = yy.list[[tt]], type = tt, 
                          distname = dd.tmp,
                          theta = theta.list[[tt]],
                          return.negative = FALSE)
    ll.neg <- loglik_LambertW(y = yy.list[[tt]], type = tt, 
                          distname = dd.tmp,
                          theta = theta.list[[tt]],
                          return.negative = TRUE)
    
    expect_equal(ll$loglik.LambertW, -ll.neg)
  }
})


test_that("input loglikelihood can be evaluated on user supplied functions", {
   # evaluating the Gaussian log-likelihood
   yy <- rnorm(100)
   built.in <- loglik_input(beta = c(0, 1), x = yy, distname = "normal") # built-in version
   # or pass your own log pdf function
   user.supplied.log.dX <- 
    loglik_input(beta = c(0, 1), x = yy, distname = "user", 
                 log.dX = function(xx, beta = beta) { 
                   dnorm(xx, mean = beta[1], sd = beta[2], log = TRUE)
                 })
   user.supplied.dX <- 
     loglik_input(beta = c(0, 1), x = yy, distname = "user", 
                  dX = function(x, beta = beta) { 
                    dnorm(x, mean = beta[1], sd = beta[2])
                  })
   
   expect_equal(built.in, user.supplied.log.dX)
   expect_equal(built.in, user.supplied.dX)
})
