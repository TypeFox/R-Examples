context("Logit Integrand")

test_that("Value of Logit Integrand equals", {
  # No random effect
  expect_equal(logit_integrand(0, 
                               X = cbind(rep(1, 2), c(.5, .6)), 
                               A = c(0, 1), 
                               fixed.effects = c(.1, .1), 
                               randomization = .5), 0.078757529007955)
  # With random effect
  expect_equal(logit_integrand(c(0, 1), 
                               X = cbind(rep(1, 2), c(.5, .6)), 
                               A = c(0, 1), 
                               fixed.effects = c(.1, .1), 
                               random.effects = 1,
                               randomization = .5), c(0.0787575290079550, 
                                                      0.0571307956171674))
  # Random effect with allocation within integrand
  expect_equal(logit_integrand(c(0, 1), 
                  X = cbind(rep(1, 2), c(.5, .6)), 
                  A = c(0, 1), 
                  fixed.effects = c(.1, .1), 
                  random.effects = 1,
                  randomization = .5,
                  allocation = .5,
                  integrate.allocation =TRUE ), c(0.315030116031820,
                                                  0.228523182468669))
               
})



