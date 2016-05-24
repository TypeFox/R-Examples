context("Priors")

test_that("Check moment dprior_lambda", {
  a <- 3; b <- 7; c <- 4; d <- 5; S <- 10; T <- 13
  expect_that(
    abs(sprior_lambda(a, b, c, d, S, T)$mean - sprior_mu(a,b)$mean*sprior_phi(b,c,d,S,T)$mean) < .Machine$double.eps, 
    is_true()
  )
})

test_that("Compare pprior_lambda with integration of dprior_lambda", {
  a <- 2; b <- 2; c <- 2.5; d <- 2; S <- 10; T <- 10
  p1 <- pprior_lambda(1, a, b, c, d, S, T)
  p2 <- integrate(function(x) dprior_lambda(x, a, b, c, d, S, T), lower=0, upper=1)
  expect_that(
    p2$message=="OK", 
    is_true()
  )
  expect_that(
    abs(p1-p2$value) < .Machine$double.eps^.25, 
    is_true()
  )
})

test_that("Compare prior predictive of x to summation of x_given_y", {
  a <- 3; b <- 1; c <- 10; d <- 20;  T <- 2
  p1 <- dprior_x(x=5, a=a, b=b, c=c, d=d, T=T)
  p2 <- sum(dprior_x_given_y(x=5, y=0:100, a, c, d)*dprior_y(0:100, a, b, T))
  expect_that(
    abs(p1-p2) < .Machine$double.eps^.75, 
    is_true()
  )
  # the same :
  expect_equivalent(p1, p2)
  
})
