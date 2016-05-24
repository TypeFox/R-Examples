context("Testing theta-utils \n")

theta <- list(beta = c(-1, 0.1), 
              delta = c(0.2, 0.2))
theta <- complete_theta(theta)

test_that("theta2unbounded inverse works", {

  theta.restr <- theta2unbounded(theta, distname = "normal")
  # returns again the beta and delta from above
  expect_equivalent(theta2unbounded(theta.restr, inverse = TRUE, 
                                   distname = "normal"),
                    theta)
})
