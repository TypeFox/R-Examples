context("Check is.theta and rtheta")

theta1 <- rtheta()  # Create a random correctly formatted theta

theta2 <- rtheta(d = 3, m = 5)
theta2$m <- 6  # m is now incoherent with the number of components

theta3 <- rtheta(d = 4, m = 2)
theta3$sigma$comp1[1, 2] <- 0  # Making the covariance matrix non-symmetric

theta4 <- rtheta(d = 10, m = 10)
theta4$sigma$comp1[1, 1] <- 0  # Destroy positive semi-definiteness

theta5 <- rtheta(d = 2, m = 5)
names(theta5)[3] <- "prop"

test_that("rtheta returns proper formatted output (old method)", {
  expect_that(is.logical(is.theta(theta1)), is_true())
  expect_that(length(is.theta(theta1)), equals(1))
  expect_that(is.theta(theta1), is_true())

  suppressWarnings({
    expect_that(is.theta(theta2), is_false())
    expect_that(is.theta(theta3), is_false())
    expect_that(is.theta(theta4), is_false())
    expect_that(is.theta(theta5), is_false())
  })
})

#
# Test all rtheta methods and many d and m
#

set.seed(10)

for (d in 2:10) {
  for (m in 2:10) {
    for (method in as.character(formals(rtheta)$method)[-1]) {
      theta.tmp <- rtheta(m = m, d = d, method = "EqualSpherical")
      test_that(paste("rtheta returns proper formatted output, ",
                      "m=", m, ", d=", d, ", method=", method, sep = ""), {
                        expect_that(is.theta(theta.tmp), is_true())
                      })
    }
  }
}


# Test degenerate input
