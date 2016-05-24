context("convhulln")
test_that("Convhull can run on an example with 3000 points", {
  set.seed(1)
  ps <- matrix(rnorm(3000), ncol=3)
  ps <- sqrt(3)*ps/drop(sqrt((ps^2) %*% rep(1,3)))
  ts <- convhulln(ps)
  expect_that(nrow(ts), equals(1996))
  ts.full <- convhulln(ps, "FA")
  expect_that(ts.full$area, equals(37.47065, tolerance=0.001))
  expect_that(ts.full$vol, equals(21.50165, tolerance=0.001))
})

test_that("convhulln throws an error with duplicated points", {
  load(file.path(system.file(package="geometry"), "extdata", "ordination.Rdata"))
  expect_error(out <- convhulln(ordination))
})

test_that("If the input matrix contains NAs, convhulln should return an error", {
  ps <- matrix(rnorm(999), ncol=3)
  ps <- sqrt(3)*ps/drop(sqrt((ps^2) %*% rep(1,3)))
  ps <- rbind(ps, NA)
  expect_error(convhulln(ps))
})

test_that("If there are not enough points to construct a simplex, an error is thrown", {         
  expect_error(convhulln(diag(4)))
})

test_that("Output to file works", {
  ## To prevent regression in package betapart
  unlink("vert.txt")
  tr <- rbind(c(3,1),c(2,1),c(4,3),c(4,2))
  convhulln(tr, "Fx TO 'vert.txt'")
  expect_true(file.exists("vert.txt"))
  vert <- scan("vert.txt")
  expect_equal(vert, c(4, 2, 1, 0, 3))
})

