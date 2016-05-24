test_that("SLOPE accepts manually specified noise level", {
  prob = random_problem(20, 10, sigma=2)
  result = SLOPE(prob$X, prob$y, sigma=2)
  expect_equal(result$sigma, 2)
})

test_that("SLOPE estimates sigma from data (when n-p is large)", {
  prob = random_problem(100, 50, amplitude=6, sigma=2)
  result = SLOPE(prob$X, prob$y)
  sigma.hat = result$sigma
  expect_equal(length(sigma.hat), 1)
  expect_true(1.5 < sigma.hat && sigma.hat < 2.5)
})

test_that("SLOPE iteratively estimates sigma from data (when n-p is small)", {
  skip("Known failure")
  prob = random_problem(100, 100, amplitude=6, sigma=2)
  result = SLOPE(prob$X, prob$y)
  sigma.hat = result$sigma
  expect_true(length(sigma.hat) > 1)
  expect_true(1.5 < tail(sigma.hat,1) && tail(sigma.hat,1) < 2.5)
})