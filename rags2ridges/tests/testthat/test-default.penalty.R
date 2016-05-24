context("Unit tests of the default.penalty function")

# Random number of sample, random number of classes
p <- 10
n <- replicate(sample(2:5, 1), sample(3:9, 1))
S <- createS(n = n, p = p)
G <- length(S)





test_that("default.penalty accept different input types", {

  # Handling one-way designs
  expect_that(default.penalty(0), is_a("matrix")) # Degenerate case
  expect_that(default.penalty(1), is_a("matrix")) # Degenerate case
  expect_that(default.penalty(2), is_a("matrix"))
  expect_that(default.penalty(4), is_a("matrix"))
  expect_that(default.penalty(sample(1:1000, 1)), is_a("matrix"))

  # Use the length of the list
  expect_that(default.penalty(S), is_a("matrix"))
  expect_that(default.penalty(G = S), is_a("matrix"))

  # Use a data.frame
  df0 <- expand.grid(Factor = c("lvl1", "lvl2"))
  expect_that(default.penalty(df0), is_a("matrix"))
  expect_that(default.penalty(df = df0), is_a("matrix"))

  # Some more elaborate examples
  Slist <- vector("list", 6)
  df1 <- expand.grid(DS = c("DS1", "DS2", "DS3"), ER = c("ER+", "ER-"))

  # Usage (various interface demonstrations)
  expect_that(default.penalty(6, df1, type = "Complete"), is_a("matrix"))
  expect_that(default.penalty(6, type = "CartesianEqual"), gives_warning("df"))
  expect_that(default.penalty(6, df1, type = "CartesianEqual"), is_a("matrix"))
  expect_that(default.penalty(Slist, df1, type = "CartesianEqual"), is_a("matrix"))
  expect_that(default.penalty(6, df1, type = "CartesianUnequal"), is_a("matrix"))
  expect_that(default.penalty(df1), is_a("matrix"))

  # A 2 by 2 by 2 design
  df2 <- expand.grid(A = c("A1", "A2"), B = c("B1", "B2"), C = c("C1", "C3"))
  expect_that(default.penalty(df2), is_a("matrix"))
  expect_that(default.penalty(df2, type = "CartesianEqual"), is_a("matrix"))
  expect_that(default.penalty(df2, type = "CartesianUnequal"), is_a("matrix"))

})
