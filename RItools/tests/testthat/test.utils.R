################################################################################
# Tests for utility functions
################################################################################

library("testthat")

context("Utility Functions")

test_that("xbal formula method", {
  # create a quick xBalance result
  df <- data.frame(Z = rep(c(1,0), 10),
                   X = rnorm(20))

  xb <- xBalance(Z ~ X, data = df)

  # we expect to have a no intercept formula
  expect_equal(Z ~ X, as.formula(xb))
})

test_that("Select a subset of xbal results (for printing, etc)", {
  # create a quick xBalance result
  set.seed(20121129)
  n <- 100
  df <- data.frame(Z = rep(c(1,0), n/2),
                   X = rnorm(n),
                   Y = rnorm(n),
                   W = cut(rnorm(n), breaks = 3),
                   K = as.factor(letters[sample.int(3, n, replace = T)]),
                   U = rnorm(n),
                   S = as.factor(rep(c("A", "B"), each = n/2)))


  # the data are will be grouped as (X, XY) and (all levels W by all levels
  # K), Y is not in an explicit group
  xb <- xBalance(Z ~ X * Y + W * K,
                 data = df,
                 strata = data.frame(none = factor("none"),
                     S = df$S,
                     U = cut(df$U, 3)),
                 report = 'all')

  # strata based subsetion is the easiest as it common across groups and
  # variables
  xb.noneU <- subset(xb, strata = c("none", "U"))

  expect_equivalent(dim(xb$results), c(18, 7, 3))
  expect_equivalent(dim(xb.noneU$results), c(18, 7, 2))
  expect_equal(attr(xb$results, "originals"),
               attr(xb.noneU$results, "originals"))

  # stat and test subsetion only limits the results and overall tables
  # repectively

  xb.zp <- subset(xb, stats = c("z", "p"))
  xb.chip <- subset(xb, tests = c("chisquare", "p.value"))

  expect_equivalent(xb.chip$results, xb$results)

  expect_equivalent(dim(xb.zp$results), c(18, 2, 3))
  expect_equivalent(dim(xb.chip$overall), c(3, 2))

  expect_equal(attr(xb$results, "originals"),
               attr(xb.zp$results, "originals"))

  expect_equal(attr(xb$results, "originals"),
               attr(xb.chip$results, "originals"))

  # for vars and groups, we might want the arguments to interact more directly
  # e.g. subseting vars only returns groups that include the vars
  # e.g. subseting groups only returns the variables contained in those groups

  # for vars, only limit the results table (for now)
  # use exact variable names -- possibly use regexes later
  xb.XY <- subset(xb, vars = c("X", "Y"))

  expect_equivalent(dim(xb.XY$results), c(2, 7, 3))
  expect_equal(length(attr(xb.XY$results, "originals")), 2)


})
