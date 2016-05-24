################################################################################
# Tests for xBalance function
################################################################################
## if working interactively in inst/tests you'll need
## library(RItools, lib.loc = '../../.local')
library("testthat")

context("xBalance Functions")

test_that("xBalance returns covariance of tests", {
  set.seed(20130801)
  n <- 500

  library(MASS)
  xs <- mvrnorm(n,
                mu = c(1,2,3),
                Sigma = matrix(c(1, 0.5, 0.2,
                    0.5, 1, 0,
                    0.2, 0, 1), nrow = 3, byrow = T))

  p <- plogis(xs[,1]- 0.25 * xs[,2] - 1)
  z <- rbinom(n, p = p, size = 1)
  s <- rep(c(0,1), each = n/2)

  dat <- cbind(z, xs, s)


  res <- xBalance(z ~ . - s,
                  data = as.data.frame(dat),
                  report = 'all',
                  strata = list("Unadj" = NULL,
                                "s"   = ~ s))

  tcov <- attr(res$overall, "tcov")

  expect_false(is.null(tcov))

  expect_equal(length(tcov), 2)
  expect_equal(dim(tcov[[1]]), c(3,3))

  # variance should be the squares of the reported null SDs
  expect_equal(sqrt(diag(tcov[[1]])), res$results[, "adj.diff.null.sd", 1])
  expect_equal(sqrt(diag(tcov[[2]])), res$results[, "adj.diff.null.sd", 2])

})

test_that("strata in formula", {
  set.seed(20130801)
  n <- 500

  library(MASS)
  xs <- mvrnorm(n,
                mu = c(1,2,3),
                Sigma = matrix(c(1, 0.5, 0.2,
                                 0.5, 1, 0,
                                 0.2, 0, 1), nrow = 3, byrow = T))
  colnames(xs) <- c("X1", "X2", "X3")

  p <- plogis(xs[,1]- 0.25 * xs[,2] - 1)
  z <- rbinom(n, p = p, size = 1)
  s <- rep(c(0,1), each = n/2)
  s2 <- rep(rep(c(0,1), each=n/4), 2)

  dat <- cbind(z, xs, s, s2)


  res <- xBalance(z ~ . - s - s2,
                  data = as.data.frame(dat),
                  report = 'all',
                  strata = list("Unadj" = NULL,
                                "s"     = ~ s))

  res2 <- xBalance(z ~ . - s - s2 + strata(s),
                   data = as.data.frame(dat),
                   report = 'all')

  expect_true(all.equal(res, res2, check.attributes=FALSE))



  res3 <- xBalance(z ~ . - s - s2,
                   data = as.data.frame(dat),
                   report = 'all',
                   strata = list("Unadj" = NULL,
                                 "s"     = ~ s,
                                 "s2"    = ~ s2))

  res4 <- xBalance(z ~ . - s - s2 + strata(s) + strata(s2),
                   data = as.data.frame(dat),
                   report = 'all')

  expect_true(all.equal(res3, res4, check.attributes=FALSE))


})

test_that("partial arguments to report", {
  data(nuclearplants)

  expect_error(xBalance(pr ~ ., data=nuclearplants, report = "a"), "multiple")
  expect_error(xBalance(pr ~ ., data=nuclearplants, report = "b"), "Invalid")
  expect_error(xBalance(pr ~ ., data=nuclearplants, report = "adj.mean"), "multiple")

  # just to test these don't error
  res <- xBalance(pr ~ ., data=nuclearplants, report = "adj.means")
  res <- xBalance(pr ~ ., data=nuclearplants, report = "adj.mean.diffs")
  res <- xBalance(pr ~ ., data=nuclearplants, report = "adj.mean.diffs.null.sd")

  # everything should be identical
  res.z1 <- xBalance(pr ~ ., data=nuclearplants, report = "z.scores")
  res.z2 <- xBalance(pr ~ ., data=nuclearplants, report = "z")
  expect_true(identical(res.z1, res.z2))

  res.chi1 <- xBalance(pr ~ ., data=nuclearplants, report = "chisquare.test")
  res.chi2 <- xBalance(pr ~ ., data=nuclearplants, report = "chi")
  expect_true(identical(res.chi1, res.chi2))

  res.std.d1 <- xBalance(pr ~ ., data=nuclearplants, report = "std.diffs")
  res.std.d2 <- xBalance(pr ~ ., data=nuclearplants, report = "std.d")
  expect_true(identical(res.std.d1, res.std.d2))

  res.a.m.d.n1 <- xBalance(pr ~ ., data=nuclearplants, report = "adj.mean.diffs.null.sd")
  res.a.m.d.n2 <- xBalance(pr ~ ., data=nuclearplants, report = "adj.mean.diffs.n")
  expect_true(identical(res.a.m.d.n1, res.a.m.d.n2))

  res.mult1 <- xBalance(pr ~ ., data=nuclearplants,
                        report = c("adj.means", "z.scores", "chisquare.test", "p.values", "adj.mean.diffs", "adj.mean.diffs.null.sd"))
  res.mult2 <- xBalance(pr ~ ., data=nuclearplants,
                        report = c("adj.means", "z", "chi", "p", "adj.mean.diffs", "adj.mean.diffs.null"))
  expect_true(identical(res.mult1, res.mult2))

  res.all1 <- xBalance(pr ~ ., data=nuclearplants, report = "al")
  res.all2 <- xBalance(pr ~ ., data=nuclearplants, report = "all")
  expect_true(identical(res.all1, res.all2))

  # let's make sure the outputs are what we expect

  # Only z and p (p is always returned)
  expect_true(all(colnames(res.z1$results) == c("z", "p")))
  expect_true(all(colnames(res.z2$results) == c("z", "p")))

  # `results` is empty; overall isn't
  expect_true(is.null(colnames(res.chi1$results)))
  expect_true(is.null(colnames(res.chi2$results)))
  expect_true(!is.null(colnames(res.chi1$overall)))
  expect_true(!is.null(colnames(res.chi2$overall)))

  # Only z and p (p is always returned)
  expect_true(all(colnames(res.a.m.d.n1$results) == c("adj.diff.null.sd", "p")))
  expect_true(all(colnames(res.a.m.d.n2$results) == c("adj.diff.null.sd", "p")))

  expect_true(all(colnames(res.mult1$results) == c("pr=0", "pr=1", "adj.diff", "adj.diff.null.sd", "z", "p")))
  expect_true(all(colnames(res.mult2$results) == c("pr=0", "pr=1", "adj.diff", "adj.diff.null.sd", "z", "p")))
  expect_true(!is.null(colnames(res.chi1$overall)))
  expect_true(!is.null(colnames(res.chi2$overall)))

})

test_that("Passing post.alignment.transform, #26", {
  data(nuclearplants)

  # Identity shouldn't have an effect
  res1 <- xBalance(pr ~ ., data=nuclearplants)
  res2 <- xBalance(pr ~ ., data=nuclearplants, post.alignment.transform = function(x) x)

  expect_true(all.equal(res1, res2)) ## allow for small numerical differences

  res3 <- xBalance(pr ~ ., data=nuclearplants, post.alignment.transform = rank)

  expect_true(all(dim(res1$results) == dim(res3$results)))

  expect_error(xBalance(pr ~ ., data=nuclearplants, post.alignment.transform = mean),
               "Invalid post.alignment.transform given")

  res4 <- xBalance(pr ~ ., data=nuclearplants, post.alignment.transform = rank, report="all")
  res5 <- xBalance(pr ~ ., data=nuclearplants, report="all")

  expect_false(isTRUE(all.equal(res4,res5)))

  # a wilcoxon rank sum test, asymptotic and w/o continuity correction
  res6 <- xBalance(pr ~ cost, data=nuclearplants, post.alignment.transform = rank, report="all")

  expect_equal(res6$results["cost", "p", "unstrat"],
               wilcox.test(cost~pr, data=nuclearplants, exact=FALSE, correct=FALSE)$p.value)

  # w/ one variable, chisquare p value should be same as p value on that variable
  expect_equal(res6$results["cost", "p", "unstrat"],
               res6$overall["unstrat","p.value"])

  # to dos: test combo of a transform with non-default stratum weights.

})
