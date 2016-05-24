context("Sampling from degenerate spaces")
# Regression test for degenerate sampling spaces (i.e. 0-dimensional)
# Both implicitly and explicitly degenerate spaces failed in 0.4
# See bug #5.

n.samples <- 10

test_that("the explicitly degenerate case is handled", {
  constr <- list(constr=rbind(c(-1,1,0),c(0,1,-1),c(1,1,1),-diag(3)), rhs=c(0,0,1,0,0,0), dir=c("=", "=", "=", "<=", "<=", "<="))
  chain <- hitandrun(constr, n.samples=n.samples)
  expect_equal(chain, matrix(1/3, ncol=3, nrow=n.samples))
})

test_that("the implicitly degenerate case is handled", {
  constr <- list(constr=rbind(c(-1,1,1),c(1,0,-1),c(1,1,1),-diag(3)), rhs=c(0,0,1,0,0,0), dir=c("<=", "<=", "=", "<=", "<=", "<="))
  chain <- hitandrun(constr, n.samples=n.samples)
  expect_equal(chain, matrix(c(0.5, 0, 0.5), ncol=3, nrow=n.samples, byrow=TRUE))
})
