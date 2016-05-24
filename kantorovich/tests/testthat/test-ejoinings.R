context("Joinings")

test_that("ejoinings in numeric mode returns a list of numeric matrices", {
  mu <- nu <- c(0.5, 0.5)
  x <- ejoinings(mu, nu)
  expect_is(x, "list")
  expect_true(length(x)==2)
  expect_equal(x[[1]], structure(c(0, 0.5, 0.5, 0), .Dim = c(2L, 2L), .Dimnames = list(
    c("1", "2"), c("1", "2"))))
  expect_equal(x[[2]], structure(c(0.5, 0, 0, 0.5), .Dim = c(2L, 2L), .Dimnames = list(
    c("1", "2"), c("1", "2"))))
  # named mu and nu with same names
  mu <- nu <- c(a=0.5, b=0.5)
  x <- ejoinings(mu, nu)
  expect_is(x, "list")
  expect_true(length(x)==2)
  expect_equal(x[[1]], structure(c(0, 0.5, 0.5, 0), .Dim = c(2L, 2L), .Dimnames = list(
    c("a", "b"), c("a", "b"))))
  expect_equal(x[[2]], structure(c(0.5, 0, 0, 0.5), .Dim = c(2L, 2L), .Dimnames = list(
    c("a", "b"), c("a", "b"))))
  # named mu and nu with different names
  mu <- c(a=0.5, b=0.5); nu <- c(b=0.5, a=0.5)
  x <- ejoinings(mu, nu)
  expect_is(x, "list")
  expect_true(length(x)==2)
  expect_equal(x[[2]], structure(c(0, 0.5, 0.5, 0), .Dim = c(2L, 2L), .Dimnames = list(
    c("a", "b"), c("a", "b"))))
  expect_equal(x[[1]], structure(c(0.5, 0, 0, 0.5), .Dim = c(2L, 2L), .Dimnames = list(
    c("a", "b"), c("a", "b"))))
})

test_that("ejoinings in bigq mode returns a list of character matrices", {
  mu <- nu <- as.bigq(c(0.5,0.5))
  x <- ejoinings(mu, nu)
  expect_equal(x[[1]], structure(c("0", "1/2", "1/2", "0"), .Dim = c(2L, 2L), .Dimnames = list(
    c("1", "2"), c("1", "2"))))
  expect_equal(x[[2]], structure(c("1/2", "0", "0", "1/2"), .Dim = c(2L, 2L), .Dimnames = list(
    c("1", "2"), c("1", "2"))))
  # named mu and nu with same names
  mu <- nu <- setNames(as.bigq(c(0.5,0.5)), c("a", "b"))
  x <- ejoinings(mu, nu)
  expect_equal(x[[1]], structure(c("0", "1/2", "1/2", "0"), .Dim = c(2L, 2L), .Dimnames = list(
    c("a", "b"), c("a", "b"))))
  expect_equal(x[[2]], structure(c("1/2", "0", "0", "1/2"), .Dim = c(2L, 2L), .Dimnames = list(
    c("a", "b"), c("a", "b"))))
  # named mu and nu with different names
  nu <- setNames(as.bigq(c(0.5,0.5)), c("b", "a"))
  x <- ejoinings(mu, nu)
  expect_equal(x[[2]], structure(c("0", "1/2", "1/2", "0"), .Dim = c(2L, 2L), .Dimnames = list(
    c("a", "b"), c("a", "b"))))
  expect_equal(x[[1]], structure(c("1/2", "0", "0", "1/2"), .Dim = c(2L, 2L), .Dimnames = list(
    c("a", "b"), c("a", "b"))))
})

test_that("Main example", {
  mu <- c(1/7,2/7,4/7)
  nu <- c(1/4,1/4,1/2)
  x <- ejoinings(mu, nu)
  expect_true(length(x)==15)
  expect_equal(x[[1]], structure(c(0.142857142857143, 0, 0.107142857142857, 0, 0, 0.25,
                                   0, 0.285714285714286, 0.214285714285714), .Dim = c(3L, 3L), .Dimnames = list(
                                     c("1", "2", "3"), c("1", "2", "3"))))
  #
  if(require(gmp)){
    mu <- as.bigq(c(1,2,4),7)
    nu <- as.bigq(c(1,1,1),c(4,4,2))
    x <- ejoinings(mu, nu)
    expect_true(length(x)==15)
    expect_equal(x[[1]], structure(c("1/7", "0", "3/28", "0", "0", "1/4", "0", "2/7",
                                     "3/14"), .Dim = c(3L, 3L), .Dimnames = list(c("1", "2", "3"),
                                                                                 c("1", "2", "3"))))
  }
})

test_that("Non-square example - with zeros", {
  mu <- c(2/5,3/5)
  nu <- c(1/4,1/4,1/4,1/4)
  joinings <- ejoinings(mu, nu, zeros=TRUE)
  expect_true(length(joinings)==12)
  expect_true(all(sapply(lapply(joinings, colSums), function(x) all.equal(x, nu, check.names=FALSE))))
  expect_true(all(sapply(lapply(joinings, rowSums), function(x) all.equal(x, c(mu,0,0), check.names=FALSE))))
  #
  mu <- as.bigq(c(2,3), 5)
  nu <- as.bigq(nu)
  joinings <- ejoinings(mu, nu, zeros=TRUE)
  expect_true(length(joinings)==12)
  expect_true(all(sapply(lapply(joinings, function(x) apply.bigq(t(as.matrix(as.bigq(x))), 2, sum)), function(x) all.equal(x, nu, check.names=FALSE))))
  expect_true(all(sapply(lapply(joinings, function(x) apply.bigq(as.matrix(as.bigq(x)), 1, sum)), function(x) all.equal(x, c(mu,0,0), check.names=FALSE))))
  #
  mu <- setNames(as.bigq(c(1,2,4), 7), c("a", "b", "c"))
  nu <- setNames(as.bigq(c(3,1), 4), c("b", "c"))
  joinings <- ejoinings(mu, nu, zeros=TRUE)
  expect_true(length(joinings)==4)
  expect_identical(joinings[[4]], structure(c("0", "0", "0", "0", "2/7", "13/28", "1/7", "0", "3/28"
  ), .Dim = c(3L, 3L), .Dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
})

test_that("Non-square example - without zeros", {
  mu <- c(2/5,3/5)
  nu <- c(1/4,1/4,1/4,1/4)
  joinings <- ejoinings(mu, nu)
  expect_true(length(joinings)==12)
  expect_true(all(sapply(lapply(joinings, colSums), function(x) all.equal(x, nu, check.names=FALSE))))
  expect_true(all(sapply(lapply(joinings, rowSums), function(x) all.equal(x, mu, check.names=FALSE))))
  #
  mu <- as.bigq(c(2,3), 5)
  nu <- as.bigq(nu)
  joinings <- ejoinings(mu, nu)
  expect_true(length(joinings)==12)
  expect_true(all(sapply(lapply(joinings, function(x) apply.bigq(t(as.matrix(as.bigq(x))), 2, sum)), function(x) all.equal(x, nu, check.names=FALSE))))
  expect_true(all(sapply(lapply(joinings, function(x) apply.bigq(as.matrix(as.bigq(x)), 1, sum)), function(x) all.equal(x, mu, check.names=FALSE))))
  #
  mu <- setNames(as.bigq(c(1,2,4), 7), c("a", "b", "c"))
  nu <- setNames(as.bigq(c(3,1), 4), c("b", "c"))
  joinings <- ejoinings(mu, nu)
  expect_true(length(joinings)==4)
  expect_identical(joinings[[4]], structure(c("0", "2/7", "13/28", "1/7", "0", "3/28"), .Dim = c(3L, 2L), .Dimnames = list(c("a", "b", "c"), c("b", "c"))))
})
