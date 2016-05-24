context("Kantorovich distance")

test_that("Main example - numeric mode", {
  # unnamed mu and nu
  mu <- c(1/7, 2/7, 4/7)
  nu <- c(1/4, 1/4, 1/2)
  x <- kantorovich(mu, nu)
  expect_equal(x, 0.107142857142857)
  # named mu and nu - same names
  mu <- setNames(mu, c("a","b","c"))
  nu <- setNames(nu, c("a","b","c"))
  x <- kantorovich(mu, nu)
  expect_equal(x, 0.107142857142857)
  # named mu and nu - different names
  mu <- setNames(mu, c("a","b","c"))
  nu <- c(c=1/2, a=1/4, b=1/4)
  x <- kantorovich(mu, nu)
  expect_equal(x, 0.107142857142857)
  # matrix distance - wrong names
  mu <- setNames(c(1/7, 2/7, 4/7), c("a","b","c"))
  nu <- setNames(c(1/4, 1/4, 1/2), c("a","b","c"))
  M <- matrix(1, nrow=3, ncol=3) - diag(3)
  expect_error(kantorovich(mu, nu, dist=M))
  # good names
  rownames(M) <- colnames(M) <- c("a","b","c")
  x <- kantorovich(mu, nu, dist=M)
  expect_equal(x, 0.107142857142857)
  # option details=TRUE
  mu <- c(1/7, 2/7, 4/7)
  nu <- c(1/4, 1/4, 1/2)
  x <- kantorovich(mu, nu, details=TRUE)
  bestj <- attr(x, "joinings")
  expect_equal(length(bestj), 1)
  expect_equal(bestj[[1]], structure(c(0.142857142857143, 0.0357142857142857, 0.0714285714285714,
                                       0, 0.25, 0, 0, 0, 0.5), .Dim = c(3L, 3L), .Dimnames = list(c("1",
                                                                                                    "2", "3"), c("1", "2", "3")))
               )
  #bestj_num <- bestj[[1]]
})

test_that("Main example - bigq mode", {
    # unnamed mu and nu
    mu <- as.bigq(c(1,2,4), 7)
    nu <- as.bigq(c(1,1,1), c(4,4,2))
    x <- kantorovich(mu, nu)
    expect_true(x==as.bigq(3,28))
    # named mu and nu - same names
    mu <- setNames(as.bigq(c(1,2,4), 7), c("a","b","c"))
    nu <- setNames(as.bigq(c(1,1,1), c(4,4,2)), c("a","b","c"))
    x <- kantorovich(mu, nu)
    expect_true(x==as.bigq(3,28))
    # named mu and nu - different names
    mu <- setNames(as.bigq(c(1,2,4), 7), c("a","b","c"))
    nu <- setNames(as.bigq(c(1,1,1), c(2,4,4)), c("c","a","b"))
    x <- kantorovich(mu, nu)
    expect_true(x==as.bigq(3,28))
    # option details=TRUE
    x <- kantorovich(mu, nu, details=TRUE)
    bestj <- attr(x, "joinings")
    expect_equal(length(bestj), 1)
    expect_equal(bestj[[1]], structure(c("1/7", "1/28", "1/14", "0", "1/4", "0", "0", "0", "1/2"),
                                       .Dim = c(3L, 3L),
                                       .Dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
    )
    bestj_num <- attr(kantorovich(as.numeric(mu), as.numeric(nu)[c(2,3,1)], details=TRUE), "joinings")[[1]]
    expect_equal(as.numeric(as.bigq(bestj[[1]])), as.vector(bestj_num))
})

test_that("Main example - character mode", {
  # unnamed mu and nu
  mu <- c("1/7", "2/7", "4/7")
  nu <- c("1/4", "1/4", "1/2")
  x <- kantorovich(mu, nu)
  expect_true(x==as.bigq(3,28))
  # named mu and nu - same names
  mu <- setNames(c("1/7", "2/7", "4/7"), c("a","b","c"))
  nu <- setNames(c("1/4", "1/4", "1/2"), c("a","b","c"))
  x <- kantorovich(mu, nu)
  expect_true(x==as.bigq(3,28))
  # named mu and nu - different names
  mu <- setNames(c("1/7", "2/7", "4/7"), c("a","b","c"))
  nu <- setNames(c("1/2", "1/4", "1/4"), c("c","a","b"))
  x <- kantorovich(mu, nu)
  expect_true(x==as.bigq(3,28))
  # option details=TRUE
  x <- kantorovich(mu, nu, details=TRUE)
  bestj <- attr(x, "joinings")
  expect_equal(length(bestj), 1)
  expect_equal(bestj[[1]], structure(c("1/7", "1/28", "1/14", "0", "1/4", "0", "0", "0", "1/2"),
                                     .Dim = c(3L, 3L),
                                     .Dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  )
  bestj_bigq <- attr(kantorovich(as.bigq(mu), as.bigq(nu)[c(2,3,1)], details=TRUE), "joinings")[[1]]
  expect_equal(as.numeric(as.bigq(bestj[[1]])), as.numeric(as.bigq(bestj_bigq)))
  expect_equal(rownames(bestj[[1]]), names(mu))
  expect_false(all(colnames(bestj[[1]])==names(nu)))
})

test_that("Non-square example - numeric mode", {
  # unnamed mu and nu
  mu <- c(2/5,3/5)
  nu <- c(1/4,1/4,1/4,1/4)
  x <- kantorovich(mu, nu)
  expect_equal(x, 0.5)
  # named mu and nu - same names
  mu <- setNames(mu, c("a","b"))
  nu <- setNames(nu, c("a","b","c","d"))
  x <- kantorovich(mu, nu)
  expect_equal(x, 0.5)
  # named mu and nu - different names
  mu <- setNames(mu, c("b","a"))
  x <- kantorovich(mu, nu)
  expect_equal(x, 0.5)
  # matrix distance
  mu <- setNames(c(2/5,3/5), c("a","b"))
  nu <- setNames(c(1/4,1/4,1/4,1/4), c("a","b","c","d"))
  M <- matrix(1, nrow=4, ncol=4) - diag(4)
  rownames(M) <- c("a","b","c","d"); colnames(M) <- names(nu)
  x <- kantorovich(mu, nu)
  expect_equal(x, 0.5)
})

test_that("Non-square example - bigq mode", {
    # unnamed mu and nu
    mu <- as.bigq(c(2/5,3/5))
    nu <- as.bigq(c(1/4,1/4,1/4,1/4))
    x <- kantorovich(mu, nu)
    expect_identical(x, as.bigq(0.5))
    # named mu and nu - same names
    mu <- setNames(mu, c("a","b"))
    nu <- setNames(nu, c("a","b","c","d"))
    x <- kantorovich(mu, nu)
    expect_identical(x, as.bigq(0.5))
    # named mu and nu - different names
    mu <- setNames(mu, c("b","a"))
    x <- kantorovich(mu, nu)
    expect_identical(x, as.bigq(0.5))
})
