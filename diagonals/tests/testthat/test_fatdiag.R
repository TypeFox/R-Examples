# load the package
library(diagonals)


# define context
context("fatdiag: create: scalar")

# create objects
f1 <- fatdiag(12, steps=3)
f2 <- fatdiag(12, size=4)
f3 <- fatdiag(12, size=c(3,4) )
f4 <- fatdiag(12, nrow=12, ncol=4, steps=4)

# test output format
test_that("output format matches", {
  expect_equal( dim(f1), c(12, 12) )
  expect_equal( dim(f2), c(12, 12) )
  expect_equal( dim(f3), c(9, 12) )
  expect_equal( dim(f4), c(12, 4) )
})

# test output values
test_that("output values match", {
  expect_equal( sum(f1), 48)
  expect_equal( sum(f2), 48)
  expect_equal( sum(f3), 36)
  expect_equal( sum(f4), 12)
})


# define context
context("fatdiag: create: vector")

# create objects
f5 <- fatdiag(1:48, steps=3)
f6 <- fatdiag(1:48, size=4)
f7 <- fatdiag(36:1, size=c(3,4) )
f8 <- fatdiag(36:1, nrow=12, ncol=9, steps=3)

# test output format
test_that("output format matches", {
  expect_equal( dim(f5), c(12, 12) )
  expect_equal( dim(f6), c(12, 12) )
  expect_equal( dim(f7), c(9, 12) )
  expect_equal( dim(f4), c(12, 4) )
})

# test output values
test_that("output values match", {
  expect_equal( sum(f5), 1176)
  expect_equal( sum(f6), 1176)
  expect_equal( sum(f7), 666)
  expect_equal( sum(f8), 666)
})


# define context
context("fatdiag: extract")

e1 <- fatdiag(f1, steps=3)
e2 <- fatdiag(f2, steps=4)
e3 <- fatdiag(f3, steps=3)
e4 <- fatdiag(f4, steps=3)

# test output format
test_that("output size matches", {
  expect_equal( dim(e1), NULL )
  expect_equal( dim(e2), NULL )
  expect_equal( dim(e3), NULL )
  expect_equal( dim(e4), NULL )
})

# test output values
test_that("output values match", {
  expect_equal( sum(e1), 48)
  expect_equal( sum(e2), 28)
  expect_equal( sum(e3), 36)
  expect_equal( sum(e4), 6)
})


# define context
context("fatdiag: errors")

# test output values
test_that("missing arguments", {
  #expect_error()
})


# define context
context("fatdiag<-: scalar")

# create objects
fatdiag(f1, steps = 3      ) <- 2
fatdiag(f2, size = 4       ) <- 2
fatdiag(f3, size = c(3,4)  ) <- 2

# test output values
test_that("output values match", {
  expect_equal( sum(f1), 96)
  expect_equal( sum(f6), 1176)
  expect_equal( sum(f3), 72)
  #expect_equal( sum(f8), 36)
})


# define context
context("fatdiag<-: vector")

# create objects
fatdiag(f5, steps = 3      ) <- 1:48
fatdiag(f6, size = c(4,4)  ) <- 1:48
fatdiag(f7, size = c(3,4)  ) <- 1:36
f9 <- diag(16)
fatdiag(f9) <- 64:1

# test output values
test_that("output values match", {
  expect_equal( sum(f5), 1176)
  expect_equal( sum(f6), 1176)
  expect_equal( sum(f7), 666)
  #expect_equal( sum(f8), 36)
  expect_equal( sum(f9), 2080)
})


# define context
context("fatdiag: errors")

# test output values
test_that("wrong arguments", {
  expect_error(fatdiag(10) <- 5)
  expect_error(fatdiag(1:10) <- 5)
  expect_error(fatdiag(array(1:1000, dim = c(10,10,10))) <- 5)
})
