context("iroundrobin iterator")

test_that("iroundrobin works with a single iterable object", {
  rr <- iroundrobin(1:5)
  expect_equal(take(rr, 5), as.list(1:5))
  expect_error(iterators::nextElem(rr), "StopIteration")
})

test_that("iroundrobin functions properly with two iterators of same type", {
  it <- iterators::iter(1:5)
  it2 <- iterators::iter(1:5)
  rr <- iroundrobin(it, it2)
  expect_equal(take(rr, 10), as.list(rep(1:5, each=2)))
  expect_error(iterators::nextElem(rr), "StopIteration")
})

test_that("iroundrobin functions properly with two iterators of different type", {
  it3 <- iterators::iter(levels(iris$Species))
  it4 <- iterators::iter(1:3)
  rr <- iroundrobin(it3, it4)
  expect_equal(take(rr, 6), list("setosa", 1, "versicolor", 2, "virginica", 3))
  expect_error(iterators::nextElem(rr), "StopIteration")
})

test_that("iroundrobin functions properly with two iterators of different type - part 2", {
  it5 <- iterators::iter(letters[1:3])
  it6 <- iterators::iter(c(TRUE, FALSE, FALSE))
  rr <- iroundrobin(it5, it6)
  expect_equal(take(rr, 6), list("a", TRUE, "b", FALSE, "c", FALSE))
  expect_error(iterators::nextElem(rr), "StopIteration")
})

test_that("iroundrobin functions properly with multiple iterators of same length", {
  it <- iterators::iter(1:5)
  it2 <- iterators::iter(1:5)
  it3 <- iterators::iter(1:5)
  rr <- iroundrobin(it, it2, it3)
  expect_equal(take(rr, 15), as.list(rep(1:5, each=3)))
  expect_error(iterators::nextElem(rr), "StopIteration")
})

test_that("iroundrobin works functions properly with multiple iterators of unequal length", {
  it <- iterators::iter(letters[1:3])
  it2 <- iterators::iter(LETTERS[5])
  it3 <- iterators::iter(letters[25:26])
  rr <- iroundrobin(it, it2, it3)
  expect_equal(take(rr, 6), list("a", "E", "y", "b", "z", "c"))
  expect_error(iterators::nextElem(rr), "StopIteration")
})

test_that("iroundrobin works functions properly with multiple vectors of unequal length", {
  rr <- iroundrobin(1:3, 5, 25:26)
  expect_equal(take(rr, 6), list(1, 5, 25, 2, 26, 3))
  expect_error(iterators::nextElem(rr), "StopIteration")
})
