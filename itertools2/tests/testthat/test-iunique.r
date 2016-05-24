context("iunique iterators")

test_that("iunique works with numeric vectors", {
  x <- rep(1:5, each=10)
  it_unique <- iunique(x)
  expect_equal(take(it_unique, 5), as.list(1:5))
  expect_error(iterators::nextElem(it_unique), "StopIteration")
})

test_that("iunique works with character vectors", {
  x <- as.character(gl(5, 10, labels=LETTERS[1:5]))
  it_unique <- iunique(x)
  expect_equal(take(it_unique, 5), as.list(LETTERS[1:5]))
  expect_error(iterators::nextElem(it_unique), "StopIteration")
})

test_that("iunique works with iterators from numeric vectors", {
  x <- rep(1:5, each=10)
  it <- iterators::iter(rep(x, 2))
  it_unique <- iunique(it)
  expect_equal(take(it_unique, 5), as.list(1:5))
  expect_error(iterators::nextElem(it_unique), "StopIteration")
})

test_that("iunique works with iterators from character vectors", {
  x <- as.character(gl(5, 10, labels=LETTERS[1:5]))
  it <- iterators::iter(rep(x, 2))
  it_unique <- iunique(it)
  expect_equal(take(it_unique, 5), as.list(LETTERS[1:5]))
  expect_error(iterators::nextElem(it_unique), "StopIteration")
})

test_that("iunique_justseen works with numeric vectors", {
  x <- rep(1:5, each=10)
  it_unique <- iunique_justseen(x)
  expect_equal(take(it_unique, 5), as.list(1:5))
  expect_error(iterators::nextElem(it_unique), "StopIteration")
})

test_that("iunique_justseen works with character vectors", {
  x <- as.character(gl(5, 10, labels=LETTERS[1:5]))
  it_unique <- iunique_justseen(x)
  expect_equal(take(it_unique, 5), as.list(LETTERS[1:5]))
  expect_error(iterators::nextElem(it_unique), "StopIteration")
})

test_that("iunique_justseen works with iterators from numeric vectors", {
  num_reps <- 7
  x <- rep(1:5, each=10)
  it <- iterators::iter(rep(x, num_reps))
  it_unique <- iunique_justseen(it)
  expect_equal(take(it_unique, 5 * num_reps), as.list(rep(1:5, num_reps)))
  expect_error(iterators::nextElem(it_unique), "StopIteration")
})

test_that("iunique_justseen works with iterators from character vectors", {
  num_reps <- 7
  x <- as.character(gl(5, 10, labels=LETTERS[1:5]))
  it <- iterators::iter(rep(x, num_reps))
  it_unique <- iunique_justseen(it)
  expect_equal(take(it_unique, 5 * num_reps), as.list(rep(LETTERS[1:5], num_reps)))
  expect_error(iterators::nextElem(it_unique), "StopIteration")
})
