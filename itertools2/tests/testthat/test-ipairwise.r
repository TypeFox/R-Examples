context("ipairwise iterator")

test_that("ipairwise functions properly with an iterator", {
  it <- iterators::iter(letters[1:4])
  it_ipairwise <- ipairwise(it)
  expect_equal(iterators::nextElem(it_ipairwise), list("a", "b"))
  expect_equal(iterators::nextElem(it_ipairwise), list("b", "c"))
  expect_equal(iterators::nextElem(it_ipairwise), list("c", "d"))
  expect_error(iterators::nextElem(it_ipairwise), "StopIteration")
})

test_that("ipairwise functions properly with a vector", {
  it_ipairwise <- ipairwise(letters[1:5])
  expect_equal(iterators::nextElem(it_ipairwise), list("a", "b"))
  expect_equal(iterators::nextElem(it_ipairwise), list("b", "c"))
  expect_equal(iterators::nextElem(it_ipairwise), list("c", "d"))
  expect_equal(iterators::nextElem(it_ipairwise), list("d", "e"))
  expect_error(iterators::nextElem(it_ipairwise), "StopIteration")
})
