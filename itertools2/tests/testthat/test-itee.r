context("itee iterator")

test_that("itee returns n independent numeric vectors", {
  # Iterate through each of the iterators without any order in mind
  iter_list <- itee(1:5, n=3)
  expect_equal(nextElem(iter_list[[1]]), 1)
  expect_equal(nextElem(iter_list[[1]]), 2)
  expect_equal(nextElem(iter_list[[1]]), 3)

  expect_equal(nextElem(iter_list[[2]]), 1)
  expect_equal(nextElem(iter_list[[2]]), 2)
  
  expect_equal(nextElem(iter_list[[3]]), 1)
  expect_equal(nextElem(iter_list[[3]]), 2)

  expect_equal(nextElem(iter_list[[1]]), 4)
  expect_equal(nextElem(iter_list[[1]]), 5)
  expect_error(nextElem(iter_list[[1]]), "StopIteration")

  expect_equal(nextElem(iter_list[[2]]), 3)
  expect_equal(nextElem(iter_list[[2]]), 4)
  expect_equal(nextElem(iter_list[[2]]), 5)
  expect_error(nextElem(iter_list[[2]]), "StopIteration")

  expect_equal(nextElem(iter_list[[3]]), 3)
  expect_equal(nextElem(iter_list[[3]]), 4)
  expect_equal(nextElem(iter_list[[3]]), 5)
  expect_error(nextElem(iter_list[[3]]), "StopIteration")

  # After the iterators are exhausted, ensure that they are truly exhausted
  expect_error(nextElem(iter_list[[1]]), "StopIteration")
  expect_error(nextElem(iter_list[[2]]), "StopIteration")
  expect_error(nextElem(iter_list[[3]]), "StopIteration")
})

# Based on GitHub Issue #36
test_that("itee returns n independent numeric vectors based on n iterators", {
  # Iterate through each of the iterators without any order in mind
  it <- iterators::iter(1:5)
  iter_list <- itee(it, n=3)
  expect_equal(nextElem(iter_list[[1]]), 1)
  expect_equal(nextElem(iter_list[[1]]), 2)
  expect_equal(nextElem(iter_list[[1]]), 3)

  expect_equal(nextElem(iter_list[[2]]), 1)
  expect_equal(nextElem(iter_list[[2]]), 2)
  
  expect_equal(nextElem(iter_list[[3]]), 1)
  expect_equal(nextElem(iter_list[[3]]), 2)

  expect_equal(nextElem(iter_list[[1]]), 4)
  expect_equal(nextElem(iter_list[[1]]), 5)
  expect_error(nextElem(iter_list[[1]]), "StopIteration")

  expect_equal(nextElem(iter_list[[2]]), 3)
  expect_equal(nextElem(iter_list[[2]]), 4)
  expect_equal(nextElem(iter_list[[2]]), 5)
  expect_error(nextElem(iter_list[[2]]), "StopIteration")

  expect_equal(nextElem(iter_list[[3]]), 3)
  expect_equal(nextElem(iter_list[[3]]), 4)
  expect_equal(nextElem(iter_list[[3]]), 5)
  expect_error(nextElem(iter_list[[3]]), "StopIteration")

  # After the iterators are exhausted, ensure that they are truly exhausted
  expect_error(nextElem(iter_list[[1]]), "StopIteration")
  expect_error(nextElem(iter_list[[2]]), "StopIteration")
  expect_error(nextElem(iter_list[[3]]), "StopIteration")
})
