context("listify vectors")

test_that("incorrect inputs raise errors", {
  bads <- list(
    list(1:4, 4:7),
    matrix(1:12, 4, 3)
  )
  
  expect_that(listify(bads[[1]]), throws_error())
  expect_that(listify(bads[[2]]), throws_error("A vector is required"))
})

test_that("listify works as expected", {
  num_ind = c(3, 5, 2)
  num_list = list(c(1,1,1), c(2,2,2,2,2), c(3,3))
  
  expect_that(listify(num_ind), is_equivalent_to(num_list))
})
