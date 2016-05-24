context("insert")

test_that("insert", {
  # list        
  xs1 = list(a=1, b=2)
  expect_equal(insert(xs1, list(a=99, c=5)), list(a=99, b=2, c=5))
  expect_equal(insert(xs1, list(a=list(99), c=5)), list(a=list(99), b=2, c=5))
  # vector 
  xs1 = c(a=1, b=2) 
  expect_equal(insert(xs1, c(a=99, c=5)), c(a=99, b=2, c=5))
  expect_equal(insert(xs1, c()), xs1)
})

