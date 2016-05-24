context("type checking")

test_that("type checking", {
  expect_null({
    check(1L,2L,is.integer)
    NULL
  })
  expect_null({
    check(c(1L,2L),c(2L,3L),is.integer, length = 2)
    NULL
  })
  expect_error({
    check(1L,1.5,is.integer)
  })
  expect_error({
    check(c(1,1.5),c(2,3.5),is.integer, length = 3)
    NULL
  })
})
