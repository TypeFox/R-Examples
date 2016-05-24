context("declare")

test_that("declare", {
  expect_true({
    declare(x,y)
    all(vapply(list(x,y),is.null,logical(1L)))
  })
  expect_identical({
    declare(x,y=numeric(),z=1:3)
    list(x,y,z)
  },list(NULL,numeric(),1:3))
})
