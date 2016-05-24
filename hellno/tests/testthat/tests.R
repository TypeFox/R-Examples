context("Testing: function arguments of base and hellno are the same")

test_that("data.frame()", {
  expect_true(
    all(
      names(as.list(args(hellno::data.frame))) %in%
      names(as.list(args(base::data.frame)))
    )
  )
})

test_that("as.data.frame()", {
  expect_true(
    all(
      names(as.list(args(hellno::data.frame))) %in%
      names(as.list(args(base::data.frame)))
    )
  )
})
