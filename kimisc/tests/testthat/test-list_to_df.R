context("list_to_df")

test_that("vector input", {
  inp <- 3:1
  ret <- list_to_df(inp)
  expect_true(is.data.frame(ret))
  expect_equal(nrow(ret), 3)
  expect_equal(ret$name, 1:3)
  expect_equal(ret$value, as.list(3:1))
  expect_equal(df_to_list(ret), as.list(inp))
})

test_that("named vector input", {
  inp <- setNames(nm = 3:1)
  ret <- list_to_df(inp)
  expect_true(is.data.frame(ret))
  expect_equal(nrow(ret), 3)
  expect_equal(ret$name, c("3", "2", "1"))
  expect_equal(ret$value, as.list(3:1))
  expect_equal(df_to_list(ret), as.list(inp))
})

test_that("list input", {
  inp <- list(3:1)
  ret <- list_to_df(inp)
  expect_true(is.data.frame(ret))
  expect_equal(nrow(ret), 1)
  expect_equal(ret$name, 1)
  expect_equal(ret$value, list(3:1))
  expect_equal(df_to_list(ret), as.list(inp))
})

test_that("NULL input", {
  inp <- NULL
  ret <- list_to_df(inp)
  expect_true(is.data.frame(ret))
  expect_equal(nrow(ret), 0)
  expect_equal(df_to_list(ret), as.list(inp))
})
