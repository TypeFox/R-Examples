
context("asMatrix")

test_that("asMatrix", {
  # empty
  expect_equal(
    asMatrixCols(list()),
    matrix(0, nrow = 0L, ncol = 0L)
  )
  expect_equal(
    asMatrixRows(list()),
    matrix(0, nrow = 0L, ncol = 0L)
  )

  # normal
  expect_equal(
    asMatrixCols(list(c(1, 2), c(3, 3), c(4, 4))),
    matrix(c(1, 2, 3, 3, 4, 4), nrow = 2, ncol = 3, byrow = FALSE)
  )
  expect_equal(
    asMatrixRows(list(c(1, 2), c(3, 3), c(4, 4))),
    matrix(c(1, 2, 3, 3, 4, 4), nrow = 3, ncol = 2, byrow = TRUE)
  )

  # names
  expect_equal(
    asMatrixCols(list(a = c(1, 2), b = c(3, 3), c = c(4, 4))),
    setColNames(matrix(c(1, 2, 3, 3, 4, 4), nrow = 2, ncol = 3, byrow = FALSE), c("a", "b", "c"))
  )
  expect_equal(
    asMatrixRows(list(a = c(1, 2), b = c(3, 3), c = c(4, 4))),
    setRowNames(matrix(c(1, 2, 3, 3, 4, 4), nrow = 3, ncol = 2, byrow = TRUE), c("a", "b", "c"))
  )
  expect_equal(
    asMatrixRows(list(a = c(x = 1, y = 2), b = c(3, 3), c = c(4, 4))),
    setColNames(setRowNames(matrix(c(1, 2, 3, 3, 4, 4), nrow = 3, ncol = 2, byrow = TRUE),
        c("a", "b", "c")), c("x", "y"))
  )
  # manually define rownames
  expect_equal(
    asMatrixCols(list(a = c(1, 2), b = c(3, 3), c = c(4, 4)), row.names = c("xx", "yy")),
    setColNames(
      setRowNames(matrix(c(1, 2, 3, 3, 4, 4), nrow = 2, ncol = 3, byrow = FALSE), c("xx", "yy")),
    c("a", "b", "c"))
  )
  # manually define rownames, but use ints
  expect_equal(
    asMatrixCols(list(a = c(1, 2), b = c(3, 3), c = c(4, 4)), row.names = 1:2),
    setColNames(
      setRowNames(matrix(c(1, 2, 3, 3, 4, 4), nrow = 2, ncol = 3, byrow = FALSE), 1:2),
    c("a", "b", "c"))
  )
  # manually define colnames
  expect_equal(
    asMatrixCols(list(a = c(1, 2), b = c(3, 3), c = c(4, 4)), col.names = c("xx", "yy", "zz")),
    setColNames(
      matrix(c(1, 2, 3, 3, 4, 4), nrow = 2, ncol = 3, byrow = FALSE),
    c("xx", "yy", "zz"))
  )
  expect_equal(
    asMatrixRows(list(a = c(1, 2), b = c(3, 3), c = c(4, 4)), col.names = c("xx", "yy")),
    setRowNames(
      setColNames(
        matrix(c(1, 2, 3, 3, 4, 4), nrow = 3, ncol = 2, byrow = TRUE),
        c("xx", "yy")
      ), 
      c("a", "b", "c")
    )
  )
  # manually define colnames, but use ints
  expect_equal(
    asMatrixCols(list(a = c(1, 2), b = c(3, 3), c = c(4, 4)), col.names = 1:3),
    setColNames(
      matrix(c(1, 2, 3, 3, 4, 4), nrow = 2, ncol = 3, byrow = FALSE),
    1:3)
  )
})

