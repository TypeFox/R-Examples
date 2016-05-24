context("rowLapply / rowSapply")

test_that("rowLapply", {
  df = data.frame(a = 1:10, b = 10:1)
  expect_true(all(rowLapply(df, length, unlist = TRUE) == 2))
  expect_true(all(rowLapply(df, sum, unlist = TRUE) == 11))
  expect_true(all(unlist(rowLapply(df, Negate(is.list), unlist = TRUE))))
  expect_true(all(unlist(rowLapply(df, is.list))))
  fun = function(x, y) sum(c(unlist(x), y))
  expect_equal(rowLapply(df[, 1L, drop = FALSE], fun, y = 1), as.list(2:11))
})


test_that("rowSapply", {
  df = data.frame(a = 1:10, b = 10:1)
  rownames(df) = letters[1:10]
  y1 = rep(2, nrow(df))
  y2 = setNames(y1, rownames(df))
  expect_equal(rowSapply(df, length, simplify = TRUE, use.names = FALSE), y1)
  expect_equal(rowSapply(df, length, simplify = TRUE, use.names = TRUE), y2)
  expect_equal(rowSapply(df, length, simplify = FALSE, use.names = FALSE), as.list(y1))
  expect_equal(rowSapply(df, length, simplify = FALSE, use.names = TRUE), as.list(y2))
  x1 = rowSapply(df, unlist, simplify = TRUE, use.names = TRUE)
  x2 = sapply(1:nrow(df), function(i) unlist(df[i,]), simplify = TRUE, USE.NAMES = FALSE)
  rownames(x2) = NULL; colnames(x2) = rownames(df)
  expect_equal(x1, x2)
  x1 = rowSapply(df, unlist, simplify = "rows", use.names = FALSE)
  x2 = as.matrix(data.frame(a = 1:10, b = 10:1))
  expect_equal(x1, x2)
  x1 = rowSapply(data.frame(a = 1:2), function(r) r$a, simplify = "rows", use.names = FALSE)
  x2 = matrix(1:2, nrow = 2)
  expect_equal(x1, x2)
})
