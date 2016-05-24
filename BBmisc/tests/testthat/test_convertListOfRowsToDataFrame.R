context("convertListOfRowstoDataFrame")

test_that("convertListOfRowstoDataFrame", {
  df1 = convertListOfRowsToDataFrame(list(list(x = 1, y = "a"), list(x = 2, y = "b")),
    strings.as.factors = FALSE)
  df2 = data.frame(x = 1:2, y = c("a", "b"), stringsAsFactors = FALSE)
  expect_equal(df1, df2)

  df1 = convertListOfRowsToDataFrame(list(c(x = "1", y = "a"), list(x = "2", y = "b")),
    strings.as.factors = FALSE)
  df2 = data.frame(x=c("1", "2"), y=c("a", "b"), stringsAsFactors = FALSE)
  expect_equal(df1, df2)

  df1 = convertListOfRowsToDataFrame(list(c("1", "a"), c("2", "b")), strings.as.factors = FALSE,
    col.names=c("x", "y"))
  df2 = data.frame(x=c("1", "2"), y=c("a", "b"), stringsAsFactors = FALSE)
  expect_equal(df1, df2)

  df1 = convertListOfRowsToDataFrame(list(list(a = 1, b = 1), list(b = 12)))
  df2 = convertListOfRowsToDataFrame(list(c(a = 1, b = 1), c(b = 12)))
  expect_equal(df1, df2)

  # names
  df1 = convertListOfRowsToDataFrame(list(list(x = 1, y = "a"), list(x = 2, y = "b")),
    strings.as.factors = FALSE, row.names = c("r1", "r2"))
  df2 = setRowNames(data.frame(x = 1:2, y = c("a", "b"), stringsAsFactors = FALSE), c("r1", "r2"))
  expect_equal(df1, df2)

  df1 = convertListOfRowsToDataFrame(list(list(x = 1, y = "a"), list(x = 2, y = "b")),
    strings.as.factors = FALSE, row.names = 1:2)
  df2 = setRowNames(data.frame(x = 1:2, y = c("a", "b"), stringsAsFactors = FALSE), 1:2)
  expect_equal(df1, df2)
  
  df1 = convertListOfRowsToDataFrame(list(list(x = 1, y = "a"), list(x = 2, y = "b")),
    strings.as.factors = FALSE, col.names = c("c1", "c2"))
  df2 = data.frame(c1 = 1:2, c2 = c("a", "b"), stringsAsFactors = FALSE)
  expect_equal(df1, df2)

  df1 = convertListOfRowsToDataFrame(list(list(x = 1, y = "a"), list(x = 2, y = "b")),
    strings.as.factors = FALSE, col.names = 1:2)
  df2 = setColNames(data.frame(1:2, c("a", "b"), stringsAsFactors = FALSE), 1:2)
  expect_equal(df1, df2)
})
