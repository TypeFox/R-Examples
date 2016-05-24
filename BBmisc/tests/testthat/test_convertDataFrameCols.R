context("convertDataFrameCols")

test_that("convertDataFrameCols", {
  df = data.frame(x=1:2, y=c("a", "b"), z=factor(c("x", "y")), stringsAsFactors=FALSE)
  df2 = convertDataFrameCols(df, chars.as.factor=TRUE)
  expect_true(is.numeric(df2$x))
  expect_true(is.factor(df2$y))
  expect_true(is.factor(df2$z))
  expect_equal(df$x, df2$x)
  expect_equal(df$y, as.character(df2$y))
  expect_equal(df$z, df2$z)

  df2 = convertDataFrameCols(df, factors.as.char=TRUE)
  expect_true(is.numeric(df2$x))
  expect_true(is.character(df2$y))
  expect_true(is.character(df2$z))
  expect_equal(df$x, df2$x)
  expect_equal(df$y, df2$y)
  expect_equal(as.character(df$z), df2$z)
  
  df2 = convertDataFrameCols(df, ints.as.num=TRUE)
  expect_true(is.double(df2$x))
  expect_true(is.character(df2$y))
  expect_true(is.factor(df2$z))
  expect_equal(df$x, df2$x)
  expect_equal(df$y, df2$y)
  expect_equal(df$z, df2$z)

  df2 = convertDataFrameCols(df, chars.as.factor=TRUE, factors.as.char=TRUE, ints.as.num=TRUE)
  expect_true(is.double(df2$x))
  expect_true(is.factor(df2$y))
  expect_true(is.character(df2$z))
})


