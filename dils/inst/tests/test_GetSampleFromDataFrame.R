context("GetSampleFromDataFrame")

test_that("GetSampleFromDataFrame catches bad input", {
  data(iris)  # provides test data set
  
  warn.n.1 <- c(5,6)
  warn.msg.1 <- "n has multiple values; using only the first value"
  bad.n.1 <- "n"
  bad.n.2 <- 0
  bad.n.3 <- -7
  bad.n.4 <- 1.5
  bad.n.5 <- 200
  
  bad.matrix <- as.matrix(iris[,1:4])
  bad.data <- iris[numeric(0),]
  
  expect_warning(GetSampleFromDataFrame(warn.n.1, iris), warn.msg.1)
  expect_error(GetSampleFromDataFrame(bad.n.1, iris))
  expect_error(GetSampleFromDataFrame(bad.n.2, iris))
  expect_error(GetSampleFromDataFrame(bad.n.3, iris))
  expect_error(GetSampleFromDataFrame(bad.n.4, iris))
  expect_error(GetSampleFromDataFrame(bad.n.5, iris))
  expect_error(GetSampleFromDataFrame(10, bad.matrix))
  expect_error(GetSampleFromDataFrame(10, bad.data))
})
