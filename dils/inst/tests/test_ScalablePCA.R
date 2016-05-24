context("ScalablePCA")

test_that("ScalablePCA catches bad input", {
  data(iris)  # provides test data set
  bad.x.matrix <- as.matrix(iris[,1:4])
  
  # ADD TESTS FOR filename
  # ADD TESTS FOR db
  
  good.subsample <- 5
  warn.subsample.1 <- c(10, 11)
  warn.subsample.msg.1 <- "subsample has multiple values; using only the first value"
  bad.subsample.msg <- "subsample must be a positive whole number or FALSE"
  bad.subsample.1 <- TRUE
  bad.subsample.2 <- 0
  bad.subsample.3 <- -5
  bad.subsample.4 <- 10.5
  bad.subsample.5 <- 200
  bad.subsample.6 <- "clam"
  
  good.n.subsamples <- 5
  warn.n.subsamples.1 <- c(25, 30)
  warn.n.subsamples.msg.1 <- "n.subsamples has multiple values; using only the first value"
  bad.n.subsamples.msg <- "n.subsamples must be a positive whole number"
  bad.n.subsamples.1 <- "clam"
  bad.n.subsamples.2 <- 0
  bad.n.subsamples.3 <- -5
  
  good.ignore.cols.1 <- 5
  bad.ignore.cols.msg <- "ignore.cols must be numeric with positive whole values"
  bad.ignore.cols.1 <- c(0, 5)
  bad.ignore.cols.2 <- c(-1, 5)
  bad.ignore.cols.3 <- c(1.5, 5)
  bad.ignore.cols.4 <- "clam"

  good.use.cols <- 1:4
  bad.use.cols.msg <- "use.cols must be numeric with positive whole values"
  bad.use.cols.1 <- c(0, 1:4)
  bad.use.cols.2 <- c(5, 1:4)
  bad.use.cols.3 <- c(-1, 1:4)
  bad.use.cols.4 <- c(1.5, 1:4)
  bad.use.cols.5 <- "clam"
  
  expect_error(ScalablePCA(x=bad.x.matrix, subsample=good.subsample, 
                           n.subsamples=good.n.subsamples, 
                           use.cols=good.use.cols), 
               "x must be a data.frame")
  
  expect_warning(ScalablePCA(x=iris, subsample=warn.subsample.1,
                             n.subsamples=good.n.subsamples, 
                             use.cols=good.use.cols),
                 warn.subsample.msg.1)
  expect_error(ScalablePCA(x=iris, subsample=bad.subsample.1,
                           n.subsamples=good.n.subsamples, 
                           use.cols=good.use.cols),
               bad.subsample.msg)
  expect_error(ScalablePCA(x=iris, subsample=bad.subsample.2,
                           n.subsamples=good.n.subsamples, 
                           use.cols=good.use.cols),
               bad.subsample.msg)
  expect_error(ScalablePCA(x=iris, subsample=bad.subsample.3,
                           n.subsamples=good.n.subsamples, 
                           use.cols=good.use.cols),
               bad.subsample.msg)
  expect_error(ScalablePCA(x=iris, subsample=bad.subsample.4,
                           n.subsamples=good.n.subsamples, 
                           use.cols=good.use.cols),
               bad.subsample.msg)
  expect_error(ScalablePCA(x=iris, subsample=bad.subsample.5,
                           n.subsamples=good.n.subsamples, 
                           use.cols=good.use.cols),
               "nrow")
  expect_error(ScalablePCA(x=iris, subsample=bad.subsample.6,
                           n.subsamples=good.n.subsamples, 
                           use.cols=good.use.cols),
               bad.subsample.msg)
  
  expect_warning(ScalablePCA(x=iris, subsample=good.subsample,
                             n.subsamples=warn.n.subsamples.1, 
                             use.cols=good.use.cols),
                 warn.n.subsamples.msg.1)
  expect_error(ScalablePCA(x=iris, subsample=good.subsample,
                           n.subsamples=bad.n.subsamples.1, 
                           use.cols=good.use.cols),
               bad.n.subsamples.msg)
  expect_error(ScalablePCA(x=iris, subsample=good.subsample,
                           n.subsamples=bad.n.subsamples.2, 
                           use.cols=good.use.cols),
               bad.n.subsamples.msg)
  expect_error(ScalablePCA(x=iris, subsample=good.subsample,
                           n.subsamples=bad.n.subsamples.3, 
                           use.cols=good.use.cols),
               bad.n.subsamples.msg)

  expect_error(ScalablePCA(x=iris, subsample=good.subsample,
                           n.subsamples=good.n.subsamples, 
                           ignore.cols=bad.ignore.cols.1),
               bad.ignore.cols.msg)
  expect_error(ScalablePCA(x=iris, subsample=good.subsample,
                           n.subsamples=good.n.subsamples, 
                           ignore.cols=bad.ignore.cols.2),
               bad.ignore.cols.msg)
  expect_error(ScalablePCA(x=iris, subsample=good.subsample,
                           n.subsamples=good.n.subsamples, 
                           ignore.cols=bad.ignore.cols.3),
               bad.ignore.cols.msg)
  expect_error(ScalablePCA(x=iris, subsample=good.subsample,
                           n.subsamples=good.n.subsamples, 
                           ignore.cols=bad.ignore.cols.4),
               bad.ignore.cols.msg)
  
  expect_error(ScalablePCA(x=iris, subsample=good.subsample,
                           n.subsamples=good.n.subsamples, 
                           use.cols=bad.use.cols.1),
               bad.use.cols.msg)
  expect_error(ScalablePCA(x=iris, subsample=good.subsample,
                           n.subsamples=good.n.subsamples, 
                           use.cols=bad.use.cols.2),
               "infinite or missing values in 'x'")
  expect_error(ScalablePCA(x=iris, subsample=good.subsample,
                           n.subsamples=good.n.subsamples, 
                           use.cols=bad.use.cols.3),
               bad.use.cols.msg)
  expect_error(ScalablePCA(x=iris, subsample=good.subsample,
                           n.subsamples=good.n.subsamples, 
                           use.cols=bad.use.cols.4),
               bad.use.cols.msg)
  expect_error(ScalablePCA(x=iris, subsample=good.subsample,
                           n.subsamples=good.n.subsamples, 
                           use.cols=bad.use.cols.5),
               bad.use.cols.msg)
})

test_that("ScalablePCA returns correct values", {
  data(iris)
  correct.results <- abs(prcomp(iris[,1:4], center=FALSE, scale.=FALSE)$rotation[,1])
  est.results <- ScalablePCA(iris, subsample=10, use.cols=1:4)
  relative.errors <- abs((correct.results - est.results) / correct.results)
  
  expect_true(all(relative.errors < .05))
})
