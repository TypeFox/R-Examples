context("Base learner")
testdataset <- list(
  data = toyml$dataset[, c(1,2,3,4,13)],
  labelindex = 5,
  labelname = "y3",
  mlmethod = "test",
  mldataset = "toyml"
)
testdataset$data$y3 <- as.factor(testdataset$data$y3)

test_that("test train/prediction base learner methods", {
  methods <- c("baseSVM", "baseJ48", "baseC5.0", "baseCART", "baseRF",
               "baseNB", "baseKNN", "baseMAJORITY", "baseRANDOM")
  names(methods) <- c("svm", "J48", "C5.0", "rpart", "randomForest",
                      "naiveBayes", "baseKNN", "majorityModel", "randomModel")
  for (modelname in names(methods)) {
    class(testdataset) <- methods[modelname]
    model <- mltrain(testdataset)
    expect_is(model, modelname)
    result <- mlpredict(model, testdataset$data[11:20, 1:4])
    expect_is(result, "data.frame", label = modelname)
    expect_equal(rownames(result), as.character(11:20), label = modelname)
    expect_equal(colnames(result), c("prediction", "probability"),
                 label = modelname)
  }
  class(testdataset) <- NULL
  expect_error(mltrain(testdataset))
  expect_error(mlpredict(testdataset))
})

test_that("KNN", {
  #TODO Train with a k and predict with other k
})

test_that("Summary", {
  class(testdataset) <- "mltransformation"
  expect_equal(summary(testdataset), summary(testdataset$data))
})


