context("Transformation tests")

test_that("binary prediction", {
  set.seed(123)
  probs <- runif(20, 0, 1)
  bipartition <- as.numeric(probs > 0.5)
  result <- utiml_binary_prediction(bipartition, probs)
  expect_is(result, "binary.prediction")
  expect_false(is.null(result$bipartition))
  expect_false(is.null(result$probability))
  expect_equal(result$probability, probs)
  expect_equal(result$bipartition, bipartition)

  names(bipartition) <- names(probs) <- c(21:40)
  result <- utiml_binary_prediction(bipartition, probs)
  expect_equal(names(result$probability), names(probs))
  expect_equal(names(result$bipartition), names(bipartition))
})

test_that("multi-label prediction", {
  set.seed(1)
  probs1 <- runif(20, 0, 1)
  probs2 <- runif(20, 0, 1)
  predictions <- list(
    class1 = utiml_binary_prediction(as.numeric(probs1 > 0.5), probs1),
    class2 = utiml_binary_prediction(as.numeric(probs2 > 0.5), probs2)
  )

  result1 <- utiml_predict(predictions, TRUE)
  expect_null(rownames(result1))
  expect_equal(colnames(result1), c("class1", "class2"))
  expect_equal(result1[,"class1"], predictions$class1$probability)
  expect_equal(result1[,"class2"], predictions$class2$probability)

  result2 <- utiml_predict(predictions, FALSE)
  TP1 <- predictions$class1$bipartition == 1
  TP2 <- predictions$class2$bipartition == 1
  expect_equal(predictions$class1$bipartition[TP1], result2[,"class1"][TP1])
  expect_equal(predictions$class2$bipartition[TP2], result2[,"class2"][TP2])
  expect_true(all(result2[,"class1"][!TP1] | result2[,"class2"][!TP1]))
  expect_true(all(result2[,"class1"][TP2] | result2[,"class2"][TP2]))
  filter <- !TP1 & !TP2
  expect_true(all(result2[,"class1"][filter] != result2[,"class2"][filter]))

  expect_true(all(attr(result1, "classes") == result2))
  expect_true(all(attr(result2, "probs") == result1))
  expect_equal(attr(result1, "type"), "probability")
  expect_equal(attr(result2, "type"), "bipartition")

  bips <- as.numeric(probs1 > 0.5)
  names(probs1) <- names(bips) <- 21:40
  predictions <- lapply(utiml_rename(c('l1', 'l2', 'l3')), function (label){
    utiml_binary_prediction(bips, probs1)
  })
  result <- utiml_predict(predictions, TRUE)
  expect_equal(rownames(result), as.character(21:40))
  expect_equal(result[,"l1"], result[,"l2"])
  expect_equal(result[,"l1"], result[,"l3"])
  result <- utiml_predict(predictions, FALSE)
  expect_equal(rownames(result), as.character(21:40))
})

test_that("prepare data", {
  set.seed(123)
  mydata <- data.frame(
    class1 = runif(10, min = 0, max = 1),
    class2 = factor(as.numeric(runif(10, min = 0, max = 1) > 0.5),
                    levels = c("0", "1"))
  )

  dataset <- utiml_prepare_data(mydata, "testDataset", "mlbase", "test", "SVM")
  expect_is(dataset, "testDataset")
  expect_is(dataset, "baseSVM")
  expect_is(dataset, "mltransformation")

  expect_equal(dataset$data, mydata)
  expect_equal(dataset$labelname, "class2")
  expect_equal(dataset$labelindex, 2)
  expect_equal(dataset$mldataset, "mlbase")
  expect_equal(dataset$mlmethod, "test")

  dataset <- utiml_prepare_data(mydata, "onlytest", "mlbase", "test", "XYZ",
                               extra1="abc", extra2=1:10)
  expect_is(dataset, "onlytest")
  expect_is(dataset, "baseXYZ")
  expect_is(dataset, "mltransformation")

  expect_equal(dataset$data, mydata)
  expect_equal(dataset$labelname, "class2")
  expect_equal(dataset$labelindex, 2)
  expect_equal(dataset$mldataset, "mlbase")
  expect_equal(dataset$mlmethod, "test")
  expect_equal(dataset$extra1, "abc")
  expect_equal(dataset$extra2, 1:10)
})

test_that("create model and predict binary model", {
  set.seed(123)
  mydata <- data.frame(
    attr = runif(10, min = 0, max = 1),
    class2 = factor(as.numeric(runif(10, min = 0, max = 1) > 0.5),
                    levels = c("0", "1")),
    row.names = seq(1, to = 20, by = 2)
  )
  dataset <- utiml_prepare_data(mydata, "testdata", "mlds", "br", "RANDOM")
  model <- utiml_create_model(dataset)
  expect_equal(attr(model, "label"), "class2")
  expect_equal(attr(model, "dataset"), "mlds")

  set.seed(123)
  predict1 <- utiml_predict_binary_model(model, mydata[, 1, drop = FALSE])
  expect_is(predict1, "binary.prediction")
  expect_named(predict1$bipartition, rownames(mydata))
  expect_named(predict1$probability, rownames(mydata))

  model <- utiml_create_model(dataset)
  set.seed(123)
  predict2 <- utiml_predict_binary_model(model, mydata[,1, drop = FALSE])
  expect_equal(predict1$probability, predict2$probability)
  expect_true(all(predict1$probability == predict2$probability))

  predict3 <- utiml_predict_binary_model(model, mydata[,1, drop = FALSE])
  expect_false(all(predict2$probability == predict3$probability))
})

test_that("create binary data", {
  dataset <- utiml_create_binary_data(toyml, "y1")
  expect_equal(ncol(dataset), toyml$measures$num.inputs + 1)
  expect_equal(dataset[seq(toyml$measures$num.inputs)],
               toyml$dataset[toyml$attributesIndexes])
  expect_equal(dataset["y1"], toyml$dataset["y1"])

  dataset <- utiml_create_binary_data(toyml, "y2")
  expect_equal(ncol(dataset), toyml$measures$num.inputs + 1)
  expect_equal(dataset[seq(toyml$measures$num.inputs)],
               toyml$dataset[toyml$attributesIndexes])
  expect_equal(dataset["y2"], toyml$dataset["y2"])

  one.column <- rep(1, toyml$measures$num.instances)
  dataset <- utiml_create_binary_data(toyml, "y3", one.column)
  expect_equal(ncol(dataset), toyml$measures$num.inputs + 2)
  expect_equal(dataset[, length(dataset)-1], one.column)
  expect_equal(dataset["y3"], toyml$dataset["y3"])

  extra.columns <- cbind(a=one.column, b=rnorm(toyml$measures$num.instances))
  dataset <- utiml_create_binary_data(toyml, "y4", extra.columns)
  expect_equal(ncol(dataset), toyml$measures$num.inputs + 3)
  expect_equal(dataset[c("a", "b")], as.data.frame(extra.columns))
  expect_equal(dataset["y4"], toyml$dataset["y4"])
})
