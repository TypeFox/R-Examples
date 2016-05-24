context("Sampling tests")
set.seed(1)
df <- data.frame(matrix(rnorm(1000), ncol = 10))
df$Label1 <- c(sample(c(0,1), 100, replace = TRUE))
df$Label2 <- c(sample(c(0,1), 100, replace = TRUE))
df$Label3 <- c(sample(c(0,1), 100, replace = TRUE))
df$Label4 <- as.numeric(df$Label1 == 0 | df$Label2 == 0 | df$Label3 == 0)
mdata <- mldr::mldr_from_dataframe(df, labelIndices = c(11, 12, 13, 14),
                             name = "testMLDR")
empty.mdata <- mldr::mldr_from_dataframe(df[, 1:13],
                                         labelIndices = c(11, 12, 13),
                                         name = "testMLDR")
set.seed(NULL)

testFolds <- function (kfold, original, msg) {
  real <- unlist(lapply(kfold$fold, names))
  expected <- unique(real)
  expect_true(all(expected == real), label=msg)
  comp <- sort(unlist(kfold$fold)) == seq(original$measures$num.instances)
  expect_true(all(comp), label=msg)
  expect_true(all(sort(real) == sort(rownames(original$dataset))), label=msg)
}

testEmptyIntersectRows <- function (a, b) {
  expect_equal(length(intersect(rownames(a$dataset), rownames(b$dataset))), 0)
}

testCompletude <- function (list, original) {
  names <- sort(unlist(lapply(list, function (fold)  rownames(fold$dataset))))
  expect_true(all(names == sort(rownames(original$dataset))))
}

test_that("random holdout", {
  folds <- create_holdout_partition(mdata, 0.7, "random")
  expect_equal(length(folds), 2)
  expect_is(folds[[1]], "mldr")
  expect_is(folds[[2]], "mldr")
  expect_equal(folds[[1]]$measures$num.instances, 70)
  expect_equal(folds[[2]]$measures$num.instances, 30)
  expect_equal(rownames(folds[[1]]$labels), rownames(folds[[2]]$labels))
  testEmptyIntersectRows(folds[[1]], folds[[2]])
  testCompletude(folds, mdata)

  subfolds <- create_holdout_partition(folds[[1]], 0.5, "random")
  testEmptyIntersectRows(subfolds[[1]], subfolds[[2]])
  testCompletude(subfolds, folds[[1]])

  folds <- create_holdout_partition(mdata, c("train"=0.5, "test"=0.5))
  expect_named(folds, c("train", "test"))
  expect_equal(folds$train$measures$num.instances,
               folds$test$measures$num.instances)
  testCompletude(folds, mdata)

  folds <- create_holdout_partition(empty.mdata, c("train"=0.5, "test"=0.5))
  expect_equal(folds$train$measures$num.instances,
               folds$test$measures$num.instances)
  testCompletude(folds, mdata)

  set.seed(1)
  f1 <- create_holdout_partition(mdata, c(0.5, 0.5))
  set.seed(1)
  f2 <- create_holdout_partition(mdata, c(0.5, 0.5))
  expect_equal(f1, f2)
  set.seed(NULL)

  expect_error(create_holdout_partition(mdata, NULL))
  expect_error(create_holdout_partition(mdata, c(0.5,0.8,0.1)))
})

test_that("stratified holdout", {
  f <- create_holdout_partition(mdata, c("a"=0.4, "b"=0.4, "c"=0.2),
                                "stratified")
  expect_equal(length(f), 3)
  expect_named(f, c("a", "b", "c"))
  expect_equal(f[[1]]$measures$num.instances, 40)
  expect_equal(f[[2]]$measures$num.instances, 40)
  expect_equal(f[[3]]$measures$num.instances, 20)
  expect_equal(rownames(f[[1]]$labels), rownames(f[[2]]$labels))
  expect_equal(rownames(f[[1]]$labels), rownames(f[[3]]$labels))

  testEmptyIntersectRows(f$a, f$b)
  testEmptyIntersectRows(f$a, f$c)
  testEmptyIntersectRows(f$b, f$c)
  testCompletude(f, mdata)

  folds <- create_holdout_partition(empty.mdata, c("train"=0.5, "test"=0.5),
                                    "stratified")
  expect_equal(folds$train$measures$num.instances,
               folds$test$measures$num.instances)
  testCompletude(folds, mdata)

  sf <- create_holdout_partition(f$a, c("a"=0.5, "b"=0.5), "stratified")
  expect_equal(length(sf), 2)
  testEmptyIntersectRows(sf$a, sf$b)
  testCompletude(sf, f$a)
})

test_that("iterative holdout", {
  f <- create_holdout_partition(mdata, c("a"=0.4, "b"=0.4, "c"=0.1, "d"=0.1),
                                "iterative")
  expect_equal(length(f), 4)
  expect_named(f, c("a", "b", "c", "d"))
  expect_equal(rownames(f[[1]]$labels), rownames(f[[2]]$labels))
  expect_equal(rownames(f[[1]]$labels), rownames(f[[3]]$labels))

  testEmptyIntersectRows(f$a, f$b)
  testEmptyIntersectRows(f$a, f$c)
  testEmptyIntersectRows(f$a, f$d)
  testEmptyIntersectRows(f$b, f$c)
  testEmptyIntersectRows(f$b, f$d)
  testEmptyIntersectRows(f$c, f$d)
  testCompletude(f, mdata)

  folds <- create_holdout_partition(empty.mdata, c("train"=0.5, "test"=0.5),
                                    "iterative")
  testEmptyIntersectRows(folds$train, folds$test)
  testCompletude(folds, mdata)

  sf <- create_holdout_partition(f$a, c("a"=0.5, "b"=0.5), "iterative")
  expect_equal(length(sf), 2)
  testEmptyIntersectRows(sf$a, sf$b)
  testCompletude(sf, f$a)

  folds <- create_holdout_partition(mdata, 0.6)
  sf <- create_holdout_partition(folds[[2]], c("a"=0.6, "b"=0.4), "iterative")
  testEmptyIntersectRows(sf$a, sf$b)
  testCompletude(sf, folds[[2]])
})

test_that("random kfold", {
  f <- create_kfold_partition(mdata, 10, "random")
  expect_is(f, "kFoldPartition")
  expect_equal(f$k, 10)
  expect_equal(length(f$fold), 10)
  for (i in 1:10) {
    expect_equal(length(f$fold[[i]]), 10)
  }
  testFolds(f, mdata, "f Random kfolds")

  fdata1 <- partition_fold(f, 1)
  fdata2 <- partition_fold(f, 10)

  expect_equal(rownames(fdata1$labels), rownames(fdata2$labels))
  expect_equal(fdata1$measures$num.instances, fdata2$measures$num.instances)

  set.seed(1)
  f1 <- create_kfold_partition(mdata, 4)
  testFolds(f1, mdata, "f1 Random kfolds")
  set.seed(1)
  f2 <- create_kfold_partition(mdata, 4)
  expect_equal(length(f1$fold), 4)
  expect_equal(length(f1$fold[[2]]), 25)
  expect_equal(f1, f2)
  expect_false(all(f1$fold[[1]] %in% f$fold[[1]]))
  set.seed(NULL)

  f3 <- create_kfold_partition(mdata, 3)
  testFolds(f3, mdata, "f3 Random kfolds")
  expect_equal(f3$k, 3)
  expect_equal(length(unlist(f3$fold)), 100)
  expect_more_than(length(f3$fold[[1]]), 32)
  expect_more_than(length(f3$fold[[2]]), 32)
  expect_more_than(length(f3$fold[[3]]), 32)

  ds <- create_holdout_partition(mdata, c("train" = 0.9, "test" = 0.1))
  f4 <- create_kfold_partition(ds$train, 9)
  testFolds(f4, ds$train, "f4 Random kfolds")

  folds <- create_kfold_partition(empty.mdata, 5)
  testFolds(folds, empty.mdata, "empty Random kfolds")
})

test_that("stratified kfold", {
  f <- create_kfold_partition(mdata, 10, "stratified")
  expect_is(f, "kFoldPartition")
  expect_equal(f$k, 10)
  expect_equal(length(f$fold), 10)
  for (i in 1:10) {
    expect_equal(length(f$fold[[i]]), 10)
  }

  testFolds(f, mdata, "f Stratified kfold")
  fdata1 <- partition_fold(f, 1)
  fdata2 <- partition_fold(f, 10)

  expect_equal(rownames(fdata1$labels), rownames(fdata2$labels))
  expect_equal(fdata1$measures$num.instances, fdata2$measures$num.instances)

  set.seed(1)
  f1 <- create_kfold_partition(mdata, 4, "stratified")
  testFolds(f1, mdata, "f1 Stratified kfold")
  set.seed(1)
  f2 <- create_kfold_partition(mdata, 4, "stratified")
  expect_equal(length(f1$fold), 4)
  expect_equal(length(f1$fold[[2]]), 25)
  expect_equal(f1, f2)
  expect_false(all(f1$fold[[1]] %in% f$fold[[1]]))
  set.seed(NULL)

  f3 <- create_kfold_partition(mdata, 3, "stratified")
  testFolds(f3, mdata, "f3 Stratified kfold")
  expect_equal(f3$k, 3)
  expect_equal(length(unlist(f3$fold)), 100)
  expect_more_than(length(f3$fold[[1]]), 32)
  expect_more_than(length(f3$fold[[2]]), 32)
  expect_more_than(length(f3$fold[[3]]), 32)

  ds <- create_holdout_partition(mdata, c("train" = 0.9, "test" = 0.1))
  f4 <- create_kfold_partition(ds$train, 9, "stratified")
  testFolds(f4, ds$train, "f4 Stratified kfold")

  folds <- create_kfold_partition(empty.mdata, 5, "stratified")
  testFolds(folds, empty.mdata, "empty stratified kfolds")
})

test_that("iterative kfold", {
  f <- create_kfold_partition(mdata, 10, "iterative")
  expect_is(f, "kFoldPartition")
  expect_equal(f$k, 10)
  expect_equal(length(f$fold), 10)
  for (i in 1:10) {
    expect_more_than(length(f$fold[[i]]), 7)
  }

  testFolds(f, mdata, "f Iterative kfold")
  fdata1 <- partition_fold(f, 1)
  fdata2 <- partition_fold(f, 10)

  expect_equal(rownames(fdata1$labels), rownames(fdata2$labels))
  expect_equal(fdata1$measures$num.instances, fdata2$measures$num.instances)

  set.seed(1)
  f1 <- create_kfold_partition(mdata, 4, "iterative")
  testFolds(f1, mdata, "f1 Iterative kfold")
  set.seed(1)
  f2 <- create_kfold_partition(mdata, 4, "iterative")
  expect_equal(length(f1$fold), 4)
  expect_equal(length(f1$fold[[2]]), 25)
  expect_equal(f1, f2)
  expect_false(all(f1$fold[[1]] %in% f$fold[[1]]))
  set.seed(NULL)

  f3 <- create_kfold_partition(mdata, 3, "iterative")
  testFolds(f3, mdata, "f3 Iterative kfold")
  expect_equal(f3$k, 3)
  expect_equal(length(unlist(f3$fold)), 100)
  expect_more_than(length(f3$fold[[1]]), 30)
  expect_more_than(length(f3$fold[[2]]), 30)
  expect_more_than(length(f3$fold[[3]]), 30)

  ds <- create_holdout_partition(mdata, c("train" = 0.9, "test" = 0.1))
  f4 <- create_kfold_partition(ds$train, 9, "iterative")
  testFolds(f4, ds$train, "f4 Iterative kfold")

  folds <- create_kfold_partition(empty.mdata, 5, "iterative")
  testFolds(folds, empty.mdata, "empty iterative kfolds")
})

test_that("subset and random subset", {
  rows <- 10:20
  cols <- 3:7

  data <- create_subset(mdata, rows, 1:10)
  expect_is(data, "mldr")
  expect_equal(data$measures$num.attributes, mdata$measures$num.attributes)
  expect_equal(create_subset(mdata, rows), data)

  data <- create_subset(mdata, 1:100, cols)
  expect_equal(data$measures$num.instances, mdata$measures$num.instances)
  expect_equal(data$dataset[data$labels$index],
               mdata$dataset[mdata$labels$index])

  data1 <- create_subset(mdata, rows, cols)
  data2 <- create_subset(mdata, rows, cols)
  expect_equal(data1, data2)

  data <- create_subset(mdata, seq(200), seq(30))
  expect_equal(data, mdata)

  data <- create_subset(mdata, c(1,2,3,-5,-4,-10), c(1,2,3,-5,-4,-10))
  #TODO test values

  data <- create_random_subset(mdata, 20, 5)
  expect_equal(data$measures$num.instances, 20)
  expect_equal(data$measures$num.attributes, 5 + data$measures$num.labels)
  expect_equal(data$dataset[,data$labels$index],
               mdata$dataset[rownames(data$dataset), mdata$labels$index])
})

test_that("Alternatives dataset for sampling", {
  dataset <- cbind(mdata$dataset[mdata$labels$index],
                   mdata$dataset[mdata$attributesIndexes])
  ndata <- mldr::mldr_from_dataframe(dataset, labelIndices = 1:4,
                                     name = "testMLDR")

  test <- create_holdout_partition(ndata)
  expect_equal(colnames(test[[1]]$dataset), colnames(ndata$dataset))
  expect_equal(colnames(test[[2]]$dataset), colnames(ndata$dataset))

  test <- create_holdout_partition(ndata, method="iterative")
  expect_equal(colnames(test[[1]]$dataset), colnames(ndata$dataset))
  expect_equal(colnames(test[[2]]$dataset), colnames(ndata$dataset))

  test <- create_holdout_partition(ndata, method="stratified")
  expect_equal(colnames(test[[1]]$dataset), colnames(ndata$dataset))
  expect_equal(colnames(test[[2]]$dataset), colnames(ndata$dataset))

  kf <- create_kfold_partition(ndata, 3)
  test <- partition_fold(kf, 1)
  expect_equal(colnames(test[[1]]$dataset), colnames(ndata$dataset))
  expect_equal(colnames(test[[2]]$dataset), colnames(ndata$dataset))

  test <- partition_fold(kf, 2, has.validation = T)
  expect_equal(colnames(test[[1]]$dataset), colnames(ndata$dataset))
  expect_equal(colnames(test[[2]]$dataset), colnames(ndata$dataset))
  expect_equal(colnames(test[[3]]$dataset), colnames(ndata$dataset))
})

#TODO test not complete partitions (using 90% per example)
