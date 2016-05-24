context("Evaluation methods")
set.seed(1234)
parts <- create_holdout_partition(toyml)
result <- predict(br(parts$train, "SVM"), parts$test)

test_that("Multi-label confusion matrix", {
  mlconfmat <- multilabel_confusion_matrix(parts$test, result)

  expected <- parts$test$dataset[,parts$test$labels$index]
  predicted <- as.bipartition(result)

  expect_equal(dim(mlconfmat$Z), dim(result))
  expect_equal(dim(mlconfmat$Y), dim(result))
  expect_equal(dim(mlconfmat$R), dim(result))
  expect_equal(mlconfmat$Z, predicted)
  expect_equal(mlconfmat$Y, expected)

  expect_true(all(apply(mlconfmat$R, 1, function (row) row %in% 1:5)))

  expect_equal(dim(mlconfmat$TP), dim(result))
  expect_equal(dim(mlconfmat$TN), dim(result))
  expect_equal(dim(mlconfmat$FP), dim(result))
  expect_equal(dim(mlconfmat$FN), dim(result))

  expect_equal(length(mlconfmat$Zi), parts$test$measures$num.instances)
  expect_equal(length(mlconfmat$Yi), parts$test$measures$num.instances)
  expect_equal(length(mlconfmat$Zl), parts$test$measures$num.labels)
  expect_equal(length(mlconfmat$Yl), parts$test$measures$num.labels)
  expect_equal(mlconfmat$Yi, apply(expected, 1, sum))
  expect_equal(mlconfmat$Zi, apply(predicted, 1, sum))
  expect_equal(mlconfmat$Yl, apply(expected, 2, sum))
  expect_equal(mlconfmat$Zl, apply(predicted, 2, sum))

  totals <- mlconfmat$TPi + mlconfmat$TNi + mlconfmat$FPi + mlconfmat$FNi
  expect_true(all(totals == 5))

  expect_equal(mlconfmat$TPi, apply(expected & predicted, 1, sum))
  expect_equal(mlconfmat$TNi, apply(!expected & !predicted, 1, sum))
  expect_equal(mlconfmat$FPi, apply(!expected & predicted, 1, sum))
  expect_equal(mlconfmat$FNi, apply(expected & !predicted, 1, sum))

  expect_equal(mlconfmat$TPl, apply(expected & predicted, 2, sum))
  expect_equal(mlconfmat$TNl, apply(!expected & !predicted, 2, sum))
  expect_equal(mlconfmat$FPl, apply(!expected & predicted, 2, sum))
  expect_equal(mlconfmat$FNl, apply(expected & !predicted, 2, sum))

  expect_error(multilabel_confusion_matrix(parts$train, result))
})

test_that("Bipartition measures", {
  labels <- as.matrix(parts$test$dataset[, parts$test$labels$index])
  expected <- parts$test$dataset[, parts$test$labels$index]

  #100% correct
  test.result <- multilabel_prediction(labels, labels, TRUE)
  mlconfmat <- multilabel_confusion_matrix(parts$test, test.result)
  expect_equal(utiml_measure_accuracy(mlconfmat), 1)
  expect_equal(utiml_measure_f1(mlconfmat), 1)
  expect_equal(utiml_measure_subset_accuracy(mlconfmat), 1)
  expect_equal(utiml_measure_precision(mlconfmat), 1)
  expect_equal(utiml_measure_recall(mlconfmat), 1)
  expect_equal(utiml_measure_hamming_loss(mlconfmat), 0)

  expect_equal(utiml_measure_macro_AUC(mlconfmat), 1)
  expect_equal(utiml_measure_macro_precision(mlconfmat), 1)
  expect_equal(utiml_measure_micro_precision(mlconfmat), 1)
  expect_equal(utiml_measure_macro_recall(mlconfmat), 1)
  expect_equal(utiml_measure_micro_AUC(mlconfmat), 1)
  expect_equal(utiml_measure_micro_recall(mlconfmat), 1)
  expect_equal(utiml_measure_macro_f1(mlconfmat), 1)
  expect_equal(utiml_measure_micro_f1(mlconfmat), 1)

  #100% incorrect
  for (i in seq(ncol(labels))) {
    pos <- labels[, i] == 1
    neg <- !pos
    labels[pos, i] <- 0
    labels[neg, i] <- 1
  }
  test.result <- multilabel_prediction(labels, labels, TRUE)
  mlconfmat <- multilabel_confusion_matrix(parts$test, test.result)
  expect_equal(utiml_measure_accuracy(mlconfmat), 0)
  expect_equal(utiml_measure_f1(mlconfmat), 0)
  expect_equal(utiml_measure_subset_accuracy(mlconfmat), 0)
  expect_equal(utiml_measure_precision(mlconfmat), 0)
  expect_equal(utiml_measure_recall(mlconfmat), 0)
  expect_equal(utiml_measure_hamming_loss(mlconfmat), 1)

  expect_equal(utiml_measure_macro_AUC(mlconfmat), 0)
  expect_equal(utiml_measure_macro_precision(mlconfmat), 0)
  expect_equal(utiml_measure_micro_precision(mlconfmat), 0)
  expect_equal(utiml_measure_macro_recall(mlconfmat), 0)
  expect_equal(utiml_measure_micro_AUC(mlconfmat), 0)
  expect_equal(utiml_measure_micro_recall(mlconfmat), 0)
  expect_equal(utiml_measure_macro_f1(mlconfmat), 0)
  expect_equal(utiml_measure_micro_f1(mlconfmat), 0)

  #Random
  set.seed(1234)
  for (i in seq(ncol(labels))) {
    labels[, i] <- utiml_normalize(rnorm(nrow(labels)))
  }
  labels <- fixed_threshold(labels, 0.5)
  test.result <- multilabel_prediction(labels, labels, TRUE)
  mlconfmat <- multilabel_confusion_matrix(parts$test, test.result)
  measures <- list(
    Accuracy = mean(rowSums(expected & labels) / rowSums(expected | labels)),
    FMeasure = mean(2 * rowSums(expected & labels) /
                      (rowSums(expected) + rowSums(labels))),
    SubsetAccuracy = mean(rowSums(expected == labels) == ncol(labels)),
    Precision = mean(rowSums(expected & labels) / rowSums(labels)),
    Recall = mean(rowSums(expected & labels) / rowSums(expected)),
    HammingLoss = mean(unlist(lapply(seq(nrow(labels)), function (i) {
      sum(xor(labels[i,], expected[i,])) / ncol(labels)
    }))),
    MacroPrecision = mean(
      colSums(labels == 1 & expected == 1) / colSums(labels == 1)
    ),
    MicroPrecision = sum(colSums(labels == 1 & expected == 1)) /
      sum(colSums(labels == 1)),
    MacroRecall = mean(
      colSums(labels == 1 & expected == 1) / colSums(expected == 1)
    ),
    MicroRecall = sum(colSums(labels == 1 & expected == 1)) /
      sum(colSums(expected == 1)),
    MacroFMeasure = (function (){
      prec <- colSums(labels == 1 & expected == 1) / colSums(labels == 1)
      rec <- colSums(labels == 1 & expected == 1) / colSums(expected == 1)
      mean(2 * prec * rec / (prec + rec))
    })(),
    MicroFMeasure = (function (){
      prec <- sum(colSums(labels == 1 & expected == 1)) /
        sum(colSums(labels == 1))
      rec <- sum(colSums(labels == 1 & expected == 1)) /
        sum(colSums(expected == 1))
      2 * prec * rec / (prec + rec)
    })()
  )
  expect_equal(utiml_measure_accuracy(mlconfmat), measures$Accuracy)
  expect_equal(utiml_measure_f1(mlconfmat), measures$FMeasure)
  expect_equal(utiml_measure_subset_accuracy(mlconfmat),
               measures$SubsetAccuracy)
  expect_equal(utiml_measure_precision(mlconfmat), measures$Precision)
  expect_equal(utiml_measure_recall(mlconfmat), measures$Recall)
  expect_equal(utiml_measure_hamming_loss(mlconfmat), measures$HammingLoss)

  expect_equal(utiml_measure_macro_precision(mlconfmat),
               measures$MacroPrecision)
  expect_equal(utiml_measure_micro_precision(mlconfmat),
               measures$MicroPrecision)
  expect_equal(utiml_measure_macro_recall(mlconfmat), measures$MacroRecall)
  expect_equal(utiml_measure_micro_recall(mlconfmat), measures$MicroRecall)
  expect_equal(utiml_measure_macro_f1(mlconfmat), measures$MacroFMeasure)
  expect_equal(utiml_measure_micro_f1(mlconfmat), measures$MicroFMeasure)
})

test_that("Ranking measures", {
  labels <- as.matrix(parts$test$dataset[, parts$test$labels$index])
  expected <- parts$test$dataset[, parts$test$labels$index]

  #100% correct
  test.result <- multilabel_prediction(labels, labels, TRUE)
  mlconfmat <- multilabel_confusion_matrix(parts$test, test.result)

  expect_equal(utiml_measure_one_error(mlconfmat), 0)
  expect_equal(utiml_measure_coverage(mlconfmat),
               parts$test$measures$cardinality - 1)
  expect_equal(utiml_measure_ranking_loss(mlconfmat), 0)
  expect_equal(utiml_measure_average_precision(mlconfmat), 1)
  expect_equal(utiml_measure_margin_loss(mlconfmat), 0)
  expect_equal(utiml_measure_is_error(mlconfmat, mlconfmat$R), 0)

  #100% incorrect
  for (i in seq(ncol(labels))) {
    pos <- labels[, i] == 1
    neg <- !pos
    labels[pos, i] <- 0
    labels[neg, i] <- 1
  }
  test.result <- multilabel_prediction(labels, labels, TRUE)
  mlconfmat <- multilabel_confusion_matrix(parts$test, test.result)

  expect_equal(utiml_measure_one_error(mlconfmat), 1)
  expect_equal(utiml_measure_coverage(mlconfmat), 4)
  expect_equal(utiml_measure_ranking_loss(mlconfmat), 1)
  #TODO study how to determine the worst case
  average.precision <- mean(sapply(seq(nrow(labels)), function (row) {
    Y <- mlconfmat$R[row, expected[row, ] == 1]
    sum(sapply(Y, function (y){
      sum(Y <= y) / y
    })) / length(Y)
  }))
  expect_equal(utiml_measure_average_precision(mlconfmat), average.precision)
  expect_equal(utiml_measure_margin_loss(mlconfmat), 4)
  dif.rank <- mlconfmat$R[, 5:1]
  colnames(dif.rank) <- colnames(mlconfmat$R)
  expect_equal(utiml_measure_is_error(mlconfmat, dif.rank), 1)

  #Random
  set.seed(1234)
  for (i in seq(ncol(labels))) {
    labels[, i] <- utiml_normalize(rnorm(nrow(labels)))
  }
  ranking <- labels
  idxRk <- t(apply(labels, 1, order, decreasing = TRUE))
  for (row in seq(nrow(labels))) {
    ranking[row, idxRk[row,]] <- 1:5
  }
  bipartition <- fixed_threshold(labels, 0.5)
  test.result <- multilabel_prediction(bipartition, labels, TRUE)
  mlconfmat <- multilabel_confusion_matrix(parts$test, test.result)

  measures <- list(
    OneError = mean(sapply(seq(nrow(labels)), function (row) {
        ifelse(expected[row, which.min(ranking[row, ])], 0, 1)
      })),
    Coverage = mean(sapply(seq(nrow(labels)), function (row) {
        max(ranking[row, expected[row,] == 1]) - 1
      })),
    RankingLoss = mean(sapply(seq(nrow(labels)), function (row) {
        idxY <- expected[row, ] == 1
        Y1 <- ranking[row, idxY]
        Y2 <- ranking[row, !idxY]
        sum(sapply(Y1, function (y) sum(y > Y2))) / (length(Y1) * length(Y2))
      })),
    AvgPrecision = mean(sapply(seq(nrow(labels)), function (row) {
      Y <- ranking[row, expected[row, ] == 1]
      sum(sapply(Y, function (y){
        sum(Y <= y) / y
      })) / length(Y)
    })),
    MarginLoss = mean(sapply(seq(nrow(labels)), function (row) {
      idxY <- expected[row, ] == 1
      max(0, max(ranking[row, idxY]) - min(ranking[row, !idxY]))
    })),
    IsError = mean(sapply(seq(nrow(labels)), function (row) {
      rowSums(t(apply(ranking - dif.rank, 1, abs))) > 0
    }))
  )

  expect_equal(utiml_measure_one_error(mlconfmat), measures$OneError)
  expect_equal(utiml_measure_coverage(mlconfmat), measures$Coverage)
  expect_equal(utiml_measure_ranking_loss(mlconfmat), measures$RankingLoss)
  expect_equal(utiml_measure_average_precision(mlconfmat),
               measures$AvgPrecision)
  expect_equal(utiml_measure_margin_loss(mlconfmat), measures$MarginLoss)
  expect_equal(utiml_measure_is_error(mlconfmat, dif.rank), measures$IsError)
})

test_that("Measures names", {
  expect_equal(utiml_measure_names("abc"), c("abc"))

  rankings <- c("average-precision", "coverage", "margin-loss", "one-error",
  "ranking-loss")

  expect_equal(utiml_measure_names("ranking"), rankings)
  expect_equal(utiml_measure_names(c("ranking", "ranking")), rankings)
  expect_equal(utiml_measure_names(c("ranking", "xyz")),
               c(rankings, "xyz"))

  macro <- c("macro-AUC", "macro-F1", "macro-precision", "macro-recall")
  micro <- c("micro-AUC", "micro-F1", "micro-precision", "micro-recall")
  expect_equal(utiml_measure_names("macro-based"), macro)
  expect_equal(utiml_measure_names("micro-based"), micro)
  expect_equal(utiml_measure_names("label-based"), sort(c(micro, macro)))

  example <- sort(c("accuracy", "F1", "hamming-loss", "precision", "recall",
               "subset-accuracy"))
  expect_equal(utiml_measure_names("example-based"), example)

  expect_equal(utiml_measure_names(),
               sort(c(rankings, example, macro, micro)))
})

test_that("Evaluate", {
  expect_equal(length(multilabel_evaluate(parts$test, result, "accuracy")), 1)
  measures <- multilabel_evaluate(parts$test, result, "example-based")
  expect_equal(length(measures), 6)
  expect_true(all(measures >= 0 & measures <= 1))
  expect_named(measures, utiml_measure_names("example-based"))

  mlconfmat <- multilabel_confusion_matrix(parts$test, result)
  expect_equal(measures, multilabel_evaluate(mlconfmat, "example-based"))

  measures <- multilabel_evaluate(mlconfmat, c("hamming-loss", "macro-accuracy",
                                               "micro-accuracy"))
  expect_equal(length(measures), 3)
  expect_true(measures["macro-accuracy"] == measures["micro-accuracy"])
  expect_true(measures["macro-accuracy"] + measures["hamming-loss"] == 1)

  expect_error(multilabel_evaluate(parts$test))
  expect_error(multilabel_evaluate(parts$test, as.matrix(result)))
  expect_error(multilabel_evaluate(parts$test, result, "mymeasure"))
})

test_that("Mulan Measures", {
  expected <- read.csv("../testfiles/flags-expected.csv")
  bipartition <- read.csv("../testfiles/flags-bipartition.csv")
  probability <- read.csv("../testfiles/flags-scores.csv")
  measures <- read.csv("../testfiles/flags-measures.csv")
  mnames <- measures[, 1]
  measures <- measures[, 2]
  names(measures) <- mnames

  dataset <- cbind(attr1=rep(1, nrow(expected)), expected)
  indexes <- seq(ncol(expected)) + 1
  mdata <- mldr::mldr_from_dataframe(dataset, indexes, name="flags")
  mlresult <- multilabel_prediction(bipartition, probability, TRUE)

  evaluation <- multilabel_evaluate(mdata, mlresult)
  expect_equal(evaluation["accuracy"], measures["Accuracy"],
               check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["hamming-loss"], measures["Hamming-Loss"],
               check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["F1"], measures["F-Measure"],
               check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["precision"], measures["Precision"],
               check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["recall"], measures["Recall"],
               check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["subset-accuracy"], measures["Subset-Accuracy"],
               check.attributes = FALSE, tolerance = 1e-3)

  #expect_equal(evaluation["macro-AUC"], measures["Macro-AUC"],
  #             check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["macro-recall"], measures["Macro-Recall"],
               check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["micro-F1"], measures["Micro-F-Measure"],
               check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["macro-F1"], measures["Macro-F-Measure"],
               check.attributes = FALSE, tolerance = 1e-3)
  #expect_equal(evaluation["micro-AUC"], measures["Micro-AUC"],
  #             check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["micro-precision"], measures["Micro-Precision"],
               check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["macro-precision"], measures["Macro-Precision"],
               check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["micro-recall"], measures["Micro-Recall"],
               check.attributes = FALSE, tolerance = 1e-3)


  expect_equal(evaluation["average-precision"], measures["Average-Precision"],
               check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["coverage"], measures["Coverage"],
               check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["one-error"], measures["OneError"],
               check.attributes = FALSE, tolerance = 1e-3)
  expect_equal(evaluation["ranking-loss"], measures["Ranking-Loss"],
               check.attributes = FALSE, tolerance = 1e-2)
  #ranking-loss has a lower precision because there are differences in the
  #ranking order when there are ties

  #expect_equal(evaluation["is-error"], measures["IsError"],
  #             check.attributes = FALSE, tolerance = 1e-3)
  #ranking_error
})

test_that("Sum mlconfmat", {
  mlconfmat <- multilabel_confusion_matrix(parts$test, result)
  dml <- mlconfmat + mlconfmat
  ndim <- dim(mlconfmat$Z)
  ndim[1] <- ndim[1] * 2

  expect_equal(dim(dml$Z), ndim)
  expect_equal(dim(dml$Y), ndim)
  expect_equal(dim(dml$R), ndim)
  expect_equal(dim(dml$TP), ndim)
  expect_equal(dim(dml$FP), ndim)
  expect_equal(dim(dml$TN), ndim)
  expect_equal(dim(dml$FN), ndim)
  expect_equal(length(dml$Zi), ndim[1])
  expect_equal(length(dml$Yi), ndim[1])
  expect_equal(dml$Zl, mlconfmat$Zl * 2)
  expect_equal(dml$Yl, mlconfmat$Yl * 2)
  expect_equal(length(dml$TPi), ndim[1])
  expect_equal(length(dml$FPi), ndim[1])
  expect_equal(length(dml$TNi), ndim[1])
  expect_equal(length(dml$FNi), ndim[1])
  expect_equal(dml$TPl, mlconfmat$TPl * 2)
  expect_equal(dml$FPl, mlconfmat$FPl * 2)
  expect_equal(dml$TNl, mlconfmat$TNl * 2)
  expect_equal(dml$FNl, mlconfmat$FNl * 2)

  ip1 <- seq(1, 10)
  ip2 <- seq(11, parts$test$measures$num.instances)
  part1 <- create_subset(parts$test, ip1)
  part2 <- create_subset(parts$test, ip2)
  dml <- multilabel_confusion_matrix(part1, result[ip1, ]) +
    multilabel_confusion_matrix(part2, result[ip2, ])
  expect_equal(dml, mlconfmat)

  expected <- read.csv("../testfiles/flags-expected.csv")
  bipartition <- read.csv("../testfiles/flags-bipartition.csv")
  probability <- read.csv("../testfiles/flags-scores.csv")
  dataset <- cbind(attr1=rep(1, nrow(expected)), expected)
  mdata <- mldr::mldr_from_dataframe(dataset, seq(ncol(expected)) + 1,
                                     name="flags")
  mlresult <- multilabel_prediction(bipartition, probability, TRUE)
  mlconfmat2 <- multilabel_confusion_matrix(mdata, mlresult)
  expect_error(mlconfmat + mlconfmat2)
})
