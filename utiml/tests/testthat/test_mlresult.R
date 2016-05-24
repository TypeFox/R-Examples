context("ML Result")

test_that("multilabel_prediction", {
  set.seed(1)
  colrows <- list(as.character(11:20), paste("lbl", 1:10, sep=''))
  probs <- matrix(stats::runif(100), ncol = 10, dimnames = colrows)
  classes <- ifelse(probs >= 0.5, 1, 0)
  mres1 = multilabel_prediction(classes, probs, TRUE)
  mres2 = multilabel_prediction(classes, probs, FALSE)

  expect_is(mres1, "mlresult")
  expect_is(mres2, "mlresult")
  expect_true(is.probability(mres1))
  expect_true(is.bipartition(mres2))
  expect_false(is.probability(mres2))
  expect_false(is.bipartition(mres1))

  expect_equal(as.bipartition(mres1),  classes)
  expect_equal(as.bipartition(mres1),  as.bipartition(mres2))
  expect_equal(as.probability(mres1),  probs)
  expect_equal(as.probability(mres1),  as.probability(mres2))
  expect_equal(dimnames(mres1), colrows)
  expect_equal(dimnames(mres1), dimnames(mres2))

  #TODO test prediction with all labels lower than 0.5
  colrows <- list(as.character(11:20), paste("lbl", 1:5, sep=''))
  probs <- matrix(stats::runif(50, min = 0, max = 0.4), ncol = 5,
                  dimnames = colrows)
  classes <- ifelse(probs >= 0.5, 1, 0)
  mres3 = multilabel_prediction(classes, probs, FALSE)
  expect_true(all(rowSums(mres3) == 1))

  #TODO test prediction with all labels lower than 0.5 and max value repeted
  probs[1:5, c(1,2)] <- 0.45
  probs[6:10, -c(1,2)] <- 0.47
  mres4 = multilabel_prediction(classes, probs, FALSE)
  expect_true(all(rowSums(mres4[1:5, c(1,2)]) == 2))
  expect_true(all(rowSums(mres4[6:10, c(1,2)]) == 0))
  expect_true(all(rowSums(mres4[1:5, -c(1,2)]) == 0))
  expect_true(all(rowSums(mres4[6:10, -c(1,2)]) == 3))
})

test_that("as.ranking", {
  rowcol <- list(as.character(1:5), c("lA", "lB"))
  mlres <- as.mlresult(matrix(seq(0.1, 1, by = 0.1), ncol=2, byrow = TRUE,
                              dimnames = rowcol))
  rk <- as.ranking(mlres)
  exp.rk <- matrix(c(rep(2, 5), rep(1, 5)), ncol = 2, dimnames = rowcol)
  expect_is(rk, "matrix")
  expect_equal(rk, exp.rk)
  expect_equal(dimnames(rk), rowcol)

  mlres2 <- as.mlresult(cbind(mlres, 1 - mlres))
  rk <- as.ranking(mlres2)
  expect_equal(rk[1, ], rk[2, ])
  expect_equal(rk[4, ], rk[5, ])
  expect_equal(rowSums(rk), c(10, 10, 9, 10, 10), check.names = FALSE)

  rk <- as.ranking(mlres2, ties = "average")
  expect_equal(rowSums(rk), c(10, 10, 10, 10, 10), check.names = FALSE)

  rk <- as.ranking(mlres2, ties = "max")
  expect_equal(rowSums(rk), c(10, 10, 11, 10, 10), check.names = FALSE)
})

test_that("as.mlresult", {
  set.seed(1234)
  predictions <- matrix(stats::runif(100), ncol = 10)
  colnames(predictions) <- paste('label', 1:10, sep='')

  mlresult <- as.mlresult(predictions)
  mlresult2 <- as.mlresult(mlresult, prob = FALSE)
  mlresult3 <- as.mlresult(predictions, threshold = 0.3)

  expect_is(mlresult, "mlresult")
  expect_is(mlresult2, "mlresult")
  expect_is(mlresult3, "mlresult")
  expect_equal(as.probability(mlresult), predictions)
  expect_equal(as.mlresult(as.data.frame(predictions)), mlresult)
  expect_equal(as.mlresult(mlresult), mlresult)
  expect_true(attr(mlresult, "type") == attr(mlresult3, "type"))
  expect_false(attr(mlresult, "type") == attr(mlresult2, "type"))

  expect_equal(as.bipartition(mlresult),
               as.bipartition(fixed_threshold(predictions, 0.5)))
  expect_equal(as.bipartition(mlresult3),
               as.bipartition(fixed_threshold(predictions, 0.3)))
})

test_that("Filter ML Result", {
  set.seed(1234)
  labels <- matrix(stats::runif(150), ncol = 10)
  colnames(labels) <- paste("label", 1:10, sep='')
  mlres <- as.mlresult(labels)
  bipartition <- as.bipartition(mlres)

  mlres1 <- multilabel_prediction(bipartition, labels, TRUE)
  mlres2 <- multilabel_prediction(bipartition, labels, FALSE)

  expect_is(mlres1[1:3, ], "mlresult")
  expect_is(mlres1[1:3], "mlresult")
  expect_equal(mlres1[1:3, ], mlres1[1:3])
  expect_is(mlres1[1:3, 1:3], "matrix")
  expect_is(mlres1[, 1:5], "matrix")
  expect_is(mlres1[, 1, drop = FALSE], "matrix")
  expect_is(mlres1[, 1], "numeric")

  expect_true(is.probability(mlres1[1:3, ]))
  expect_true(is.bipartition(mlres2[1:3, ]))

  expect_equal(mlres1[, c("label1", "label3")], labels[, c("label1", "label3")])
  expect_equal(mlres2[, c("label1", "label3")],
               bipartition[, c("label1", "label3")])
  expect_equal(mlres1[, 3:6], labels[, 3:6])
  expect_equal(mlres2[, 2:8], bipartition[, 2:8])
  expect_equal(mlres1[1:5, 1:5], labels[1:5, 1:5])
  expect_equal(mlres2[2:8, 2:8], bipartition[2:8, 2:8])

  mlres3 <- mlres1[1:8]
  expect_equal(as.probability(mlres3), labels[1:8, ])
  expect_equal(as.bipartition(mlres3), bipartition[1:8, ])
})

