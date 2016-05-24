context("Ensemble tests")

test_that("Majority votes", {
  probs <- matrix(
    c(1, 1, 1, 1, 0.6, 0.1, 0.8, 0.2, 0.8, 0.3, 0.4, 0.1),
    ncol = 3
  )
  preds <- matrix(
    unlist(as.numeric(probs > 0.5)),
    ncol = 3
  )

  result <- utiml_ensemble_majority(preds, probs)
  expect_equal(result$bipartition, c(1,0,1,0))
  expect_equal(result$probability, c(0.8,0.2,0.9,0.15))

  result2 <- utiml_ensemble_majority(preds, probs)
  expect_equal(result, result2)
})

test_that("Other votes", {
  probs <- matrix(
    c(1, 1, 1, 1, 0.6, 0.1, 0.8, 0.2, 0.8, 0.3, 0.4, 0.1),
    ncol = 3
  )
  preds <- matrix(
    unlist(as.numeric(probs > 0.5)),
    ncol = 3
  )

  result <- utiml_ensemble_maximum(preds, probs)
  expect_equivalent(result$bipartition, c(1,1,1,1))
  expect_equivalent(result$probability, c(1,1,1,1))

  result <- utiml_ensemble_minimum(preds, probs)
  expect_equivalent(result$bipartition, c(1,0,0,0))
  expect_equivalent(result$probability, c(0.6,0.1,0.4,0.1))

  result <- utiml_ensemble_average(preds, probs)
  expect_equivalent(result$bipartition, c(1,0,1,0))
  expect_equivalent(result$probability, c(0.8, 1.4/3, 2.2/3, 1.3/3))
})

test_that("Check votes and method", {
  valids <- c(avg  = "utiml_ensemble_average", maj  = "utiml_ensemble_majority",
    max  = "utiml_ensemble_maximum", min  = "utiml_ensemble_minimum")

  for (schema in names(valids)) {
    expect_equal(utiml_ensemble_method(schema), valids[[schema]], schema)
    expect_true(utiml_ensemble_check_voteschema(schema))
  }
  expect_true(utiml_ensemble_check_voteschema(NULL, TRUE))
  expect_true(utiml_ensemble_check_voteschema(NULL))

  expect_equal(utiml_ensemble_method("other_method"), "other_method")
  expect_equal(utiml_ensemble_method("xyz"), "xyz")
  expect_error(utiml_ensemble_check_voteschema("other_method"))
  expect_error(utiml_ensemble_check_voteschema("xyz"))
  expect_error(utiml_ensemble_check_voteschema(NULL, FALSE))
})

test_that("Compute ensemble votes", {
  probs <- matrix(
    c(1, 1, 1, 1, 0.6, 0.1, 0.8, 0.2, 0.8, 0.3, 0.4, 0.1),
    ncol = 3
  )
  preds <- matrix(
    unlist(as.numeric(probs > 0.5)),
    ncol = 3
  )

  vmethod <- utiml_ensemble_method("maj")
  result <- utiml_compute_ensemble(preds, probs, vmethod, c(1,2,3,4))
  expect_is(result, "binary.prediction")
  expect_named(result$bipartition, as.character(c(1,2,3,4)))
  expect_named(result$probability, as.character(c(1,2,3,4)))
  expect_equal(result$bipartition, c(1,0,1,0), check.names = FALSE)
  expect_equal(result$probability, c(0.8,0.2,0.9,0.15), check.names = FALSE)

  result <- utiml_compute_ensemble(preds, probs, vmethod, c(4,3,2,1))
  expect_named(result$bipartition, as.character(c(4,3,2,1)))
  expect_named(result$probability, as.character(c(4,3,2,1)))
})

test_that("Binary ensemble predictions", {
  probs <- matrix(
    c(1, 1, 1, 1, 0.6, 0.1, 0.8, 0.2, 0.8, 0.3, 0.4, 0.1),
    ncol = 3
  )
  preds <- matrix(
    unlist(as.numeric(probs > 0.5)),
    ncol = 3
  )

  bpreds <- lapply(seq(ncol(preds)), function (i){
    bipartition <- preds[, i]
    probability <- probs[, i]
    names(bipartition) <- names(probability) <- 1:4
    utiml:::utiml_binary_prediction(bipartition, probability)
  })
  result <- utiml_predict_binary_ensemble(bpreds, "maj")
  expect_is(result, "binary.prediction")
  expect_named(result$bipartition, as.character(c(1,2,3,4)))
  expect_named(result$probability, as.character(c(1,2,3,4)))
  expect_equal(result$bipartition, c(1,0,1,0), check.names = FALSE)
  expect_equal(result$probability, c(0.8,0.2,0.9,0.15), check.names = FALSE)

  bpreds <- lapply(3, function (i){
    bipartition <- preds[, i]
    probability <- probs[, i]
    names(bipartition) <- names(probability) <- 4:1
    utiml_binary_prediction(bipartition, probability)
  })
  result <- utiml_predict_binary_ensemble(bpreds, "maj")
  expect_named(result$bipartition, as.character(c(4,3,2,1)))
  expect_named(result$probability, as.character(c(4,3,2,1)))
})

test_that("Ensemble predict", {
  pred1 <- utiml_predict(list(
    class1 = utiml_binary_prediction(c(1, 1, 1, 1), c(1, 1, 1, 1)),
    class2 = utiml_binary_prediction(c(1, 0, 1, 0), c(0.6,0.1,0.8,0.2)),
    class3 = utiml_binary_prediction(c(1, 0, 0, 0), c(0.8,0.3,0.4,0.1))
  ), TRUE)
  pred2 <- utiml_predict(list(
    class1 = utiml_binary_prediction(c(1, 1, 1, 1), c(1, 1, 1, 1)),
    class2 = utiml_binary_prediction(c(1, 0, 1, 0), c(0.6,0.1,0.8,0.2)),
    class3 = utiml_binary_prediction(c(1, 0, 0, 0), c(0.8,0.3,0.4,0.1))
  ), TRUE)
  pred3 <- utiml_predict(list(
    class1 = utiml_binary_prediction(c(1, 1, 1, 1), c(0.5,0.5,0.5,0.5)),
    class2 = utiml_binary_prediction(c(1, 1, 1, 1), c(0.6,0.6,0.6,0.6)),
    class3 = utiml_binary_prediction(c(1, 1, 1, 1), c(0.7,0.7,0.7,0.7))
  ), TRUE)
  preds <- list(pred1, pred2, pred3)

  result1 <- compute_multilabel_predictions(preds, "maj")
  result2 <- compute_multilabel_predictions(preds, "maj", FALSE)
  expect_is(result1, "mlresult")
  expect_is(result2, "mlresult")
  expect_true(all(result1 == attr(result2, "probs")))
  expect_true(all(result2 == attr(result1, "classes")))

  expect_equal(result2[,1], c(1,1,1,1))
  expect_equal(result2[,2], c(1,0,1,0))
  expect_equal(result2[,3], c(1,0,0,0))
  expect_equal(utiml_predict_ensemble(preds, NULL), preds)

  rownames(pred1) <- rownames(pred2) <- rownames(pred3) <- c(11:14)
  preds <- list(pred1, pred2, pred3)
  result <- compute_multilabel_predictions(preds, "max")
  expect_equal(dimnames(result), dimnames(pred1))
  expected <- matrix(c(1, 0.6, 0.8, 1, 0.6, 0.7, 1, 0.8, 0.7, 1, 0.6, 0.7),
                     ncol = 3, byrow = T)
  dimnames(expected) <- dimnames(pred1)
  expect_equal(result[,1:3], expected)
  expect_equal(utiml_predict_ensemble(preds, NULL), preds)

  rownames(pred3) <- c(1:4)
  expect_error(compute_multilabel_predictions(c(preds, pred3), "maj"))
})
