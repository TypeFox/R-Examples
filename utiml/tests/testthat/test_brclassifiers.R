context("BR based classifiers")
train <- toyml
test <- toyml$dataset[10:40, toyml$attributesIndexes]

predictionTest <- function (model) {
  set.seed(123)
  pred <- predict(model, test)
  expect_is(pred, "mlresult")
  expect_equal(nrow(pred), nrow(test))
  expect_equal(ncol(pred), toyml$measures$num.labels)
  expect_equal(colnames(pred), rownames(toyml$labels))
  expect_equal(rownames(pred), rownames(test))

  set.seed(123)
  pred1 <- predict(model, test, prob = FALSE)
  expect_is(pred1, "mlresult")
  expect_equal(as.matrix(pred1), attr(pred, "classes"))
  expect_equal(as.matrix(pred), attr(pred1, "probs"))

  pred
}

baseTest <- function (model, expected.class) {
   expect_is(model, expected.class)
   expect_equal(names(model$models), rownames(toyml$labels))

   predictionTest(model)
}

ensembleTest <- function (model, expected.class) {
  expect_is(model, expected.class)
  expect_equal(length(model$models), model$rounds)

  predictionTest(model)
}

test_that("Binary Relevance", {
  model <- br(train, "RANDOM")
  baseTest(model, "BRmodel")
})

test_that("BR Plus", {
  model <- brplus(train, "RANDOM")
  pred1 <- baseTest(model, "BRPmodel")

  expect_is(model$initial, "BRmodel")
  pred0 <- predict(model$initial, test)

  pred2 <- predict(model, test, strategy="NU")
  expect_equal(colnames(pred2), rownames(train$labels))

  pred3 <- predict(model, test, "Stat")
  expect_equal(colnames(pred3), rownames(train$labels))

  new.chain <- c("y3", "y4", "y1", "y2", "y5")
  pred4 <- predict(model, test, "Ord", new.chain)

  expect_error(predict(model, test, "xay"))
  expect_error(predict(model, test, "Ord"))
  expect_error(predict(model, test, "Ord", new.chain[1:2]))
  expect_error(predict(model, test, "Ord", c(new.chain, "extra")))
  expect_error(predict(model, test, "Ord", c("a", "b", "c")))
})

test_that("Classifier Chain", {
  model <- cc(train, "RANDOM")
  pred <- baseTest(model, "CCmodel")
  mpred <- as.matrix(pred)

  set.seed(123)
  pred1 <- predict(model, test, prob = FALSE)
  expect_is(pred1, "mlresult")
  expect_equal(as.matrix(pred1), attr(pred, "classes"))
  expect_equal(as.matrix(pred), attr(pred1, "probs"))

  new.chain <- c("y5", "y4", "y3", "y2", "y1")
  model2 <- cc(train, "RANDOM", new.chain)
  expect_equal(model2$chain, new.chain)

  set.seed(123)
  pred2 <- predict(model2, test)
  expect_equal(colnames(pred2), rownames(train$labels))

  set.seed(123)
  pred3 <- predict(model2, test)
  expect_false(isTRUE(all.equal(pred3, pred1)))
  expect_equal(pred3, pred2)

  expect_error(cc(train, "RANDOM", chain=c("a", "b", "c", "d", "e")))
  expect_error(cc(train, "RANDOM", chain=c(new.chain, "extra")))
})

test_that("CTRL", {
  model <- ctrl(train, "RANDOM")
  pred1 <- baseTest(model, "CTRLmodel")
  baseTest(ctrl(train, "RANDOM", validation.threshold = 1), "CTRLmodel")

  model2 <- ctrl(train, "RANDOM", m = 2,
                 validation.size = 0.2, validation.threshold = 0)
  pred2 <- baseTest(model2, "CTRLmodel")
  expect_equal(model2$rounds, 2)

  expect_error(ctrl(train, "RANDOM", 0))
  expect_error(ctrl(train, "RANDOM", validation.size=0))
  expect_error(ctrl(train, "RANDOM", validation.size=1))
  expect_error(ctrl(train, "RANDOM", validation.threshold=1.1))
  expect_error(predict(model, test, "ABC"))
  expect_error(predict(model, test, NULL))
})

test_that("EBR", {
  model1 <- ebr(train, "RANDOM")
  pred1 <- ensembleTest(model1, "EBRmodel")

  model2 <- ebr(train, "RANDOM", m=3, subsample=0.5, attr.space=0.40)
  expect_equal(model2$nrow, 50)
  expect_equal(model2$ncol, 4)
  pred2 <- ensembleTest(model2, "EBRmodel")

  predictions <- predict(model1, test, vote.schema = NULL)
  expect_equal(length(predictions), 10)
  predictions <- predict(model2, test, vote.schema = NULL)
  expect_equal(length(predictions), 3)

  expect_error(ebr(train, "RANDOM", subsample=0))
  expect_error(ebr(train, "RANDOM", attr.space=0))
  expect_error(ebr(train, "RANDOM", m=0))
  expect_error(predict(model1, test, "ABC"))
})

test_that("ECC", {
  model1 <- ecc(train, "RANDOM")
  pred1 <- ensembleTest(model1, "ECCmodel")

  model2 <- ecc(train, "RANDOM", m=3, subsample=0.5, attr.space=0.40)
  expect_equal(model2$nrow, 50)
  expect_equal(model2$ncol, 4)
  pred2 <- ensembleTest(model2, "ECCmodel")

  predictions <- predict(model1, test, vote.schema = NULL)
  expect_equal(length(predictions), 10)
  predictions <- predict(model2, test, vote.schema = NULL)
  expect_equal(length(predictions), 3)

  expect_error(ecc(train, "RANDOM", subsample=0))
  expect_error(ecc(train, "RANDOM", attr.space=0))
  expect_error(ecc(train, "RANDOM", m=0))
  expect_error(predict(model1, test, "ABC"))
})

test_that("DBR", {
  model <- dbr(train, "RANDOM")
  pred <- baseTest(model, "DBRmodel")
  expect_is(model$estimation, "BRmodel")

  estimative <- predict(model$estimation, test, prob = FALSE)
  pred1 <- predict(model, test, estimative)
  expect_equal(pred1, pred)

  model <- dbr(train, "RANDOM", estimate = FALSE)
  expect_error(predict(model, test))
})

test_that("Meta-BR", {
  model <- mbr(train, "RANDOM")
  pred <- baseTest(model, "MBRmodel")
  expect_is(model$basemodel, "BRmodel")

  model <- mbr(train, "RANDOM", folds=2, phi=0.3)
  pred <- baseTest(model, "MBRmodel")

  expect_error(mbr(train, "RANDOM", folds=0))
  expect_error(mbr(train, "RANDOM", phi=1.1))
})

test_that("Nestest Stack", {
  model <- ns(train, "RANDOM")
  pred <- baseTest(model, "NSmodel")
  mpred <- as.matrix(pred)

  set.seed(123)
  pred1 <- predict(model, test, prob = FALSE)
  expect_is(pred1, "mlresult")
  expect_equal(as.matrix(pred1), attr(pred, "classes"))
  expect_equal(as.matrix(pred), attr(pred1, "probs"))

  new.chain <- c("y5", "y4", "y3", "y2", "y1")
  model2 <- ns(train, "RANDOM", new.chain)
  expect_equal(model2$chain, new.chain)

  set.seed(123)
  pred2 <- predict(model2, test)
  expect_equal(colnames(pred2), rownames(train$labels))

  set.seed(123)
  pred3 <- predict(model2, test)
  expect_false(isTRUE(all.equal(pred3, pred1)))
  expect_equal(pred3, pred2)

  expect_error(ns(train, "RANDOM", chain=c("a", "b", "c", "d", "e")))
  expect_error(ns(train, "RANDOM", chain=c(new.chain, "extra")))
})

test_that("Prudent", {
  model <- prudent(train, "RANDOM")
  #TODO pred <- baseTest(model, "PruDentmodel")
  expect_is(model$basemodel, "BRmodel")

  model <- prudent(train, "RANDOM", phi=0.3)
  #TODO pred <- baseTest(model, "PruDentmodel")

  expect_error(prudent(train, "RANDOM", phi=1.1))
})

test_that("RDBR", {
  model <- rdbr(train, "RANDOM")
  pred <- baseTest(model, "RDBRmodel")
  expect_is(model$estimation, "BRmodel")

  estimative <- predict(model$estimation, test, prob = FALSE)
  pred1 <- predict(model, test, estimative)
  expect_equal(pred1, pred)

  model <- dbr(train, "RANDOM", estimate = FALSE)
  expect_error(predict(model, test))
})
