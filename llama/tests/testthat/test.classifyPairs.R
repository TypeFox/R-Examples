test_that("classifyPairs classifies", {
    res = classifyPairs(classifier=idtestclassifier, d)
    expect_equal(unique(res$predictions$id), 11:20)
    by(res$predictions, res$predictions$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("b", "c")))
        expect_equal(ss$score, c(1, 0))
    })
})

test_that("classifyPairs returns predictor", {
    res = classifyPairs(classifier=idtestclassifier, d)
    fold$id = 1:10
    preds = res$predictor(fold)
    expect_equal(unique(preds$id), 1:10)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("b", "c")))
        expect_equal(ss$score, c(1, 0))
    })
})

test_that("classifyPairs returns predictor that works without IDs", {
    res = classifyPairs(classifier=idtestclassifier, d)
    fold$id = 1:10
    preds = res$predictor(fold[d$features])
    expect_equal(unique(preds$id), 1:10)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("b", "c")))
        expect_equal(ss$score, c(1, 0))
    })
})

test_that("classifyPairs raises error without classifier", {
    expect_error(classifyPairs())
})

test_that("classifyPairs raises error without data", {
    expect_error(classifyPairs(testclassifier))
})

test_that("classifyPairs raises error without train/test split", {
    expect_error(classifyPairs(testclassifier, dnosplit), "Need data with train/test split!")
})

test_that("classifyPairs respects minimize", {
    res = classifyPairs(classifier=idtestclassifier, d)
    expect_equal(unique(res$predictions$id), 11:20)
    by(res$predictions, res$predictions$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("b", "c")))
        expect_equal(ss$score, c(1, 0))
    })

    fold$id = 1:10
    preds = res$predictor(fold)
    expect_equal(unique(preds$id), 1:10)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("b", "c")))
        expect_equal(ss$score, c(1, 0))
    })
})

test_that("classifyPairs allows combination function", {
    res = classifyPairs(classifier=idtestclassifier, e, combine=othertestclassifier)
    expect_equal(unique(res$predictions$id), 11:20)
    by(res$predictions, res$predictions$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("a", "b")))
        expect_equal(ss$score, c(1, 0))
    })

    fold$id = 1:10
    preds = res$predictor(fold)
    expect_equal(unique(preds$id), 1:10)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("a", "b")))
        expect_equal(ss$score, c(1, 0))
    })
})

test_that("classifyPairs works with NA predictions", {
    res = classifyPairs(classifier=natestclassifier, d)
    expect_equal(unique(res$predictions$id), 11:20)
    by(res$predictions, res$predictions$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("b", "c")))
        expect_equal(ss$score, as.numeric(c(NA, NA)))
    })
    fold$id = 1:10
    preds = res$predictor(fold)
    expect_equal(unique(preds$id), 1:10)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("b", "c")))
        expect_equal(ss$score, as.numeric(c(NA, NA)))
    })

    res = classifyPairs(classifier=natestclassifier, combine=natestclassifier, d)
    expect_equal(unique(res$predictions$id), 11:20)
    by(res$predictions, res$predictions$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("b", "c")))
        expect_equal(ss$score, as.numeric(c(NA, NA)))
    })
})

test_that("classifyPairs works with one class train data", {
    # when run with --as-cran, this fails because the llama package that provides the classifier isn't installed
    skip.expensive()
    res = classifyPairs(classifier=makeLearner("classif.rpart"), one)
    expect_equal(unique(res$predictions$id), 11:20)
    by(res$predictions, res$predictions$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("a")))
        expect_equal(ss$score, c(1))
    })
})

test_that("classifyPairs works with probabilities", {
    res = classifyPairs(classifier=probtestclassifier, d)
    expect_equal(unique(res$predictions$id), 11:20)
    by(res$predictions, res$predictions$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("c", "b")))
        expect_equal(ss$score, c(.2, .1))
    })
    fold$id = 1:10
    preds = res$predictor(fold)
    expect_equal(unique(preds$id), 1:10)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("c", "b")))
        expect_equal(ss$score, c(.2, .1))
    })

    res = classifyPairs(classifier=probtestclassifier, combine=probtestclassifier, d)
    expect_equal(unique(res$predictions$id), 11:20)
    by(res$predictions, res$predictions$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("c", "b")))
        expect_equal(ss$score, c(.2, .1))
    })
    fold$id = 1:10
    preds = res$predictor(fold)
    expect_equal(unique(preds$id), 1:10)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("c", "b")))
        expect_equal(ss$score, c(.2, .1))
    })
})
