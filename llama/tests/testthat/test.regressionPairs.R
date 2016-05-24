test_that("regressionPairs predicts", {
    res = regressionPairs(regressor=testregressor, d)
    expect_equal(unique(res$predictions$id), 11:20)
    by(res$predictions, res$predictions$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("c", "b")))
        expect_equal(ss$score, c(1, -1))
    })
})

test_that("regressionPairs returns predictor", {
    res = regressionPairs(regressor=testregressor, d)
    fold$id = 1:10
    preds = res$predictor(fold)
    expect_equal(unique(preds$id), 1:10)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("c", "b")))
        expect_equal(ss$score, c(1, -1))
    })
})

test_that("regressionPairs returns predictor that works without IDs", {
    res = regressionPairs(regressor=testregressor, d)
    fold$id = 1:10
    preds = res$predictor(fold[d$features])
    expect_equal(unique(preds$id), 1:10)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("c", "b")))
        expect_equal(ss$score, c(1, -1))
    })
})

test_that("regressionPairs raises error without regressor", {
    expect_error(regressionPairs())
})

test_that("regressionPairs raises error without data", {
    expect_error(regressionPairs(testregressor))
})

test_that("regressionPairs raises error without train/test split", {
    fold = data.frame(a=rep.int(0, 10), b=c(rep.int(1, 5), rep.int(0, 5)), c=c(rep.int(0, 5), rep.int(1, 5)), best=c(rep.int("c", 5), rep.int("b", 5)))
    d = list(data=rbind(fold, fold), features=c("a"), minimize=T, performance=c("b", "c"))
    class(d) = "llama.data"
    expect_error(regressionPairs(testregressor, d))
})

test_that("regressionPairs respects minimize", {
    fold = data.frame(a=rep.int(0, 10), best=rep.int("b", 10), foo=rep.int(2, 10), bar=rep.int(1, 10))
    d = list(data=rbind(cbind(fold, id=1:10), cbind(fold, id=11:20)),
        train=list(1:nrow(fold)), test=list(1:nrow(fold) + nrow(fold)),
        features=c("a"), performance=c("foo", "bar"), minimize=F, ids=c("id"))
    class(d) = "llama.data"
    attr(d, "hasSplits") = TRUE
    res = regressionPairs(regressor=testregressor, d)
    expect_equal(unique(res$predictions$id), 11:20)
    by(res$predictions, res$predictions$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("foo", "bar")))
        expect_equal(ss$score, c(1, -1))
    })

    fold$id = 1:10
    preds = res$predictor(fold)
    expect_equal(unique(preds$id), 1:10)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, factor(c("foo", "bar")))
        expect_equal(ss$score, c(1, -1))
    })
})

test_that("regressionPairs allows combine classifier", {
    res = regressionPairs(regressor=testregressor, d, combine=idtestclassifier)
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

test_that("regression works with NA predictions", {
    res = regressionPairs(natestregressor, d)
    expect_equal(unique(res$predictions$id), 11:20)
    by(res$predictions, res$predictions$id, function(ss) {
        expect_equal(ss$algorithm, NA)
        expect_equal(ss$score, Inf)
    })
    fold$id = 1:10
    preds = res$predictor(fold)
    expect_equal(unique(preds$id), 1:10)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, NA)
        expect_equal(ss$score, Inf)
    })

    res = regressionPairs(natestregressor, d, combine=natestclassifier)
    expect_equal(unique(res$predictions$id), 11:20)
    by(res$predictions, res$predictions$id, function(ss) {
        expect_equal(ss$algorithm, NA)
        expect_equal(ss$score, Inf)
    })
    fold$id = 1:10
    preds = res$predictor(fold)
    expect_equal(unique(preds$id), 1:10)
    by(preds, preds$id, function(ss) {
        expect_equal(ss$algorithm, NA)
        expect_equal(ss$score, Inf)
    })
})
