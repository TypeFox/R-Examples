test_that("cvFolds splits", {
    d = list(data=data.frame(a=rep.int(0, 10)), best=rep.int(0, 10))
    class(d) = "llama.data"

    dfs = cvFolds(d)
    expect_equal(length(dfs$train), 10)
    expect_equal(length(dfs$test), 10)
    expect_true(all(sapply(dfs$train, length) == 9))
    expect_true(all(sapply(dfs$test, length) == 1))
})

test_that("cvFolds complains when there's not enough data", {
    d = list(data=data.frame(a=rep.int(0, 5)), best=rep.int(0, 5))
    class(d) = "llama.data"

    expect_error(cvFolds(d), "Requested 10 folds, but cannot produce this many.")
})

test_that("cvFolds splits with best list", {
    d = list(data=data.frame(a=rep.int(0, 10)))
    d$best = list(0, 0, 0, c(0,1), 0, 0, c(0,1), 0, 0, 0)
    class(d) = "llama.data"

    dfs = cvFolds(d, nfolds=2L)
    expect_equal(length(dfs$train), 2)
    expect_equal(length(dfs$test), 2)
    expect_true(all(sapply(dfs$train, length) == 5))
    expect_true(all(sapply(dfs$test, length) == 5))
})

test_that("cvFolds allows to specify number of folds", {
    d = list(data=data.frame(a=rep.int(0, 10)), best=rep.int(0, 10))
    class(d) = "llama.data"

    dfs = cvFolds(d, nfolds=2L)
    expect_equal(length(dfs$train), 2)
    expect_equal(length(dfs$test), 2)
    expect_true(all(sapply(dfs$train, length) == 5))
    expect_true(all(sapply(dfs$test, length) == 5))
})

test_that("cvFolds allows -1 for leave-one-out", {
    d = list(data=data.frame(a=rep.int(0, 10)), best=rep.int(0, 10))
    class(d) = "llama.data"

    dfs = cvFolds(d, nfolds=-1L)
    expect_equal(length(dfs$train), 10)
    expect_equal(length(dfs$test), 10)
    expect_true(all(sapply(dfs$train, length) == 9))
    expect_true(all(sapply(dfs$test, length) == 1))
})

test_that("cvFolds stratifies", {
    d = list(data=data.frame(a=rep.int(0, 10)), best=c(rep.int(0, 5), rep.int(1, 5)))
    class(d) = "llama.data"

    dfs = cvFolds(d, nfolds=5L, T)
    expect_equal(length(dfs$train), 5)
    expect_equal(length(dfs$test), 5)
    expect_true(all(sapply(dfs$train, length) == 8))
    expect_true(all(sapply(dfs$test, length) == 2))

    expect_true(all(sapply(dfs$train, function(x) { sum(d$best[x]==0) }) == 4))
    expect_true(all(sapply(dfs$train, function(x) { sum(d$best[x]==1) }) == 4))
    expect_true(all(sapply(dfs$test, function(x) { sum(d$best[x]==0) }) == 1))
    expect_true(all(sapply(dfs$test, function(x) { sum(d$best[x]==1) }) == 1))
})

test_that("cvFolds replaces existing splits", {
    d = list(data=data.frame(a=rep.int(0, 10)), best=rep.int(0, 10), train=1, test=2)
    class(d) = "llama.data"

    dfs = cvFolds(d)
    expect_equal(length(dfs$train), 10)
    expect_equal(length(dfs$test), 10)
    expect_true(all(sapply(dfs$train, length) == 9))
    expect_true(all(sapply(dfs$test, length) == 1))
})

test_that("cvFolds class and attributes", {
    d = list(data=data.frame(a=rep.int(0, 10)), best=rep.int(0, 10))
    class(d) = "llama.data"

    dfs = cvFolds(d)
    expect_equal(class(dfs), "llama.data")
    expect_true(attr(dfs, "hasSplits"))
})
