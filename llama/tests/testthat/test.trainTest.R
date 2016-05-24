test_that("trainTest splits", {
    d = list(data=data.frame(a=rep.int(0, 10)), best=rep.int(0, 10))
    class(d) = "llama.data"

    dtt = trainTest(d)
    expect_equal(length(dtt$train), 1)
    expect_equal(length(dtt$test), 1)
    expect_equal(length(dtt$train[[1]]), 6)
    expect_equal(length(dtt$test[[1]]), 4)
})

test_that("trainTest splits with best list", {
    d = list(data=data.frame(a=rep.int(0, 10)))
    d$best = list(0, 0, 0, c(0, 1), 0, 0, 0, 0, 0, 0)
    class(d) = "llama.data"

    dtt = trainTest(d)
    expect_equal(length(dtt$train), 1)
    expect_equal(length(dtt$test), 1)
    expect_equal(length(dtt$train[[1]]), 6)
    expect_equal(length(dtt$test[[1]]), 4)
})

test_that("trainTest allows to specify split ratio", {
    d = list(data=data.frame(a=rep.int(0, 10)), best=rep.int(0, 10))
    class(d) = "llama.data"

    dtt = trainTest(d, trainpart=0.1)
    expect_equal(length(dtt$train), 1)
    expect_equal(length(dtt$test), 1)
    expect_equal(length(dtt$train[[1]]), 1)
    expect_equal(length(dtt$test[[1]]), 9)
})

test_that("trainTest stratifies", {
    d = list(data=data.frame(a=rep.int(0, 10)), best=c(rep.int(0, 5), rep.int(1, 5)))
    class(d) = "llama.data"

    dtt = trainTest(d, stratify = T)
    expect_equal(length(dtt$train), 1)
    expect_equal(length(dtt$test), 1)
    expect_equal(length(dtt$train[[1]]), 6)
    expect_equal(length(dtt$test[[1]]), 4)

    expect_equal(sum(d$best[dtt$train[[1]]]==0), 3)
    expect_equal(sum(d$best[dtt$train[[1]]]==1), 3)
    expect_equal(sum(d$best[dtt$test[[1]]]==0), 2)
    expect_equal(sum(d$best[dtt$test[[1]]]==1), 2)
})

test_that("trainTest replaces existing splits", {
    d = list(data=data.frame(a=rep.int(0, 10)), best=rep.int(0, 10), train=c(1,1), test=c(1,1))
    class(d) = "llama.data"

    dtt = trainTest(d)
    expect_equal(length(dtt$train), 1)
    expect_equal(length(dtt$test), 1)
    expect_equal(length(dtt$train[[1]]), 6)
    expect_equal(length(dtt$test[[1]]), 4)
})
