test_that("imputeCensored raises error without data", {
    expect_error(imputeCensored())
})

test_that("imputeCensored leaves successes alone", {
    data = data.frame(a=rep.int(0, 10), foo=rep.int(0, 10), bar=rep.int(T, 10))
    d = list(data=data, features=c("a"), performance=c("foo"), success=c("bar"), best=rep.int("b", 10))
    class(d) = "llama.data"
    res = imputeCensored(d, testregressor)
    expect_identical(res$data, res$original_data)
})

test_that("imputeCensored fails with no non-censored values", {
    data = data.frame(a=rep.int(0, 10), foo=rep.int(0, 10), bar=rep.int(F, 10))
    d = list(data=data, features=c("a"), performance=c("foo"), success=c("bar"), best=rep.int("b", 10))
    class(d) = "llama.data"
    expect_error(imputeCensored(d, testregressor), "Cannot impute for  foo , no non-censored values!")
})

test_that("imputeCensored makes everything successes", {
    data = data.frame(a=rep.int(0, 10), foo=rep.int(0, 10), bar=c(rep.int(T, 5), rep.int(F, 5)))
    d = list(data=data, features=c("a"), performance=c("foo"), success=c("bar"), best=rep.int("b", 10))
    class(d) = "llama.data"
    res = imputeCensored(d, testregressor)
    expect_false(all(res$original_data$bar))
    expect_true(all(res$data$bar))
})

test_that("imputeCensored imputes non-successes", {
    data = data.frame(a=rep.int(0, 10), foo=rep.int(0, 10), bar=c(rep.int(T, 5), rep.int(F, 5)))
    d = list(data=data, features=c("a"), performance=c("foo"), success=c("bar"), best=rep.int("b", 10))
    class(d) = "llama.data"
    res = imputeCensored(d, footestregressor)
    expect_identical(res$data$foo[1:5], res$original_data$foo[1:5])
    expect_true(all(res$data$foo[6:10] == 1))
})
