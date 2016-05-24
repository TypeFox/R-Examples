test_that("tuneModel tunes", {
    ps = makeParamSet(makeIntegerParam("test", lower = 1, upper = 100))
    design = generateRandomDesign(10, ps)
    res = tuneModel(d, classify, testclassifier, design, misclassificationPenalties, nfolds = 2L, quiet = TRUE)

    expect_equal(class(res), "llama.model")
    expect_equal(attr(res, "type"), "classify")

    expect_equal(dim(res$predictions), c(20, 4))
    expect_true(res$parvals$test >= 1 && res$parvals$test <= 100)
    expect_true(all(sapply(res$inner.parvals, function(x) (x$test >= 1 && x$test <= 100))))
})
