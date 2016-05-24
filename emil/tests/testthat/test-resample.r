context("Resampling schemes")

y <- y_na <- factor(rep(1:2, c(10, 40)))
y_na[c(12, 14, 16)] <- NA

test_that("resampling general", {
    for(method in c("bootstrap", "crossvalidation", "holdout")){
        scheme <- resample(method, y)
        expect_that(scheme, is_a("resample"))
        expect_that(scheme, is_a(method))

        sub_scheme <- subresample(scheme[[1]], y)
        expect_identical(attributes(scheme[[1]]),
                         attributes(sub_scheme[[1]]))

        expect_true(all(is.na(sub_scheme[index_test(scheme[[1]]),])))
        expect_true(all(!is.na(sub_scheme[index_fit(scheme[[1]]),])))

        scheme <- resample(method, y_na)
        expect_true(all(is.na(scheme[is.na(y_na),])))
    }
})

test_that("repeated holdout", {
    ho <- resample("holdout", y, test_fraction=1/5, nfold=5)

    ho.tab <- lapply(ho, table, y)
    expect_true(all(sapply(ho.tab[-1], all.equal, ho.tab[[1]])))
})

test_that("cross-validation", {
    cv <- resample("crossvalidation", y, nfold=5, nrepeat=3)

    cv.tab <- lapply(cv, table, y)
    expect_true(all(sapply(cv.tab[-1], all.equal, cv.tab[[1]])))

    y[1] <- NA
    cv <- resample("crossvalidation", y, nfold=5, nrepeat=3)
    cv.tab <- lapply(cv, table, y)
    expect_that(range(sapply(cv.tab, "[", 1)), is_equivalent_to(1:2))
    expect_that(range(sapply(cv.tab, "[", 2)), is_equivalent_to(7:8))
    expect_true(all(sapply(cv.tab[-1], function(x) all.equal(x[,2], cv.tab[[1]][,2]))))
})

