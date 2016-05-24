context("Variable summaries")

with_mock_HTTP({
    ds <- loadDataset("test ds")
    gen <- ds$gender

    test_that("table 'method' dispatch", {
        expect_identical(table(1:5), base::table(1:5))
        expect_identical(table(useNA="ifany", 1:5),
            base::table(useNA="ifany", 1:5))
        expect_identical(table(useNA="ifany", c(NA, 1:5)),
            base::table(useNA="ifany", c(NA, 1:5)))
    })

    test_that("unsupported table methods", {
        expect_error(table(gen, 1:5),
            "Cannot currently tabulate Crunch variables with non-Crunch vectors")
        expect_error(table(1:5, gen),
            "Cannot currently tabulate Crunch variables with non-Crunch vectors")
        expect_error(table(), "nothing to tabulate")
    })

    test_that("unsupported aggregation methods", {
        expect_error(mean(ds$textVar),
            paste(dQuote("mean"), "is not defined for TextVariable"))
        expect_error(sd(ds$textVar),
            paste(dQuote("sd"), "is not defined for TextVariable"))
        expect_error(median(ds$textVar),
            paste(dQuote("median"), "is not defined for TextVariable"))
        expect_error(min(ds$textVar),
            paste(dQuote("min"), "is not defined for TextVariable"))
        expect_error(max(ds$textVar),
            paste(dQuote("max"), "is not defined for TextVariable"))
    })

    test_that("supported summary methods for numeric", {
        expect_equivalent(min(ds$birthyr), 1920)
        expect_equivalent(max(ds$birthyr), 1995)
        expect_equivalent(mean(ds$birthyr), 1964.951)
        expect_equivalent(median(ds$birthyr), 1968)
        expect_identical(round(sd(ds$birthyr), 2), 15.17)
    })
})

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("can fetch variable summaries", {
                summ <- getSummary(ds$v1)
                expect_true(is.list(summ))
                expect_equivalent(summ$mean, mean(df$v1, na.rm=TRUE))
                expect_equivalent(summ$stddev, sd(df$v1, na.rm=TRUE))
            })
            test_that("method dispatch", {
                expect_identical(mean(ds$v1), mean(df$v1))
                expect_equivalent(mean(ds$v1, na.rm=TRUE),
                    mean(df$v1, na.rm=TRUE))
                expect_identical(sd(ds$v1), sd(ds$v1))
                expect_equivalent(sd(ds$v1, na.rm=TRUE),
                    sd(ds$v1, na.rm=TRUE))
                expect_identical(median(ds$v1), median(ds$v1))
                expect_identical(median(ds$v1, na.rm=TRUE),
                    median(ds$v1, na.rm=TRUE))
            })
            test_that("table", {
                expect_equivalent(table(ds$v4), table(df$v4))
                expect_equivalent(table(ds$v4, ds$v3), table(df$v4, df$v3))
            })
            test_that("table works with CrunchExpr", {
                expect_equivalent(table(ds$v4[ds$v3 < 10]),
                    table(df$v4[df$v3 < 10]))
            })
            test_that("table throws error if not equally filtered", {
                expect_error(table(ds$v4, ds$v2[ds$v3 < 10]),
                    "Filter expressions in variables must be identical")
            })
            test_that("summary", {
                expect_equivalent(round(unclass(summary(ds$v1)), 2),
                    round(unclass(summary(df$v1)), 2))
                expect_equivalent(as.numeric(summary(ds$v4)),
                    summary(df$v4))
            })
        })
    })
}
