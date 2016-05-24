context("Variable types")

test_that("Validation on type setting", {
    expect_error(castVariable(, "foo"),
        paste(dQuote("foo"),
        "is not a Crunch variable type that can be assigned."))
})

with_mock_HTTP({
    test.ds <- loadDataset("test ds")

    test_that("Variable type method", {
        expect_identical(type(test.ds[["birthyr"]]), "numeric")
        expect_identical(type(test.ds$gender), "categorical")
    })
})

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("type casting and 'as'", {
                testvar <- ds$v1
                type(testvar) <- "text"
                expect_true(is.Text(testvar))
                expect_true(is.Numeric(castVariable(testvar, "numeric")))
                expect_true(is.Text(castVariable(testvar, "text")))

                type(testvar) <- "numeric"
                expect_true(is.Numeric(testvar))
                expect_true(is.Numeric(ds$v1)) ## same remote object
            })
        })
    })
}
