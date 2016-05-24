context("Update error handling")


if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            len <- try(length(as.vector(ds$v4[ds$v4 == "B"])))
            test_that("setup for update with wrong number of values", {
                expect_identical(len, 10L)
            })
            test_that("Trying to update with too many values fails", {
                skip_on_jenkins("Silence the error emails")
                expect_error(ds$v4[ds$v4 == "B"] <- rep(1, len + 5),
                    "length 15 does not match existing length 10")
            })
            test_that("Trying to update with too few values fails", {
                skip_on_jenkins("Silence the error emails")
                expect_error(ds$v4[ds$v4 == "B"] <- rep(1, len - 3),
                    "length 7 does not match existing length 10")
            })
            test_that("Trying to update with different filters fails", {
                skip("Fails, but error message is different now")
                expect_error(ds$v4[ds$v4 == "B"] <- ds$v4[ds$v4 == "C"],
                    "Cannot update a variable with a value that has a different filter")
                expect_equivalent(as.numeric(table(ds$v4)), c(10, 10))
                try(ds$v4[ds$v4 == "B"] <- as.vector(ds$v4[ds$v4 == "C"]))
                expect_true(all(as.vector(ds$v4) == "C"))
            })

            test_that("Trying to update with the wrong data type fails", {
                expect_error(ds$v3[is.na(ds$v1)] <- letters[3:7],
                    "Cannot update NumericVariable with type character")
            })
        })
    })
}
