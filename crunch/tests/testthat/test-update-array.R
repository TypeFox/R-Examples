context("Update array variables")


if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(mrdf), {
            ds <- mrdf.setup(ds)
            test_that("Subvariable values before trying to update", {
                expect_equivalent(as.vector(ds$CA$mr_1),
                    as.factor(c("1.0", "0.0", "1.0", NA)))
                expect_equivalent(as.vector(ds$CA$mr_2),
                    as.factor(c("0.0", "0.0", "1.0", NA)))
                expect_equivalent(as.vector(ds$CA$mr_3),
                    as.factor(c("0.0", "0.0", "1.0", NA)))
            })
            try(ds$CA[ds$v4 == "B"] <- c("1.0"))
            test_that("Can update array subvariables", {
                expect_equivalent(as.vector(ds$CA$mr_1),
                    as.factor(c("1.0", "0.0", "1.0", NA)))
                expect_equivalent(as.vector(ds$CA$mr_2),
                    as.factor(c("1.0", "0.0", "1.0", NA)))
                expect_equivalent(as.vector(ds$CA$mr_3),
                    as.factor(c("1.0", "0.0", "1.0", NA)))
            })
        })

        with(test.dataset(newDatasetFromFixture("apidocs")), {
            test_that("Can update an individual subvariable conditionally", {
                expect_equivalent(as.vector(ds$allpets$Cat, mode="id")[1:5],
                    c(1, 9, 1, 1, 9))
                ds$allpets$Cat[2] <- 1
                expect_equivalent(as.vector(ds$allpets$Cat, mode="id")[1:5],
                    c(1, 1, 1, 1, 9))
            })
        })
    })
}
