context("Update a dataset")

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("Can update numeric variable with values", {
                ds$v3 <- 9:28
                test <- as.vector(ds$v3) - df$v3
                expect_true(all(test == 1))
            })

            ds$v3 <- 1
            test_that("Value recycling on insert is consistent with R", {
                expect_true(all(as.vector(ds$v3) == 1))
            })

            ds$v3[1:10] <- 2
            test_that("Update numeric with R numeric filter and values", {
                expect_equivalent(mean(ds$v3), 1.5)
            })
            ds$v3[ds$v3 == 1] <- 3
            test_that("Update numeric with LogicalExpression filter", {
                expect_equivalent(mean(ds$v3), 2.5)
            })
            ds[ds$v3 == 2, "v3"] <- 4
            test_that("Update with LogicalExpression within dataset", {
                expect_equivalent(mean(ds$v3), 3.5)
            })
            ds$v3 <- c(rep(5, 10), rep(7, 10))
            test_that("Just update the values", {
                expect_equivalent(mean(ds$v3), 6)
            })

            test_that("Can update numeric variable with expresssion", {
                ds$v3 <- ds$v3 + 2
                expect_equivalent(as.vector(ds$v3), c(rep(7, 10), rep(9, 10)))
            })

            test_that("Can filter on is.na", {
                ds$v3[is.na(ds$v2)] <- 0
                expect_equivalent(as.vector(ds$v3),
                    c(rep(7, 10), rep(9, 5), rep(0, 5)))
            })

            test_that("Can update text", {
                ds$v2[is.na(ds$v1)] <- "z"
                expect_identical(as.vector(ds$v2)[1:8],
                    c(rep("z", 5), "f", "g", "h"))
                ds[ds$v2 %in% "z", "v2"] <- "y"
                expect_identical(as.vector(ds$v2)[1:8],
                    c(rep("y", 5), "f", "g", "h"))
            })

            test_that("Can update datetime", {
                newvals <- as.Date(0:12, origin="1985-10-26")
                ds$v5[ds$v5 >= as.Date("1955-11-12")] <- newvals
                expect_identical(max(ds$v5), as.Date("1985-11-07"))
            })

            date.before <- rep(c("2014-04-15", "2014-08-15"), 2)
            date.after <- c("2014-04-15", "2014-09-15", "2014-04-15",
                "2014-09-15")
            date.df <- data.frame(wave=as.Date(date.before))
            with(test.dataset(date.df, "date.ds"), {
                test_that("Another datetime update", {
                    expect_identical(as.vector(date.ds$wave),
                        as.Date(date.before))
                    date.ds$wave[date.ds$wave == as.Date("2014-08-15")] <- as.Date("2014-09-15")
                    expect_identical(as.vector(date.ds$wave),
                        as.Date(date.after))
                })
            })

            ## Categorical
            ds$v4[is.na(ds$v2)] <- "B"
            test_that("Can update categorical variables with character", {
                expect_equivalent(table(ds$v4)["B"], c(B=13L))
            })
            ds$v4[is.na(ds$v2)] <- factor("C")
            test_that("Can update categorical with factor", {
                expect_equivalent(table(ds$v4)["C"], c(C=12L))
            })
            ds$v4[is.na(ds$v2)] <- c(2,1,2,1,2)
            test_that("Can update categorical with numeric (ids)", {
                expect_equivalent(table(ds$v4), table(df$v4))
            })
            test_that("Validation on categorical update", {
                expect_error(ds$v4[is.na(ds$v2)] <- as.factor(LETTERS[1:5]),
                    "Input values A, D, and E are not present in the category names of variable")
            })
        })
    })
}
