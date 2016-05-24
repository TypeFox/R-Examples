context("Weights")

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("Can set weight variable", {
                expect_identical(weight(ds), NULL)
                weight(ds) <- ds$v3
                expect_equivalent(weight(ds), ds$v3)
                ds <- refresh(ds)
                expect_equivalent(weight(ds), ds$v3)
                weight(ds) <- NULL
                expect_identical(weight(ds), NULL)
            })
            test_that("Errors are properly handled when setting weight", {
                expect_error(weight(ds) <- "a",
                    "Weight must be a Variable or NULL")
                ## test error handling when trying to set non-numeric
            })

            test_that("If weight is set, computations are weighted", {
                expect_equivalent(table(ds$v4),
                    structure(c(B=10, C=10), class="table"))
                weight(ds) <- ds$v3
                expect_equivalent(table(ds$v4),
                    structure(c(B=sum(seq(8, 26, 2)), C=sum(seq(9, 27, 2))),
                    class="table"))
            })

            test_that("If weight is set, dim() is still unweighted", {
                weight(ds) <- NULL
                expect_identical(nrow(ds), 20L)
                weight(ds) <- ds$v3
                expect_identical(nrow(ds), 20L)
                ds <- refresh(ds)
                expect_identical(nrow(ds), 20L)
            })

        })
        with(test.dataset(df), {
            test_that("We have a clean dataset", {
                expect_identical(weight(ds), NULL)
            })
            test_that("Reverting to old version rolls back weight variables", {
                ds <- saveVersion(ds, "Before w")
                ds$w <- 1:20
                weight(ds) <- ds$w
                expect_equivalent(weight(ds), ds$w)
                ds <- restoreVersion(ds, "Before w")
                expect_identical(weight(ds), NULL)
            })
            test_that("And I can add new weights because weight_variables is valid", {
                ds <- refresh(ds)
                ds$w2 <- 2:21
                weight(ds) <- ds$w2
                expect_equivalent(weight(ds), ds$w2)
            })
        })

        with(test.dataset(df), {
            test_that("I can delete my weight variable and add a new one", {
                ds$w <- 1:20
                weight(ds) <- ds$w
                expect_equivalent(weight(ds), ds$w)
                expect_true(is.Numeric(ds$w))
                expect_equivalent(as.array(crtabs(~ v4, data=ds)),
                    array(c(100, 110), dim=2L, dimnames=list(v4=c("B", "C"))))
                ## Delete that variable. Confirm that it is gone and we are
                ## unweighted
                with(consent(), ds$w <- NULL)
                expect_identical(weight(ds), NULL)
                expect_identical(ds$w, NULL)
                expect_equivalent(as.array(crtabs(~ v4, data=ds)),
                    array(c(10, 10), dim=2L, dimnames=list(v4=c("B", "C"))))
                ## Now add another weight and repeat. Confirm that we can
                ## and that calculations are weighted
                ds$w <- 20:1
                weight(ds) <- ds$w
                expect_equivalent(weight(ds), ds$w)
                expect_true(is.Numeric(ds$w))
                expect_equivalent(as.array(crtabs(~ v4, data=ds)),
                    array(c(110, 100), dim=2L, dimnames=list(v4=c("B", "C"))))
                ## Now force the dataset to drop on the server and reload it
                ## to confirm that our changes were persisted correctly
                ds <- releaseAndReload(ds)
                expect_equivalent(weight(ds), ds$w)
                expect_true(is.Numeric(ds$w))
                expect_equivalent(as.array(crtabs(~ v4, data=ds)),
                    array(c(110, 100), dim=2L, dimnames=list(v4=c("B", "C"))))
            })
        })
    })
}
