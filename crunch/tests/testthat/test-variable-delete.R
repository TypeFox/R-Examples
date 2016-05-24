context("Deleting variables")

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("deleteVariable(s)", {
                v1 <- ds$v1
                expect_true(all(c("v1", "v4") %in% names(ds)))
                ds <- deleteVariable(ds, c("v1", "v4"))
                expect_true(!any(c("v1", "v4") %in% names(ds)))
                expect_error(refresh(v1),
                    "Variable not found. It may have been deleted.")
            })
        })

        with(test.dataset(df), {
            test_that("deleteVariables with consent", {
                with(temp.option(crunch.require.confirmation=TRUE), {
                    expect_true("v3" %in% names(ds))
                    expect_error(ds <- deleteVariable(ds, "v3"),
                        "Must confirm deleting variable")
                    with(consent(), {
                        ds <- deleteVariable(ds, "v3")
                        expect_true(is.null(ds$v3))
                        expect_true(is.null(refresh(ds)$v3))
                    })
                })
            })
        })

        with(test.dataset(df), {
            test_that("Delete variable by assigning NULL", {
                with(temp.option(crunch.require.confirmation=TRUE), {
                    expect_true("v3" %in% names(ds))
                    expect_identical(names(ds)[2], "v2")
                    expect_error(ds$v3 <- NULL,
                        "Must confirm deleting variable")
                    expect_error(ds[[2]] <- NULL,
                        "Must confirm deleting variable")
                    with(consent(), {
                        ds$v3 <- NULL
                        expect_true(is.null(ds$v3))
                        expect_true(is.null(refresh(ds)$v3))
                        ds[[2]] <- NULL
                        expect_true(is.null(ds$v2))
                        expect_true(is.null(refresh(ds)$v2))
                    })
                })
            })

            test_that("Assigning NULL doesn't ask you about deleting 0 variables", {
                expect_message(ds$NOTAVARIABLE <- df$NOTAVARIABLE,
                    paste(dQuote("NOTAVARIABLE"), "is not a variable; nothing to delete by assigning NULL"))
            })
        })

        with(test.dataset(newDatasetFromFixture("apidocs")), {
            test_that("Array variables are fully deleted", {
                expect_true("petloc" %in% names(ds))
                expect_false("petloc_home" %in% names(ds))
                ds$petloc <- NULL
                expect_false("petloc" %in% names(ds))
                expect_false("petloc_home" %in% names(ds))
            })

            test_that("Deleting array/MR subvariables", {
                expect_true("allpets" %in% names(ds))
                expect_true(is.MR(ds$allpets))
                expect_false("allpets_2" %in% names(ds))
                expect_true("allpets_2" %in% aliases(subvariables(ds$allpets)))
                # ds$allpets <- deleteSubvariable(ds$allpets, "Dog") ## Cannot overwrite one Variable with another
                deleteSubvariable(ds$allpets, "Dog")
                ds <- refresh(ds)
                expect_true("allpets" %in% names(ds))
                expect_true(is.MR(ds$allpets))
                expect_false("allpets_2" %in% names(ds))
                expect_false("allpets_2" %in% aliases(subvariables(ds$allpets)))
                expect_identical(names(subvariables(ds$allpets)),
                    c("Cat", "Bird"))
            })
        })

        with(test.dataset(newDatasetFromFixture("apidocs")), {
            ## Attempt to reproduce @persephonet's error
            cats <- categories(ds$allpets)
            cats[[3]]$selected <- TRUE
            categories(ds$allpets) <- cats
            test_that("Delete MR subvar with multiple 'selected' attributes", {
                expect_identical(names(Filter(is.selected, categories(ds$allpets))),
                    c("selected", "not asked")) ## There are two selected
                deleteSubvariable(ds$allpets, "Dog")
                ds <- refresh(ds)
                expect_true("allpets" %in% names(ds))
                expect_true(is.MR(ds$allpets))
                expect_false("allpets_2" %in% names(ds))
                expect_false("allpets_2" %in% aliases(subvariables(ds$allpets)))
                expect_identical(names(subvariables(ds$allpets)),
                    c("Cat", "Bird"))
            })
        })

        with(test.dataset(newDatasetFromFixture("apidocs")), {
            test_that("delete method on variable", {
                expect_true("q1" %in% names(ds))
                d <- try(delete(ds$q1))
                ds <- refresh(ds)
                expect_false("q1" %in% names(ds))
            })
            test_that("delete deletes all subvariables too", {
                expect_true("petloc" %in% names(ds))
                d <- try(delete(ds$petloc))
                ds <- refresh(ds)
                expect_false("petloc" %in% names(ds))
                expect_false("petloc_home" %in% names(ds))
            })
        })
    })
}
