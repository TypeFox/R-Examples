context("Versions")

with_mock_HTTP({
    ds <- loadDataset("test ds")
    test_that("Version catalog exists", {
        expect_true(inherits(versions(ds), "VersionCatalog"))
    })

    test_that("Version catalog attributes (and sorting)", {
        expect_identical(length(versions(ds)), 2L)
        expect_identical(names(versions(ds)),
            c("another version", "initial load"))
        expect_identical(descriptions(versions(ds)),
            c("another version", "initial load"))
        expect_identical(strftime(timestamps(versions(ds)), "%d %b %Y"),
            c("15 Feb 2015", "12 Feb 2015"))
    })

    test_that("Version catalog print method", {
        expect_identical(formatVersionCatalog(versions(ds),
                            from=strptime("2015-02-17", "%Y-%m-%d")),
            data.frame(Name=c("another version", "initial load"),
                Timestamp=c("1 day ago", "4 days ago"),
                stringsAsFactors=FALSE))
    })

    test_that("saveVersion makes the right request", {
        expect_error(saveVersion(ds, "Today"),
            'POST /api/datasets/dataset1/savepoints.json \\{"description":"Today"\\}')
    })

    test_that("restoreVersion makes the right request", {
        expect_error(restoreVersion(ds, "initial load"),
            'POST /api/datasets/dataset1/savepoints/v2/revert/')
        expect_error(restoreVersion(ds, 2),
            'POST /api/datasets/dataset1/savepoints/v2/revert/')
        expect_error(restoreVersion(ds, "not a version"),
            paste0(dQuote("not a version"),
            " does not match any available versions"))
    })
})


if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("Dataset imported correctly", {
                validImport(ds)
            })

            ## Release and re-lease
            ds <- releaseAndReload(ds)

            test_that("There is an initial version", {
                expect_identical(names(versions(ds)), "initial import")
                expect_identical(length(versions(ds)), 1L)
            })

            ## Make changes:
            # 1. Edit variable metadata
            names(categories(ds$v4))[1:2] <- c("d", "e")
            name(ds$v2) <- "Variable Two"
            description(ds$v3) <- "The third variable in the dataset"

            # 2. Edit dataset metadata
            description(ds) <- "A dataset for testing"

            # 3. Reorder variables
            ordering(ds) <- VariableOrder(VariableGroup("Even", ds[c(2,4,6)]),
                VariableGroup("Odd", ds[c(1,3,5)]))

            # 4. Derive variable
            ds$v7 <- ds$v3 - 6

            # 5. Add non-derived variable
            ds$v8 <- rep(1:5, 4)

            ## Assert those things
            test_that("The edits are made", {
                expect_identical(names(na.omit(categories(ds$v4))),
                    c("d", "e"))
                expect_identical(name(ds$v2), "Variable Two")
                expect_identical(description(ds$v3),
                    "The third variable in the dataset")
                expect_identical(description(ds), "A dataset for testing")
                expect_identical(as.vector(ds$v7), df$v3 - 6)
                expect_equivalent(as.vector(ds$v8), rep(1:5, 4))
                expect_identical(aliases(variables(ds)),
                    paste0("v", c(2,4,6,1,3,5,7,8)))
            })

            ## Release and re-lease
            ds <- releaseAndReload(ds)

            ## Assert those things again
            test_that("The edits made are still there after releasing", {
                expect_identical(names(na.omit(categories(ds$v4))),
                    c("d", "e"))
                expect_identical(name(ds$v2), "Variable Two")
                expect_identical(description(ds$v3),
                    "The third variable in the dataset")
                expect_identical(description(ds), "A dataset for testing")
                expect_identical(as.vector(ds$v7), df$v3 - 6)
                expect_equivalent(as.vector(ds$v8), rep(1:5, 4))
                expect_identical(aliases(variables(ds)),
                    paste0("v", c(2,4,6,1,3,5,7,8)))
            })

            ## Save a version
            try(saveVersion(ds, "My changes"))
            test_that("There are now two versions", {
                expect_identical(length(versions(ds)), 2L)
                expect_identical(names(versions(ds))[1], "My changes")
            })

            ## Release and re-lease
            ds <- releaseAndReload(ds)

            ## Revert to the first version
            ds <- try(restoreVersion(ds, "initial import"))
            test_that("Restoring restored correctly", {
                expect_identical(length(versions(ds)), 1L)
                validImport(ds)
            })

            ## Release and re-lease
            ds <- releaseAndReload(ds)

            test_that("Added variables are really removed by rolling back", {
                ## This was user-reported: Order was reverted but derived
                ## variables persisted, and by assigning an empty order, you can
                ## recover them.
                ordering(ds) <- VariableOrder()
                expect_true(setequal(names(ds), names(df)))
            })

            test_that("And now we can add variables again that we added and reverted", {
                expect_true(is.null(ds$v7))
                ## This would error if "v7" were still lurking somewhere
                ds$v7 <- ds$v3 - 7
                expect_identical(as.vector(ds$v7), df$v3 - 7)

                expect_true(is.null(ds$v8))
                ds$v8 <- rep(6:10, 4)
                expect_equivalent(as.vector(ds$v8), rep(6:10, 4))
            })
        })
    })
}
